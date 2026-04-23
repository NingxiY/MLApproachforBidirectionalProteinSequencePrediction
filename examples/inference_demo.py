"""
inference_demo.py — Koina fragment intensity prediction wrapper

Supports two fragmentation modes:
  CID -> Prosit_2020_intensity_CID      (2 inputs: sequence + charge)
  HCD -> AlphaPeptDeep_ms2_generic      (4 inputs: sequence + charge + CE + instrument)

Usage:
    from inference_demo import predict_single, predict_batch

    result = predict_single("PEPTIDE", charge=2, fragmentation="CID")
    results = predict_batch([{"sequence": "PEPTIDE", "charge": 2}], fragmentation="HCD")
"""

import json
import requests

# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------
DEFAULT_ENDPOINT = "http://localhost:8502"
PUBLIC_ENDPOINT  = "https://koina.proteomicsdb.org"

# ---------------------------------------------------------------------------
# Model routing table
# ---------------------------------------------------------------------------
MODEL_ROUTES = {
    "CID": {
        "model":              "Prosit_2020_intensity_CID",
        "has_ce":             False,
        "has_instrument":     False,
    },
    "HCD": {
        "model":              "AlphaPeptDeep_ms2_generic",
        "has_ce":             True,
        "has_instrument":     True,
        "default_ce":         30,
        "default_instrument": "Lumos",
    },
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def predict_single(
    sequence: str,
    charge: int,
    fragmentation: str,
    endpoint: str = DEFAULT_ENDPOINT,
    ce: int = 30,
) -> dict:
    """
    Predict fragment intensities for a single peptide.

    Args:
        sequence:      Amino acid sequence (uppercase, standard residues)
        charge:        Precursor charge state
        fragmentation: "CID" or "HCD"
        endpoint:      Koina server base URL
        ce:            Collision energy (HCD only; ignored for CID)

    Returns:
        {
          "sequence":      str,
          "charge":        int,
          "fragmentation": str,
          "ions": [
              {"annotation": "y3+1", "mz": float, "intensity": float},
              ...                        # z=1 b/y ions only, sorted by type then position
          ]
        }

    Raises:
        ValueError:          Unknown fragmentation type
        RuntimeError:        annotation field missing in API response
        requests.HTTPError:  Non-2xx HTTP response from Koina
    """
    results = predict_batch(
        [{"sequence": sequence, "charge": charge}],
        fragmentation=fragmentation,
        endpoint=endpoint,
        ce=ce,
    )
    return results[0]


def predict_batch(
    peptides: list,
    fragmentation: str,
    endpoint: str = DEFAULT_ENDPOINT,
    ce: int = 30,
) -> list:
    """
    Predict fragment intensities for a batch of peptides (single API call).

    Args:
        peptides:      [{"sequence": str, "charge": int}, ...]
        fragmentation: "CID" or "HCD"
        endpoint:      Koina server base URL
        ce:            Collision energy for all peptides (HCD only)

    Returns:
        List of result dicts (same order and length as input):
        [{"sequence": str, "charge": int, "fragmentation": str, "ions": [...]}, ...]

    Raises:
        ValueError:          Unknown fragmentation type or empty peptide list
        RuntimeError:        annotation field missing in API response
        requests.HTTPError:  Non-2xx HTTP response from Koina
    """
    if fragmentation not in MODEL_ROUTES:
        raise ValueError(
            f"Unknown fragmentation '{fragmentation}'. Must be one of: {list(MODEL_ROUTES)}"
        )
    if not peptides:
        raise ValueError("peptides list is empty")

    route     = MODEL_ROUTES[fragmentation]
    sequences = [p["sequence"] for p in peptides]
    charges   = [p["charge"]   for p in peptides]
    payload   = _build_payload(sequences, charges, fragmentation, ce)
    url       = f"{endpoint.rstrip('/')}/v2/models/{route['model']}/infer"

    response = requests.post(
        url,
        headers={"Content-Type": "application/json"},
        data=json.dumps(payload),
        timeout=60,
    )
    response.raise_for_status()

    ion_lists = _parse_response(response.json(), n=len(peptides))

    return [
        {
            "sequence":      peptide["sequence"],
            "charge":        peptide["charge"],
            "fragmentation": fragmentation,
            "ions":          ions,
        }
        for peptide, ions in zip(peptides, ion_lists)
    ]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _build_payload(
    sequences: list,
    charges: list,
    fragmentation: str,
    ce: int,
) -> dict:
    """
    Build a Triton HTTP inference request payload.

    Data format: flat list with shape descriptor [n, 1].
    Triton maps flat data + shape to the correct tensor layout.
    This matches the format used in all existing eval_koina_*.py scripts.
    """
    route = MODEL_ROUTES[fragmentation]
    n = len(sequences)

    inputs = [
        {
            "name":     "peptide_sequences",
            "shape":    [n, 1],
            "datatype": "BYTES",
            "data":     sequences,          # flat list: ["SEQ1", "SEQ2", ...]
        },
        {
            "name":     "precursor_charges",
            "shape":    [n, 1],
            "datatype": "INT32",
            "data":     charges,            # flat list: [2, 3, ...]
        },
    ]

    if route["has_ce"]:
        inputs.append({
            "name":     "collision_energies",
            "shape":    [n, 1],
            "datatype": "FP32",
            "data":     [float(ce)] * n,
        })

    if route["has_instrument"]:
        inputs.append({
            "name":     "instrument_types",
            "shape":    [n, 1],
            "datatype": "BYTES",
            "data":     [route["default_instrument"]] * n,
        })

    return {"id": "koina_request", "inputs": inputs}


def _parse_response(response_json: dict, n: int) -> list:
    """
    Extract z=1 b/y ions from a Triton inference response.

    Returns a list of length n. Each element is a list of ion dicts:
        [{"annotation": "y3+1", "mz": float, "intensity": float}, ...]
    Ions are sorted by type (b before y) then by position number.

    Notes on model differences handled here:
      - Prosit CID: fixed 174-element output; invalid ions have intensity = -1.0
      - APD HCD:    dynamic output (only valid ions); no -1.0 sentinels; FP64 intensities

    Raises:
        RuntimeError: if 'annotation' output is absent from the response.
    """
    outputs = {o["name"]: o for o in response_json.get("outputs", [])}

    if "annotation" not in outputs:
        raise RuntimeError(
            "Koina response is missing the 'annotation' output field. "
            f"Present outputs: {list(outputs)}. "
            "Cannot determine ion identity without annotations."
        )

    raw_annotations = outputs["annotation"]["data"]
    raw_intensities = outputs["intensities"]["data"]
    raw_mz          = outputs["mz"]["data"]

    total = len(raw_annotations)
    if total % n != 0:
        raise RuntimeError(
            f"Output length {total} is not evenly divisible by batch size {n}."
        )
    stride = total // n

    results = []
    for i in range(n):
        s = i * stride
        e = s + stride

        ions = []
        for ann, intensity, mz in zip(
            raw_annotations[s:e],
            raw_intensities[s:e],
            raw_mz[s:e],
        ):
            if float(intensity) < 0:        # Prosit sentinel for inapplicable ions
                continue
            if not ann.endswith("+1"):      # keep z=1 only
                continue
            if not (ann.startswith("b") or ann.startswith("y")):
                continue
            ions.append({
                "annotation": ann,
                "mz":         float(mz),
                "intensity":  float(intensity),
            })

        ions.sort(key=lambda x: (x["annotation"][0], _ion_position(x["annotation"])))
        results.append(ions)

    return results


def _ion_position(annotation: str) -> int:
    """Parse numeric fragment position from annotation, e.g. 'y12+1' -> 12."""
    try:
        return int(annotation[1:annotation.index("+")])
    except (ValueError, IndexError):
        return 0
