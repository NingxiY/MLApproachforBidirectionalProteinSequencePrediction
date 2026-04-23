"""
parsers.py — Input file parsers for ProteinBridge

parse_fasta:    FASTA bytes  -> peptide dicts for inference_demo.predict_batch()
save_temp_file: uploaded bytes -> temp file path for subprocess tools (Casanovo)
parse_mztab:    Casanovo .mztab path -> list of sequence result dicts
"""

import os
import tempfile


def parse_fasta(file_bytes: bytes) -> list[dict]:
    """
    Parse FASTA-formatted bytes into a list of peptide dicts.

    Args:
        file_bytes: raw bytes of a .fasta / .fa / .txt file

    Returns:
        [{"sequence": str, "charge": int}, ...]
        charge is fixed at 2 (standard default for demo purposes).
        Sequences are uppercased and have all whitespace removed.

    Raises:
        ValueError: if the input contains no valid (non-empty) sequences
    """
    text = file_bytes.decode("utf-8", errors="replace")

    peptides = []
    current_lines = []

    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            # Save the sequence accumulated so far before starting a new record.
            _flush(current_lines, peptides)
            current_lines = []
        elif line:
            current_lines.append(line)

    # Flush the final record.
    _flush(current_lines, peptides)

    if not peptides:
        raise ValueError(
            "No valid sequences found in the uploaded file. "
            "Make sure it is a FASTA file with at least one non-empty sequence."
        )

    return peptides


def _flush(lines: list, peptides: list) -> None:
    """Join accumulated sequence lines and append to peptides if non-empty."""
    sequence = "".join(lines).upper().replace(" ", "")
    if sequence:
        peptides.append({"sequence": sequence, "charge": 2})


# ---------------------------------------------------------------------------
# Temp file helper
# ---------------------------------------------------------------------------

def save_temp_file(file_bytes: bytes, suffix: str) -> str:
    """
    Write uploaded bytes to a named temp file and return its path.

    The file is NOT deleted automatically. The caller must call
    os.unlink(path) after the subprocess has finished with it.

    Args:
        file_bytes: raw bytes from the uploaded file
        suffix:     file extension including the dot, e.g. ".mgf" or ".mzML"
                    Casanovo inspects the extension to choose its parser.

    Returns:
        Absolute path to the temp file on disk.
    """
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    try:
        tmp.write(file_bytes)
    finally:
        tmp.close()
    return tmp.name


# ---------------------------------------------------------------------------
# Casanovo mzTab output parser
# ---------------------------------------------------------------------------

# Column names used by Casanovo — taken directly from the PSH header line in
# casanovo/casanovo_20260329194013.mztab (the only real output in this repo).
_COL_SEQUENCE  = "sequence"
_COL_SCORE     = "search_engine_score[1]"
_COL_AA_SCORES = "opt_ms_run[1]_aa_scores"


def parse_mztab(path: str) -> list[dict]:
    """
    Parse a Casanovo mzTab output file into a list of result dicts.

    Only PSM rows are read. MTD (metadata) and PSH (header) rows are used
    for setup only. Column positions are resolved dynamically from the PSH
    header so the parser does not break if Casanovo reorders columns.

    Args:
        path: absolute or relative path to the .mztab file

    Returns:
        [
          {
            "sequence":  str,          # predicted peptide sequence
            "score":     float,        # Casanovo search_engine_score[1]
            "aa_scores": list[float],  # per-residue confidence scores
          },
          ...
        ]

    Raises:
        FileNotFoundError: if path does not exist
        ValueError: if no PSH header is found, required columns are missing,
                    or no PSM rows are present
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"mzTab file not found: {path}")

    col: dict[str, int] = {}   # column name -> index within a PSM/PSH row
    results = []

    with open(path, encoding="utf-8") as fh:
        for lineno, raw_line in enumerate(fh, start=1):
            line = raw_line.rstrip("\n")
            if not line:
                continue

            row_type, _, rest = line.partition("\t")

            if row_type == "PSH":
                # Build the column index map from the header.
                # PSH row: PSH <tab> col0 <tab> col1 ...
                headers = rest.split("\t")
                col = {name: idx for idx, name in enumerate(headers)}
                _require_columns(col, [_COL_SEQUENCE, _COL_SCORE, _COL_AA_SCORES], path)
                continue

            if row_type != "PSM":
                continue  # skip MTD and any other row types

            if not col:
                raise ValueError(
                    f"{path}: PSM row at line {lineno} appeared before the PSH header."
                )

            fields = rest.split("\t")

            try:
                sequence  = fields[col[_COL_SEQUENCE]]
                score     = float(fields[col[_COL_SCORE]])
                aa_scores = [float(x) for x in fields[col[_COL_AA_SCORES]].split(",")]
            except (IndexError, ValueError) as exc:
                raise ValueError(
                    f"{path}: could not parse PSM row at line {lineno}: {exc}"
                ) from exc

            if not sequence:
                continue  # skip blank sequence rows

            results.append({
                "sequence":  sequence,
                "score":     score,
                "aa_scores": aa_scores,
            })

    if not results:
        raise ValueError(f"{path}: no PSM rows found — file may be empty or malformed.")

    return results


def _require_columns(col: dict, required: list, path: str) -> None:
    """Raise ValueError if any required column name is absent from col."""
    missing = [c for c in required if c not in col]
    if missing:
        raise ValueError(
            f"{path}: PSH header is missing required columns: {missing}. "
            f"Present columns: {list(col)}"
        )
