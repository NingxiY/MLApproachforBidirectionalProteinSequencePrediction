"""
app.py — ProteinBridge backend (Phase 1, Step 1.3)

Serves index.html at GET / and accepts file uploads at POST /predict.
seq2spec: real inference via parse_fasta() + predict_batch() (Koina).
spec2seq: Casanovo 5.1.2 subprocess via casanovo_env_v5.

Run:
    uvicorn app:app --reload --port 8000
"""

import os
import shutil
import subprocess
import tempfile

from fastapi import FastAPI, File, Form, UploadFile
from fastapi.responses import FileResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware

from parsers import parse_fasta, save_temp_file, parse_mztab
from examples.inference_demo import predict_batch, PUBLIC_ENDPOINT

# ---------------------------------------------------------------------------
# Casanovo paths — resolved relative to this file so the server can be
# launched from any working directory.
# ---------------------------------------------------------------------------
_HERE      = os.path.dirname(os.path.abspath(__file__))
# Use the wrapper script so that DYLD_FRAMEWORK_PATH / DYLD_LIBRARY_PATH are
# set inside the spawned process itself.  When uvicorn runs under the macOS
# system Python, those env vars are stripped from child processes by dyld
# unless the child sets them itself.
_CASANOVO  = os.path.join(_HERE, "casanovo_env_v5", "bin", "casanovo_run.sh")
_CKPT      = os.path.join(_HERE, "casanovo", "casanovo.ckpt")
_CONFIG    = os.path.join(_HERE, "casanovo", "casanovo.yaml")

app = FastAPI(title="ProteinBridge API")

# Allow the browser to call /predict when index.html is opened directly
# from the filesystem (file://) during development.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["POST"],
    allow_headers=["*"],
)


@app.get("/")
def serve_index():
    """Serve the frontend."""
    return FileResponse("index.html")


@app.post("/predict")
async def predict(
    file: UploadFile = File(...),
    mode: str = Form(...),
):
    """
    Accept a file upload and a mode string, run inference, return results.

    Args:
        file: uploaded FASTA (.fasta/.fa/.txt) for seq2spec,
              or .mgf / .mzML spectrum file for spec2seq
        mode: "seq2spec" or "spec2seq"

    Returns:
        seq2spec: {"mode": "seq2spec", "results": [{sequence, charge,
                   fragmentation, ions: [{annotation, mz, intensity}]}]}
        spec2seq: {"mode": "spec2seq", "results": [{sequence, score, aa_scores}]}
    """
    if mode not in ("seq2spec", "spec2seq"):
        return JSONResponse(
            status_code=400,
            content={"error": f"Unknown mode '{mode}'. Must be 'seq2spec' or 'spec2seq'."},
        )

    contents = await file.read()

    if mode == "seq2spec":
        try:
            peptides = parse_fasta(contents)
        except ValueError as exc:
            return JSONResponse(status_code=400, content={"error": str(exc)})

        results = predict_batch(peptides, fragmentation="CID", endpoint=PUBLIC_ENDPOINT)
        return {"mode": "seq2spec", "results": results}

    # ── spec2seq — Casanovo ──────────────────────────────────────────────────

    # Validate file extension.
    filename = (file.filename or "").lower()
    if not (filename.endswith(".mgf") or filename.endswith(".mzml")):
        return JSONResponse(
            status_code=400,
            content={"error": "spec2seq requires an .mgf or .mzML file."},
        )

    if not contents:
        return JSONResponse(status_code=400, content={"error": "Uploaded file is empty."})

    for label, path in [
        ("Casanovo executable", _CASANOVO),
        ("Model checkpoint",    _CKPT),
        ("Config file",         _CONFIG),
    ]:
        if not os.path.exists(path):
            return JSONResponse(
                status_code=500,
                content={"error": f"{label} not found at {path}."},
            )

    # Write the uploaded spectrum to a temp file.
    suffix = ".mzML" if filename.endswith(".mzml") else ".mgf"
    input_path = save_temp_file(contents, suffix=suffix)
    work_dir   = None

    try:
        # Casanovo 5.x always appends ".mztab" to whatever -o name is given,
        # and on Python 3.10 it cannot handle absolute output paths (pathlib
        # glob limitation).  Run casanovo inside a dedicated temp directory and
        # use a bare basename so the output lands at {work_dir}/result.mztab.
        work_dir = tempfile.mkdtemp(prefix="pb_casanovo_")
        output_base = "result"
        output_mztab = os.path.join(work_dir, f"{output_base}.mztab")

        cmd = [
            _CASANOVO, "sequence",
            "-m", _CKPT,
            "-c", _CONFIG,
            "-o", output_base,   # bare name — Casanovo appends .mztab
            input_path,
        ]

        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
            cwd=work_dir,   # run inside the temp dir so -o is a relative path
        )

        if proc.returncode != 0:
            return JSONResponse(
                status_code=500,
                content={
                    "error": "Casanovo subprocess failed.",
                    "stderr": proc.stderr[-2000:],
                },
            )

        if not os.path.exists(output_mztab) or os.path.getsize(output_mztab) == 0:
            return JSONResponse(
                status_code=500,
                content={"error": "Casanovo produced no output. Check that the spectrum file is valid."},
            )

        try:
            results = parse_mztab(output_mztab)
        except (FileNotFoundError, ValueError) as exc:
            return JSONResponse(
                status_code=500,
                content={"error": f"Could not parse Casanovo output: {exc}"},
            )

    finally:
        try:
            os.unlink(input_path)
        except OSError:
            pass
        if work_dir:
            shutil.rmtree(work_dir, ignore_errors=True)

    return {"mode": "spec2seq", "results": results}
