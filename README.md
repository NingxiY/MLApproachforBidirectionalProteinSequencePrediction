# ProteinBridge

A web app for bidirectional peptide prediction. You can upload a protein sequence file and get a predicted mass spectrum back, or upload a mass spectrum file and get predicted peptide sequences back. Both predictions use real machine learning models — no hardcoded results.

- **Sequence to Spectrum** uses [Koina](https://koina.proteomicsdb.org), a public API that runs the Prosit model to predict MS/MS fragment ion intensities.
- **Spectrum to Sequence** uses [Casanovo](https://github.com/Noble-Lab/casanovo), a transformer model that performs de novo peptide sequencing directly from raw spectra.

---

## Features

| Mode | You upload | You get back |
|------|-----------|--------------|
| Sequence to Spectrum | A FASTA file (`.fasta`, `.fa`, or `.txt`) | A predicted MS/MS fragment ion spectrum (bar chart + peak stats) |
| Spectrum to Sequence | A spectrum file (`.mgf` or `.mzML`) | Predicted peptide sequences with per-residue confidence scores |

---

## Project Structure

```
app.py              Backend server. Handles file uploads, calls Koina or Casanovo,
                    returns JSON results.

index.html          Frontend. Single HTML file, no build step needed.

parsers.py          File parsers: reads FASTA files and Casanovo mzTab output.

examples/
  inference_demo.py Koina API wrapper used by app.py for seq2spec predictions.

tests/
  test_parsers.py   Unit tests for the FASTA parser (runs without internet).
  test_inference.py Integration tests against the public Koina API (needs internet).
  conftest.py       Pytest setup file.

casanovo/
  casanovo.yaml     Casanovo configuration (committed).
  casanovo.ckpt     Model weights — NOT committed, ~548 MB, must be provided locally.
  predicted_spectra_noseq.mgf  Sample spectrum file for testing spec2seq.
  predicted_spectra.mgf        Same spectra with known sequences (for reference).

sample.fasta        A one-peptide FASTA file for testing seq2spec.
DEPLOYMENT.md       Notes on running a self-hosted Koina server.
```

---

## Requirements

- **Python 3.9 or later** for running the server
- **Python 3.10** separately for the Casanovo environment (see Casanovo Setup below)
- **pip** (comes with Python)
- **Internet connection** for seq2spec (calls the public Koina API at `koina.proteomicsdb.org`)
- **~2 GB free disk space** if you set up the Casanovo environment and model weights
- macOS or Linux recommended; Windows should work but is untested

---

## Setup

### 1. Clone the repository

```bash
git clone <your-repo-url>
cd researchCapstone
```

### 2. Install server dependencies

```bash
python3 -m pip install fastapi uvicorn python-multipart requests
```

### 3. Start the server

```bash
python3 -m uvicorn app:app --reload --port 8000
```

You should see output like:

```
INFO:     Uvicorn running on http://127.0.0.1:8000
```

### 4. Open the app

Open your browser and go to:

```
http://localhost:8000
```

---

## Casanovo Setup (required for Spectrum to Sequence)

The Spectrum to Sequence feature requires two things that are not included in the repository because they are too large:

**1. The model checkpoint** (`casanovo/casanovo.ckpt`, ~548 MB)

This file contains the trained Casanovo model weights. You need to place it at `casanovo/casanovo.ckpt` inside the project directory. Obtain it from your team or download from the [Casanovo releases page](https://github.com/Noble-Lab/casanovo/releases).

**2. A Python 3.10 virtual environment** (`casanovo_env_v5/`)

Casanovo 5.1.2 requires Python 3.10 and a specific set of dependencies. The app expects a virtual environment at `casanovo_env_v5/` in the project root. To create it:

```bash
# Install Python 3.10 if you do not have it
# macOS: download from https://www.python.org/downloads/

# Create the environment
python3.10 -m venv casanovo_env_v5

# Install Casanovo from GitHub (PyPI version is outdated)
DYLD_FRAMEWORK_PATH="$HOME/Library/Frameworks" \
DYLD_LIBRARY_PATH="$HOME/Library/Frameworks/Python.framework/Versions/3.10/lib" \
casanovo_env_v5/bin/pip install "git+https://github.com/Noble-Lab/casanovo.git@v5.1.2"
```

Then create a wrapper script that the server calls:

```bash
cat > casanovo_env_v5/bin/casanovo_run.sh << 'EOF'
#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
FRAMEWORKS="$HOME/Library/Frameworks"
DYLIBS="$FRAMEWORKS/Python.framework/Versions/3.10/lib"
export DYLD_FRAMEWORK_PATH="$FRAMEWORKS${DYLD_FRAMEWORK_PATH:+:$DYLD_FRAMEWORK_PATH}"
export DYLD_LIBRARY_PATH="$DYLIBS${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"
exec "$SCRIPT_DIR/casanovo" "$@"
EOF
chmod +x casanovo_env_v5/bin/casanovo_run.sh
```

If the checkpoint or environment are missing, the Spectrum to Sequence button will return an error. The Sequence to Spectrum feature will still work without them.

---

## Testing the App

### Sequence to Spectrum

1. Open `http://localhost:8000`
2. Select **Sequence to Spectrum**
3. Upload `sample.fasta` (included in the repo)
4. Click **Run Prediction**
5. A spectrum bar chart appears showing predicted fragment ions for PEPTIDEK

### Spectrum to Sequence

Requires the Casanovo setup above.

1. Select **Spectrum to Sequence**
2. Upload `casanovo/predicted_spectra_noseq.mgf` (included in the repo)
3. Click **Run Prediction** (takes ~5–10 seconds on CPU)
4. Results show predicted peptide sequences with confidence scores for each spectrum

---

## Running Tests

**Unit tests** (no internet or Casanovo required):

```bash
python3 -m pip install pytest
python3 -m pytest tests/test_parsers.py -v
```

**Koina integration tests** (requires internet):

```bash
python3 tests/test_inference.py
```

---

## Known Limitations

- **Casanovo is local-only.** The model weights and environment are not distributed with the repo and must be set up manually. There is no cloud fallback for spec2seq.
- **Koina uses the public API.** seq2spec calls `koina.proteomicsdb.org`, which is rate-limited and may be slow for large batches. For private data or high throughput, run a self-hosted Koina server (see `DEPLOYMENT.md`).
- **seq2spec shows only the first sequence.** If you upload a FASTA file with multiple sequences, only the first peptide's spectrum is displayed in the chart. All sequences are processed by the backend.
- **Charge is fixed at 2+.** The FASTA parser assumes charge state 2 for all peptides. This is a reasonable default for tryptic peptides but not universally correct.
- **Not production-ready.** There is no authentication, rate limiting, or persistent storage. This is a capstone demo.
