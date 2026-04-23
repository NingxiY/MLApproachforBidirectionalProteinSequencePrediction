# Koina Self-Hosted Deployment Guide

This guide explains how to run a local Koina inference server and use `inference_demo.py` to predict MS/MS fragment intensities for bacterial proteomics.

---

## Why Self-Hosted?

The public Koina API (`https://koina.proteomicsdb.org`) works for testing, but self-hosting is recommended for:
- **Data privacy** — bacterial protein sequences do not leave your infrastructure
- **Throughput** — no rate limiting; run large batches locally
- **Reproducibility** — pin a specific model version

---

## Environment Requirements

- Linux (x86_64)
- NVIDIA GPU with CUDA support (tested: L40, A100)
- One of: Docker 20.10+ with NVIDIA Container Toolkit, or Apptainer 1.0+
- Python 3.8+ with `requests` library (`pip install requests`)

---

## Option A: Docker (recommended for workstations)

```bash
# Pull the image
docker pull ghcr.io/wilhelm-lab/koina:latest

# Start the server (GPU, port 8502)
docker run --rm --gpus all \
  -p 8502:8502 \
  ghcr.io/wilhelm-lab/koina:latest
```

Verify the server is ready:
```bash
curl http://localhost:8502/v2/health/ready
# Expected: {}
```

---

## Option B: Apptainer (HPC clusters, e.g. NEU Explorer)

```bash
# Convert Docker image to SIF (one-time, ~10-30 min)
apptainer pull koina.sif docker://ghcr.io/wilhelm-lab/koina:latest

# Submit as a SLURM job (see test_koina.sh)
sbatch test_koina.sh
```

Example SLURM script (`test_koina.sh`):
```bash
#!/bin/bash
#SBATCH --partition=sharing
#SBATCH --gres=gpu:l40:1
#SBATCH --time=01:00:00

apptainer run --nv koina.sif &
sleep 60   # wait for server startup

curl http://localhost:8502/v2/health/ready
```

---

## Quick Start (Python)

```python
from inference_demo import predict_single, predict_batch, PUBLIC_ENDPOINT

# Single peptide — CID
result = predict_single("PEPTIDEK", charge=2, fragmentation="CID")
for ion in result["ions"][:5]:
    print(f"  {ion['annotation']:8s}  mz={ion['mz']:.4f}  intensity={ion['intensity']:.4f}")

# Single peptide — HCD (self-hosted)
result = predict_single(
    "PEPTIDEK", charge=2, fragmentation="HCD",
    endpoint="http://localhost:8502", ce=30
)

# Batch — CID via public API
peptides = [
    {"sequence": "AAGIT",    "charge": 2},
    {"sequence": "PEPTIDEK", "charge": 2},
]
results = predict_batch(peptides, "CID", endpoint=PUBLIC_ENDPOINT)
```

---

## Model Routing

| Fragmentation | Model | Inputs | Notes |
|---------------|-------|--------|-------|
| `CID` | `Prosit_2020_intensity_CID` | sequence, charge | No collision energy |
| `HCD` | `AlphaPeptDeep_ms2_generic` | sequence, charge, CE, instrument | Default CE=30, instrument=Lumos |

Use `CID` for bacterial data (LTQ/Orbitrap CID instruments). Use `HCD` for human/HCD data.

---

## API Format (Triton HTTP v2)

**Request** (`POST /v2/models/{model_name}/infer`):
```json
{
  "id": "request_1",
  "inputs": [
    {"name": "peptide_sequences",  "shape": [2, 1], "datatype": "BYTES",  "data": ["AAGIT", "PEPTIDEK"]},
    {"name": "precursor_charges",  "shape": [2, 1], "datatype": "INT32",  "data": [2, 2]},
    {"name": "collision_energies", "shape": [2, 1], "datatype": "FP32",   "data": [30.0, 30.0]},
    {"name": "instrument_types",   "shape": [2, 1], "datatype": "BYTES",  "data": ["Lumos", "Lumos"]}
  ]
}
```
*(CID model: omit `collision_energies` and `instrument_types`)*

**Response** (outputs array):
```json
{
  "outputs": [
    {"name": "intensities", "data": [0.95, 0.43, ...]},
    {"name": "mz",          "data": [175.12, 246.15, ...]},
    {"name": "annotation",  "data": ["y1+1", "y2+1", ...]}
  ]
}
```

**Annotation format**: `{ion_type}{position}+{charge}`, e.g. `y3+1`, `b2+2`.
`inference_demo.py` filters to z=1 b/y ions only and discards intensity = -1.0 sentinels.

---

## Example curl Commands

**Check server health:**
```bash
curl http://localhost:8502/v2/health/ready
```

**CID single peptide:**
```bash
curl -X POST http://localhost:8502/v2/models/Prosit_2020_intensity_CID/infer \
  -H "Content-Type: application/json" \
  -d '{
    "id": "test",
    "inputs": [
      {"name": "peptide_sequences", "shape": [1,1], "datatype": "BYTES",  "data": ["AAGIT"]},
      {"name": "precursor_charges", "shape": [1,1], "datatype": "INT32",  "data": [2]}
    ]
  }'
```

**HCD single peptide:**
```bash
curl -X POST http://localhost:8502/v2/models/AlphaPeptDeep_ms2_generic/infer \
  -H "Content-Type: application/json" \
  -d '{
    "id": "test",
    "inputs": [
      {"name": "peptide_sequences",  "shape": [1,1], "datatype": "BYTES", "data": ["AAGIT"]},
      {"name": "precursor_charges",  "shape": [1,1], "datatype": "INT32", "data": [2]},
      {"name": "collision_energies", "shape": [1,1], "datatype": "FP32",  "data": [30.0]},
      {"name": "instrument_types",   "shape": [1,1], "datatype": "BYTES", "data": ["Lumos"]}
    ]
  }'
```

---

## Run Tests

```bash
# Against public API (no local server required)
python test_inference.py

# Against local self-hosted server
python test_inference.py --local
```
