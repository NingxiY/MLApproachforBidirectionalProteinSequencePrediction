"""
test_inference.py — Validate inference_demo.py against the public Koina API.

Runs 3 example peptides through both CID and HCD models.
Uses the public endpoint by default; pass --local to test localhost:8502.

Usage:
    python test_inference.py               # public API
    python test_inference.py --local       # self-hosted on port 8502
"""

import sys
import argparse
from examples.inference_demo import predict_single, predict_batch, PUBLIC_ENDPOINT, DEFAULT_ENDPOINT

# ---------------------------------------------------------------------------
# Test peptides (bacterial proteins, diverse length and charge)
# ---------------------------------------------------------------------------
TEST_PEPTIDES = [
    {"sequence": "AAGIT",      "charge": 2, "note": "short, z=2"},
    {"sequence": "PEPTIDEK",   "charge": 2, "note": "medium, z=2, lysine C-term"},
    {"sequence": "ACDEFGHIK",  "charge": 3, "note": "medium, z=3, contains Cys"},
]


def run_tests(endpoint: str) -> None:
    print(f"\nEndpoint: {endpoint}")
    print("=" * 60)

    all_passed = True

    # ── Test 1: predict_single, CID ─────────────────────────────────────────
    print("\n[Test 1] predict_single — CID (Prosit_2020_intensity_CID)")
    for p in TEST_PEPTIDES:
        try:
            result = predict_single(
                p["sequence"], p["charge"], "CID", endpoint=endpoint
            )
            n_ions = len(result["ions"])
            top = result["ions"][:3] if result["ions"] else []
            top_str = ", ".join(f"{i['annotation']}({i['intensity']:.3f})" for i in top)
            print(f"  {p['sequence']:12s} z={p['charge']}  ions={n_ions:3d}  top: {top_str}")
            assert n_ions > 0, "No ions returned"
        except Exception as e:
            print(f"  FAIL {p['sequence']}: {e}")
            all_passed = False

    # ── Test 2: predict_single, HCD ─────────────────────────────────────────
    print("\n[Test 2] predict_single — HCD (AlphaPeptDeep_ms2_generic, CE=30, Lumos)")
    for p in TEST_PEPTIDES:
        try:
            result = predict_single(
                p["sequence"], p["charge"], "HCD", endpoint=endpoint, ce=30
            )
            n_ions = len(result["ions"])
            top = result["ions"][:3] if result["ions"] else []
            top_str = ", ".join(f"{i['annotation']}({i['intensity']:.3f})" for i in top)
            print(f"  {p['sequence']:12s} z={p['charge']}  ions={n_ions:3d}  top: {top_str}")
            assert n_ions > 0, "No ions returned"
        except Exception as e:
            print(f"  FAIL {p['sequence']}: {e}")
            all_passed = False

    # ── Test 3: predict_batch, CID ───────────────────────────────────────────
    print("\n[Test 3] predict_batch — CID, all 3 peptides in one call")
    try:
        results = predict_batch(TEST_PEPTIDES, "CID", endpoint=endpoint)
        assert len(results) == len(TEST_PEPTIDES), "Batch result count mismatch"
        for p, r in zip(TEST_PEPTIDES, results):
            assert r["sequence"] == p["sequence"], "Sequence mismatch in batch result"
            assert len(r["ions"]) > 0, f"No ions for {p['sequence']}"
            print(f"  {p['sequence']:12s}  ions={len(r['ions'])}")
        print("  Batch order preserved: OK")
    except Exception as e:
        print(f"  FAIL: {e}")
        all_passed = False

    # ── Test 4: error handling — invalid fragmentation type ──────────────────
    print("\n[Test 4] Error handling — invalid fragmentation type")
    try:
        predict_single("PEPTIDE", 2, "INVALID", endpoint=endpoint)
        print("  FAIL: should have raised ValueError")
        all_passed = False
    except ValueError as e:
        print(f"  OK: ValueError raised as expected: {e}")

    # ── Summary ──────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    if all_passed:
        print("All tests PASSED")
    else:
        print("Some tests FAILED — see output above")
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--local", action="store_true",
        help="Use local self-hosted Koina (localhost:8502) instead of public API"
    )
    args = parser.parse_args()

    endpoint = DEFAULT_ENDPOINT if args.local else PUBLIC_ENDPOINT
    run_tests(endpoint)
