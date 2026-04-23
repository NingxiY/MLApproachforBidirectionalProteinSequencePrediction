"""
test_parsers.py — pytest tests for parsers.parse_fasta()
"""

import pytest
from parsers import parse_fasta


# ── Valid input cases ──────────────────────────────────────────────────────────

def test_single_record():
    data = b">sp|P12345|MYG_HUMAN Myoglobin\nMGLSDGEWQLVL\n"
    result = parse_fasta(data)
    assert result == [{"sequence": "MGLSDGEWQLVL", "charge": 2}]


def test_multiple_records():
    data = b">seq1\nPEPTIDEK\n>seq2\nAAAKPEAAAK\n"
    result = parse_fasta(data)
    assert result == [
        {"sequence": "PEPTIDEK",   "charge": 2},
        {"sequence": "AAAKPEAAAK", "charge": 2},
    ]


def test_multiline_sequence_is_concatenated():
    data = b">seq1\nMGLSD\nGEWQL\nVLNVW\n"
    result = parse_fasta(data)
    assert result == [{"sequence": "MGLSDGEWQLVLNVW", "charge": 2}]


def test_mixed_case_is_uppercased():
    data = b">seq1\nPeptideK\n"
    result = parse_fasta(data)
    assert result[0]["sequence"] == "PEPTIDEK"


def test_blank_lines_between_records_are_ignored():
    data = b">seq1\n\nACDEF\n\n>seq2\n\nGHIKL\n"
    result = parse_fasta(data)
    assert result == [
        {"sequence": "ACDEF", "charge": 2},
        {"sequence": "GHIKL", "charge": 2},
    ]


def test_charge_is_always_2():
    data = b">seq1\nPEPTIDEK\n>seq2\nAAAKPEAAAK\n"
    result = parse_fasta(data)
    assert all(p["charge"] == 2 for p in result)


def test_no_header_plain_text_sequence():
    # A plain .txt file with just a sequence and no ">" header
    # has no record delimiters, so the whole content becomes one sequence.
    data = b"PEPTIDEK\n"
    result = parse_fasta(data)
    assert result == [{"sequence": "PEPTIDEK", "charge": 2}]


# ── Error cases ────────────────────────────────────────────────────────────────

def test_empty_input_raises():
    with pytest.raises(ValueError, match="No valid sequences found"):
        parse_fasta(b"")


def test_header_only_raises():
    with pytest.raises(ValueError, match="No valid sequences found"):
        parse_fasta(b">header one\n>header two\n")


def test_headers_with_only_blank_lines_raises():
    with pytest.raises(ValueError, match="No valid sequences found"):
        parse_fasta(b">header\n\n   \n\n>another\n\n")
