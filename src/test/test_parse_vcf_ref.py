import os
import subprocess

import pytest
from src.vcf2fasta import (
    parse_fasta,
    translate_fasta,
    translate_vcf)

# ------------------------------
# parse_fasta tests
# ------------------------------
def test_parse_fasta_haploid(tmp_path):
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(">seq1\nACGT\n>seq2\nGGTA\n")

    result = parse_fasta(str(fasta_file))
    assert result == {"seq1": "ACGT", "seq2": "GGTA"}

def test_parse_fasta_nd16(tmp_path):
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(">seq1\n0123\n>seq2\n0321\n")

    result = parse_fasta(str(fasta_file))
    assert result == {"seq1": "0123", "seq2": "0321"}

def test_parse_fasta_binary(tmp_path):
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(">seq1\n0101\n>seq2\n1010\n")

    result = parse_fasta(str(fasta_file))
    assert result == {"seq1": "0101", "seq2": "1010"}

# allow identical chromosome with same sequences
def test_parse_fasta_duplicate(tmp_path):
    fasta_file = tmp_path / "test_dup.fasta"
    fasta_file.write_text(">seq1\nACGT\n>seq1\nACGT\n")
    result = parse_fasta(str(fasta_file))
    assert result["seq1"] == "ACGT"

# raise value error if conflict sequences for same chromosome
def test_parse_fasta_conflicting_duplicate(tmp_path):
    fasta_file = tmp_path / "test_conflict.fasta"
    fasta_file.write_text(">seq1\nACGT\n>seq1\nGGTA\n")
    with pytest.raises(ValueError):
        parse_fasta(str(fasta_file))

# raise file not found error for cannot find file
def test_parse_fasta_no_file():
    with pytest.raises(FileNotFoundError):
        parse_fasta("nonexistent.fasta")

# ------------------------------
# translate_fasta tests
# ------------------------------
def test_translate_fasta_nucleotide():
    raw_seq = "ACGT"
    assert translate_fasta(raw_seq, "nucleotide", "nd16") == "05af"

def test_translate_fasta_nd16():
    raw_seq = "05af"
    assert translate_fasta(raw_seq, "nd16", "nd16") == "05af"

# ------------------------------
# translate_vcf tests
# ------------------------------
def test_translate_vcf():
    # test nd16
    nd16 = translate_vcf("./data/example.vcf", "nd16")
    for sample_name in ["0", "1", "2", "3"]:
        assert sample_name in nd16

    for pos in ["chr1,22", "chr1,51", "chr1,100"]:
        assert any(pos in nd16[sample] for sample in nd16)

    sample3_bases = nd16["3"]
    assert sample3_bases["chr1,22"] == "K"
    assert sample3_bases["chr1,51"] == "M"
    assert sample3_bases["chr1,100"] == "f"

    # test binary
    binary = translate_vcf("./data/example.vcf", "binary")
    assert binary["3"]["chr1,22"] == "1"
    assert binary["3"]["chr1,51"] == "1"
    assert binary["3"]["chr1,100"] == "1"

# ------------------------------
# parse_vcf_ref tests
# ------------------------------
def test_parse_vcf_ref(tmp_path):
    fasta_file = "./data/chr1.fna"
    vcf_file = "./data/example.vcf"
    out_file = tmp_path / "out.fasta"  # write to a temp file

    cmd = [
        "python3", "../vcf2fasta.py",
        str(vcf_file),
        str(out_file),
        "--ref", str(fasta_file),
        "--refType", "nd16"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"Script failed: {result.stderr}"

    assert os.path.exists(out_file)

    with open(out_file) as f:
        lines = f.read().splitlines()

    headers = [line for line in lines if line.startswith(">")]
    assert len(headers) == 4, f"Expected 4 sequences, found {len(headers)}"

    for i, line in enumerate(lines):
        if line.startswith(">"):
            # Next line should be sequence
            seq_line = lines[i + 1]
            assert len(seq_line) == 100, f"Sequence after {line} is {len(seq_line)} chars, expected 100"