import pytest
from src.vcf2fasta import parse_fasta, translate_fasta, translate_vcf, vcf2fasta

@pytest.mark.parametrize(
    "fasta_content,expected",
    [
        (">seq1\nACGT\n>seq2\nGGTA\n", {"seq1": "ACGT", "seq2": "GGTA"}),  # haploid
        (">seq1\n0123\n>seq2\n0321\n", {"seq1": "0123", "seq2": "0321"}),  # nd16
        (">seq1\n0101\n>seq2\n1010\n", {"seq1": "0101", "seq2": "1010"}),  # binary
    ]
)
def test_parse_fasta(tmp_path, fasta_content, expected):
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)

    result = parse_fasta(str(fasta_file))
    assert result == expected


@pytest.mark.parametrize(
    "fasta_content,expected",
    [
        (">seq1\nACGT\nACGT\n>seq2\nGGTA\nACGT\n", {"seq1": "ACGTACGT", "seq2": "GGTAACGT"})
    ]
)
def test_parse_fasta_multiple_lines_sequence(tmp_path, fasta_content, expected):
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)

    result = parse_fasta(str(fasta_file))
    assert result == expected


@pytest.mark.parametrize(
    "fasta_content,should_raise",
    [
        (">seq1\nACGT\n>seq1\nACGT\n", False),
        (">seq1\nACGT\n>seq1\nGGTA\n", True),
    ]
)
def test_parse_fasta_identical_names(tmp_path, fasta_content, should_raise):
    """
    Test parsing a fasta files with identical chromosome names.
    - allow identical chromosome with same sequences, will return one sequence for each chromosome
    - raise value error if conflict sequences for same chromosome name

    :param tmp_path:
    :param fasta_content:
    :param should_raise:
    :return:
    """
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)

    if should_raise:
        with pytest.raises(ValueError):
            parse_fasta(str(fasta_file))
    else:
        result = parse_fasta(str(fasta_file))
        assert list(result.values())[0] in ["ACGT"]


@pytest.mark.parametrize(
    "fasta_path,error",
    [
        ("notExisting.fasta", FileNotFoundError),
    ]
)
def test_parse_fasta_no_file(fasta_path, error):
    with pytest.raises(error):
        parse_fasta(fasta_path)


@pytest.mark.parametrize(
    "fasta_content,error",
    [
        ("", FileNotFoundError),
    ]
)
def test_parse_fasta_no_seq(tmp_path, fasta_content, error):
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)

    with pytest.raises(error):
        parse_fasta(str(fasta_file))


@pytest.mark.parametrize(
    "seq,src_type,target_type,expected",
    [
        ("ACGT", "nucleotide", "nd16", "05af"),
        ("05af", "nd16", "nd16", "05af"),
    ]
)
def test_translate_fasta(seq, src_type, target_type, expected):
    """
    Test translate fasta file with nucleotide and nd16 dataType.
    - Nucleotide given, translate to nd16
    - Nd16 given, copy directly
    :param seq:
    :param src_type:
    :param target_type:
    :param expected:
    :return:
    """
    assert translate_fasta(seq, src_type, target_type) == expected


@pytest.mark.parametrize(
    "output_type,expected",
    [
        ("nd16", {"chr1,22": "K", "chr1,51": "M", "chr1,100": "f"}),
        ("binary", {"chr1,22": "1", "chr1,51": "1", "chr1,100": "1"})
    ]
)
def test_translate_vcf_param(output_type, expected):
    """
    Test translate vcf param with correct output type.
    :param output_type:
    :param expected:
    :return:
    """
    result = translate_vcf("data/example.vcf", output_type)
    for pos, base in expected.items():
        assert result["3"][pos] == base


def check_sequences(sequences, expect_num, expect_length):
    assert len(sequences) == expect_num
    for seq in sequences.values():
        assert len(seq) == expect_length


def test_vcf2fasta_single_ref(tmp_path):
    """
    Test vcf2fasta() for single reference sequence and verify sequence translation.
    """
    fasta_file = "./data/chr1.fna"
    vcf_file = "./data/example.vcf"
    out_file = tmp_path /"out.fasta"

    vcf2fasta(vcf_file, out_file, encoding="nd16", ref=fasta_file, refType="nd16", debug=False)

    sequences = parse_fasta(out_file)
    check_sequences(sequences, 4, 100)


def test_vcf2fasta_multiple_refs(tmp_path):
    """
    Test vcf2fasta() for a VCF file with two reference sequences.
    - check correct number of sequences in output file
    - check each sequence length
    - check specific genotypes at given positions
    """
    fasta_file = "./data/twoSeqRef.fna"
    vcf_file = "./data/twoSeqExample.vcf"
    out_file = tmp_path /"output.fasta"

    vcf2fasta(vcf_file, out_file, ref=fasta_file)

    sequences = parse_fasta(out_file)
    check_sequences(sequences, 8, 100)

    expected_positions = {
        ("0", "chr1"): {14: "?"},
        ("0", "chr2"): {20: "?"},
        ("3", "chr1"): {49: "M"},
        ("3", "chr2"): {87: "f"},
    }

    def check_genotype_at_position(seq, position, expected_base):
        actual = seq[position]
        assert actual == expected_base

    for seq_name, seq in sequences.items():
        parts = seq_name.split("_")

        if len(parts) < 2:
            raise ValueError(f"Invalid sequence name format: {seq_name}")

        sample, chrom = parts[0], parts[1]
        key = (sample, chrom)

        if key in expected_positions:
            for pos, base in expected_positions[key].items():
                check_genotype_at_position(seq, pos, base)