# '/usr/bin/env python3
"""Parse a multi-sample VCF file and produce a FASTA file."""
import argparse
import textwrap
from typing import Tuple, List, Dict, Optional

import pysam

nd16phased = {
    "AA": "0",
    "AC": "1",
    "AG": "2",
    "AT": "3",
    "CA": "4",
    "CC": "5",
    "CG": "6",
    "CT": "7",
    "GA": "8",
    "GC": "9",
    "GG": "a",
    "GT": "b",
    "TA": "c",
    "TC": "d",
    "TG": "e",
    "TT": "f",
    # Phased with one unknown
    ".A": "g",
    ".C": "h",
    ".G": "i",
    ".T": "j",
    "A.": "n",
    "C.": "o",
    "G.": "p",
    "T.": "q",
    "..": "?",
    # haploid states are homozygous genotypes
    "a": "0",
    "A": "0",
    "c": "5",
    "C": "5",
    "g": "a",
    "G": "a",
    "t": "f",
    "T": "f"
}

nd16unphased = {
    "AA": "0",
    "CC": "5",
    "GG": "a",
    "TT": "f",
    "AC": "M",
    "CA": "M",
    "AG": "R",
    "GA": "R",
    "AT": "W",
    "TA": "W",
    "CG": "S",
    "GC": "S",
    "CT": "Y",
    "TC": "Y",
    "GT": "K",
    "TG": "K",
    "..": "?",
    # haploid states are homozygous genotypes
    "a": "0",
    "A": "0",
    "c": "5",
    "C": "5",
    "g": "a",
    "G": "a",
    "t": "f",
    "T": "f"
}

binary = {
    "00": "0",
    "11": "1",
    "01": "1",
    "10": "1",
    ".1": "1",
    "1.": "1",
    "0.": "?",
    ".0": "?",
    "..": "?"
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser("Transform a multi-sample VCF file into a Fasta format")
    parser.add_argument("vcf", type=str, help="A multi-sample VCF file")
    parser.add_argument("fasta", type=str, help="A path to an output fasta file")
    parser.add_argument("--readable", action="store_true",
                        help="""
                        The way to write fasta file. 
                        When flag is on, wrap fasta file 80 characters per line. Suits human and some editors reading.
                        Default off, no breaks within a single sequence for bioinformatic tools reading.
                        """)
    parser.add_argument("--encoding", type=str, choices=["nd16", "binary"], default="nd16",
                        help="""
    A type of encoding used during to translate a diploid variant into a single character.
    Currently available encodings are:
        binary         - reference allele is 0 and any non-reference allele is mutation encoded as 1
        nd16 (default) - full nucleotide diploid model, available in phased and unphased form which
                         is chosen based on the phasing of the VCF file.
    """)
    parser.add_argument("--ref", type=str,
                        help="A reference Fasta file to fill in genotypes not present in the VCF",
                        default="none")
    parser.add_argument("--ref_type", type=str,
                        choices=["nucleotide", "nd16"],
                        help="The sequence type of reference Fasta file. "
                             "Supports 'nd16' and 'nucleotide' (default).",
                        default="nucleotide")
    parser.add_argument("--debug", action="store_true",
                        help="Debug mode (default: off)")
    args = parser.parse_args()
    return args


def vcf2fasta(vcf, fasta, readable = False, encoding="nd16", ref="none",
              ref_type="nucleotide", debug=False) -> None:
    """
    The main function of vcf2fasta. If reference is given,
    then map the full-length sequence in fasta file.
    Otherwise, the fasta file only contains variable sites.
    :param vcf:
    :param fasta:
    :param readable:
    :param encoding:
    :param ref:
    :param ref_type:
    :param debug:
    :return:
    """
    if ref.lower() == "none":
        names, sequences = parse_vcf(vcf, encoding)
        write_fasta(names, sequences, fasta, readable)
    else:
        parse_vcf_ref(vcf, fasta, encoding, ref, ref_type, debug)


def write_fasta(names, sequences, file, readable=False) -> None:
    """
    Write fasta file, default generating a sequence per line for bioinformatic tools reading.
    When readable flag is on, wrap fasta file 80 characters per line for editors and human reading.
    :param names:
    :param sequences:
    :param file:
    :param readable:
    :return:
    """
    if readable:
        items = [
            "\n".join([">" + name, textwrap.fill(seq, width=80)])
            for name, seq in zip(names, sequences)
        ]
    else:
        items = [
            "\n".join([">" + name, seq])
            for name, seq in zip(names, sequences)
        ]

    with open(file, "w") as fasta:
        fasta.write("\n\n".join(items))
        fasta.write("\n")


def parse_vcf(file, encoding) -> Tuple[List[str], List[str]]:
    """
    Parse the vcf file and return the sample names and sequences representing variable sites.
    :param file:
    :param encoding:
    :return:
    """
    with pysam.VariantFile(file) as vcf:
        names = list(vcf.header.samples)
        sequences = get_sequences(vcf.fetch(), encoding)
    return names, sequences


def read_fasta(fasta_path: str) -> List[Tuple[str, str]]:
    """
    Read a fasta file and return a list of sequences.
    Allow a single sequences for multiple lines.
    Raises error when file does not exist or empty.
    """
    sequences: List[Tuple[str, str]] = []
    seq_id: Optional[str] = None
    seq_lines: List[str] = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id:
                    sequences.append((seq_id, "".join(seq_lines)))
                seq_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)

        if seq_id:
            sequences.append((seq_id, "".join(seq_lines)))

    if not sequences:
        raise FileNotFoundError(f"No sequences found in fasta file: {fasta_path}")

    return sequences


def validate_sequences(seq_list: List[Tuple[str, str]]) -> Dict[str, str]:
    """
    Validate sequences: join lines into full sequence
    Allow identical chromosome names
    - extract one sequence when sequences are the same,
    - raise error when sequences are not the same.
    """
    validated: Dict[str, str] = {}
    for seq_id, seq in seq_list:
        if seq_id in validated:
            if validated[seq_id] != seq:
                raise ValueError(f"Duplicate sequence ID '{seq_id}' with conflicting sequences")
        else:
            validated[seq_id] = seq
    return validated


def parse_fasta(fasta_path: str) -> Dict[str, str]:
    """
    Read and validate the pass in fasta file, return a list of reference sequences.
    :param fasta_path:
    :return:
    """
    seq = read_fasta(fasta_path)
    ref_seq = validate_sequences(seq)
    return ref_seq


def translate_fasta(seq: str, ref_type: str, encoding: str) -> str:
    """
    Translate reference sequence into encoded type.
    :param seq:
    :param ref_type:
    :param encoding:
    :return:
    """
    if ref_type == "nd16":
        return seq

    translated_list: List[str] = [translate_genome(allele, encoding) for allele in seq]
    translated_seq: str = "".join(translated_list)

    return translated_seq


def translate_vcf(vcf_file, encoding) -> Dict[str, Dict[str, str]]:
    """
    Translate variants in vcf file into encoded type,
    store the chrom name, sample name, positions in the sample map.
    :param vcf_file:
    :param encoding:
    :return:
    """
    with pysam.VariantFile(vcf_file) as vcf:
        samples_map: Dict[str, Dict[str, str]] = {}
        for variant in vcf.fetch():
            chrom = variant.chrom  # store chromosome in pos_key
            for sample in variant.samples.values():
                base = sample2base(sample, encoding)
                sample_name = sample.name
                pos_key = f"{chrom},{variant.pos}"
                if sample_name not in samples_map:
                    samples_map[sample_name] = {}
                samples_map[sample_name][pos_key] = base
    return samples_map


def parse_vcf_ref(vcf, fasta, encoding, ref, ref_type, debug=False, readable=False) -> None:
    """
    Read and translate variant genotype in vcf file and reference sequences.
    Map the missing sites in vcf with reference bases.
    :param vcf:
    :param fasta:
    :param encoding:
    :param ref:
    :param ref_type:
    :param debug:
    :param readable:
    :return: output fasta file
    """
    ref_dict = parse_fasta(ref)

    translated_ref_dict = {
        chrom: translate_fasta(seq, ref_type, encoding)
        for chrom, seq in ref_dict.items()
    }

    samples_map = translate_vcf(vcf, encoding)

    fasta_entries = map_ref_sites(samples_map, translated_ref_dict)

    write_fasta(
        [name for name, _ in fasta_entries],
        [seq for _, seq in fasta_entries],
        fasta,
        readable
    )


def map_ref_sites(samples_map: Dict[str, Dict[str, str]],
                  translated_ref_dict: Dict[str, str]) -> List[Tuple[str, str]]:
    """
    Map the missing sites in vcf with reference bases.
    :param samples_map:
    :param translated_ref_dict:
    :return:
    """
    fasta_entries: List[Tuple[str, str]] = []
    for sample_name, variants in samples_map.items():
        for chrom, ref_seq in translated_ref_dict.items():
            seq_list = list(ref_seq)

            # apply variants for this chromosome
            for pos_key, base in variants.items():
                chrom_v, pos = pos_key.split(",")
                if chrom_v != chrom:
                    continue
                idx: int = int(pos) - 1  # VCF positions are 1-based
                if 0 <= idx < len(seq_list):
                    seq_list[idx] = base

            final_seq = "".join(seq_list)
            fasta_entries.append((f"{sample_name}_{chrom}", final_seq))
    return fasta_entries


def get_sequences(variants, encoding="nd16") -> List[str]:
    """
    Return list of sequences representing the variable sites.
    :param variants: variants from vcf file
    :param encoding:
    :return: a list of sequences corresponding to each sample
    """
    sequences_row: List[List[str]] = [variant2bases(variant, encoding) for variant in variants]
    sequences_column = transpose_list(sequences_row)
    sequences: List[str] = ["".join(seq) for seq in sequences_column]
    return sequences


def transpose_list(lst: List[List[str]]) -> List[List[str]]:
    """
    Reconstruct the lists by columns (original rows)
    :param lst: a list of lists, each element represents a vcf row
    :return: a list of lists, each element represents a vcf column (each sample)
    """
    return list(map(list, zip(*lst)))


def variant2bases(variant, encoding="nd16") -> List[str]:
    """
    Convert all variants genotype to encoded bases.
    :param variant:
    :param encoding:
    :return:
    """
    bases = [sample2base(sample, encoding) for sample in variant.samples.values()]
    return bases


def sample2base(sample, encoding="nd16") -> str:
    """
    Convert a single sample's genotype to the encoded type
    :param sample:
    :param encoding:
    :return:
    """
    if encoding == "binary":
        alleles = [
            1 if allele is not None and allele > 1 else allele
            for allele in sample.allele_indices
        ]
    else:
        alleles = sample.alleles
    genome_list = ["." if allele is None else str(allele) for allele in alleles]
    genome: str = "".join(genome_list)

    translated_genome = translate_genome(genome, encoding, sample.phased)
    return translated_genome


def translate_genome(genome, encoding="nd16", phased=True) -> str:
    """
    Translate a genome sequence into encoded type.
    :param genome:
    :param encoding:
    :param phased:
    :return:
    """
    genome = genome.upper()
    if len(genome) not in [1,2]:
        raise ValueError(f"Genome must be diploid. The ploidy is: {len(genome)}")
    if encoding == "binary":
        return binary[genome]
    if encoding == "nd16" and phased:
        return nd16phased[genome]
    if encoding == "nd16" and not phased:
        return nd16unphased[genome]
    else:
        raise ValueError(f"Encoding \"{encoding}\" is not supported.")



if __name__ == "__main__":
    args = parse_args()
    vcf2fasta(**vars(args))
