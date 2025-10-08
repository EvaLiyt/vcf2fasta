# '/usr/bin/env python3
"""Parse a multi-sample VCF file and produce a FASTA file."""
import os
import textwrap
import argparse
import pysam
import random
import logging

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

def setup_logger(debug=False, output_path=None):
    logger = logging.getLogger("my_parser")
    logger.setLevel(logging.DEBUG)

    if logger.hasHandlers():
        logger.handlers.clear()

    # Console handler – show only INFO or higher
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if debug and output_path:
        log_file = output_path + ".log"
        fh = logging.FileHandler(log_file, mode="w")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        logger.debug(f"Debug logging enabled. Writing to file: {log_file}")

    return logger


def parse_args():
    parser = argparse.ArgumentParser("Transform a multi-sample VCF file into a Fasta format")
    parser.add_argument("vcf", type=str, help="A multi-sample VCF file")
    parser.add_argument("fasta", type=str, help="A path to an output fasta file")
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
    parser.add_argument("--refType", type=str,
                        choices=["nd16", "nucleotide","binary"],
                        help="The sequence type of reference Fasta file. Supports 'nd16' and 'nucleotide' (default).",
                        default="nucleotide")
    parser.add_argument("--debug", type=bool,
                        help="Debug mode (True or False, default is False)",
                        default=False)
    args = parser.parse_args()
    return args


def vcf2fasta(vcf, fasta, encoding="nd16", ref="none", refType="nucleotide", debug=False):
    if ref.lower() == "none":
        names, sequences = parse_vcf(vcf, encoding)
        write_fasta(names, sequences, fasta)
    else:
        parse_vcf_ref(vcf, fasta, encoding, ref, refType, debug)


def write_fasta(names, sequences, file):
    items = [
        "\n".join([">" + name, seq])
        for name, seq in zip(names, sequences)
    ]
    with open(file, "w") as fasta:
        fasta.write("\n\n".join(items))
        fasta.write("\n")


def parse_vcf(file, encoding):
    with pysam.VariantFile(file) as vcf:
        names = list(vcf.header.samples)
        sequences = get_sequences(vcf.fetch(), encoding)
    return names, sequences


def parse_fasta(fasta_path):
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Could not find FASTA file: {fasta_path}")

    ref_dict = {}
    seq_id = None
    seq_lines = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id:
                    seq = "".join(seq_lines)
                    if seq_id in ref_dict and ref_dict[seq_id] != seq:
                        raise ValueError(f"Duplicate sequence ID '{seq_id}' with conflicting sequences")
                    ref_dict[seq_id] = seq
                # Start new sequence
                seq_id = line[1:].split()[0]  # take first word as ID
                seq_lines = []
            else:
                seq_lines.append(line)

        # Save the last sequence
        if seq_id:
            seq = "".join(seq_lines)
            if seq_id in ref_dict and ref_dict[seq_id] != seq:
                raise ValueError(f"Duplicate sequence ID '{seq_id}' with conflicting sequences")
            ref_dict[seq_id] = seq

    if not ref_dict:
        raise FileNotFoundError(f"No sequences found in FASTA file: {fasta_path}")

    return ref_dict


def translate_fasta(ref_sequence_raw, ref_type, encoding):
    if ref_type == "nd16":
        return ref_sequence_raw
    else:
        return "".join(translate_genome(allele, encoding) for allele in ref_sequence_raw)


def translate_vcf(vcf_file, encoding):
    with pysam.VariantFile(vcf_file) as vcf:
        samples_map = {}
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


def parse_vcf_ref(vcf, fasta, encoding, ref, refType, debug=False):
    ref_dict = parse_fasta(ref)

    translated_ref_dict = {
        chrom: translate_fasta(seq, refType, encoding)
        for chrom, seq in ref_dict.items()
    }

    samples_map = translate_vcf(vcf, encoding)

    # 4. Construct full sequences per sample
    sample_names = []
    sample_sequences = []

    for sample_name, variants in samples_map.items():
        sample_names.append(sample_name)

        full_seq = []

        for chrom, ref_seq in translated_ref_dict.items():
            seq_list = list(ref_seq)  # convert to mutable list

            # Apply variants for this chromosome
            for pos_key, base in variants.items():
                chrom_v, pos = pos_key.split(",")
                pos = int(pos) - 1  # VCF positions are 1-based
                if chrom_v != chrom:
                    continue
                seq_list[pos] = base

            full_seq.append("".join(seq_list))

        sample_sequences.append("".join(full_seq))

    # 5. Write all samples to FASTA
    write_fasta(sample_names, sample_sequences, fasta)

    if debug:
        for sample_name, variants_map in samples_map.items():
            for variant_pos in variants_map:
                chrom, pos = variant_pos.split(",")
                pos = int(pos) - 1

            # Build chromosome sequence for this sample
            chrom_seq = list(translated_ref_dict[chrom])
            for variant_pos, base in variants_map.items():
                chrom_v, pos_v = variant_pos.split(",")
                pos_v = int(pos_v) - 1
                if chrom_v == chrom:
                    chrom_seq[pos_v] = base
            sequence_str = "".join(chrom_seq)

            # Random position debug
            random_index = random.randint(0, len(chrom_seq) - 1)
            vcf_pos = random_index + 1
            ref_base = translated_ref_dict[chrom][random_index]
            sample_base = chrom_seq[random_index]

            logger.debug(f"DEBUG SAMPLE for {sample_name}_{chrom}")
            logger.debug(f"Random position: {vcf_pos} (VCF 1-based) / {random_index} (Python 0-based)")
            logger.debug(f"Ref base: {ref_base}, var base: {sample_base}")
            logger.debug(f"Ref seq start: {''.join(translated_ref_dict[chrom][:10])}")
            logger.debug(f"Sample seq start: {sequence_str[:10]}")


def get_sequences(variants, encoding="nd16"):
    sequences = [variant2bases(variant, encoding) for variant in variants]
    sequences = transpose_list(sequences)
    sequences = map("".join, sequences)
    return list(sequences)


def transpose_list(lst):
    return list(map(list, zip(*lst)))


def variant2bases(variant, encoding="nd16"):
    return [sample2base(sample, encoding) for sample in variant.samples.itervalues()]


def sample2base(sample, encoding="nd16"):
    if encoding == "binary":
        alleles = [
            1 if allele is not None and allele > 1 else allele
            for allele in sample.allele_indices
        ]
    else:
        alleles = sample.alleles
    genome = ["." if allele is None else str(allele) for allele in alleles]
    genome = "".join(genome)
    return translate_genome(genome, encoding, sample.phased)


def translate_genome(genome, encoding="nd16", phased=True):
    genome = genome.upper()
    # if len(genome) != 2:
    #     raise ValueError(f"Genome must be diploid. The ploidy is: {len(genome)}")
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
    debug_mode = args.debug
    logger = setup_logger(debug=debug_mode, output_path=args.fasta)
    if debug_mode:
        log_file = args.fasta + ".log"
        print(f"Debug logging enabled. Writing to file: {log_file}")

    vcf2fasta(**vars(args))

# file = "../data/1459.vcf"
# fasta = file.replace(".vcf", ".fasta")
# encoding = "nd16"
# ref = "../data/filtered_sequence.fna"
# parse_vcf_ref(file, fasta, encoding, ref)
