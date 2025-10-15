# VCF to Nexus or FASTA 
Convert VCF files into nexus or fasta formats for use with the BEAST phylonco package https://github.com/bioDS/beast-phylonco

## Requirements
The following tools are required for this project:

- **Python 3** – can be installed through [Homebrew](https://brew.sh/) on Linux/macOS, or through [Chocolatey](https://chocolatey.org/) on Windows.

- **Make** – already available on most Linux distributions, macOS and Windows users may need to install it with the package manager introduced. In the following instructions, Mac users should use `gmake` instead of `make`.

## Installation:
1. Clone the repository with `git clone https://github.com/bioDS/vcf2fasta.git` to the workspace.

2. Navigate to the directory and run `make venv` and `make install`. Then required package `pysam` can be automatically installed and ready to use.
 
## Example VCF files
Example VCF files can be downloaded from https://github.com/amkozlov/cellphy-paper/blob/main/data/

The examples below use the VCF file `E15.CellPhy.vcf` from [E15.zip](https://github.com/amkozlov/cellphy-paper/blob/main/data/E15.zip)

Citation for dataset: https://doi.org/10.1186/s13059-021-02583-w

## Converting VCF to nexus

**Full diploid genotype**

Replace `<vcf file>` with the path of your VCF, and `<nexus file>` with the path of your Nexus file. 

```
python3 src/vcf2nexus.py <vcf file> <nexus file> --encoding nd16
```

Example usage:
```
python3 src/vcf2nexus.py E15.CellPhy.vcf E15.CellPhy.nex --encoding nd16
```

**Binary genotype**
```
python3 src/vcf2nexus.py <vcf file> <nexus file> --encoding binary
```

Example usage:
```
python3 src/vcf2nexus.py E15.CellPhy.vcf E15.CellPhy.binary.nex --encoding binary
```

## Converting VCF to fasta

**Full diploid genotype**

Replace `<vcf file>` with the path of your VCF, and `<fasta file>` with the path of your Fasta file. 

```
python3 src/vcf2fasta.py <vcf file> <fasta file> --encoding nd16
```

Example usage:
```
python3 src/vcf2fasta.py E15.CellPhy.vcf E15.CellPhy.fasta --encoding nd16
```

**Binary genotype**
```
python3 src/vcf2fasta.py <vcf file> <fasta file> --encoding binary
```

Example usage:
```
python3 src/vcf2fasta.py E15.CellPhy.vcf E15.CellPhy.binary.fasta --encoding binary
```
