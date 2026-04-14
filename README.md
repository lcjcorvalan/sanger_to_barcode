# sanger_to_barcode

This script converts `.ab1` Sanger sequencing files into consensus FASTA sequences and performs BLAST searches against the BOLD Systems database (v4).

## Pipeline

1. Convert `.ab1` files to FASTQ
2. Quality check 
3. Sequence trimming
4. Merge forward and reverse reads
5. Reading frame selection
6. Stop codon detection
7. BLAST search against BOLD

## Dependencies

Before running the script, make sure you have the following installed:

1. Tracy (Docker version)
2. Python 3
3. Biopython
4. R
5. R package: Biostrings

> Note: The `lgbio` user environment already includes all dependencies.

## Usage

```bash
Rscript Sanger_barcode.r Ab1.csv  
```
The Ab1.csv file must be formatted as follows:

```csv
"Fwd_file","Rev_file","Note: .ab1 files must be located in the paths specified in the CSV file"
"Ab1_file/PAV296_rbcLB-F_A01_BCPlan-06-S.ab1","Ab1_file/PAV296_rbcLB-R_A01_BCPlan-07-S.ab1","Ab1_file/PAV296_rbcLB" 

```

Note: .ab1 files must be located in the paths specified in the CSV file

# Options

The settings of parameters shoub done directly on Rscrip

cutoff <- 30      # Minimum quality threshold

window <- 4       # Sliding window size

step <- 2         # Step size

AA_lib <- "Ath_rbcL_aa.fa"  # Amino acid reference library for frame detection


# Outputs
`Prefix_final.fasta` : Final fasta

`Prefix_F_trim.fq` : Trimmed Fwd fq \n

`Prefix_R_trim.fq` : Trimmed Rev fq \n

`tbl_final.csv` : Final table with all informations


