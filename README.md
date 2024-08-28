# NGS-Reads-Simulator

NGS-Reads-Simulator is a tool for simulating artificial reads for Next-Generation Sequencing (NGS). It allows you to generate simulated reads with various genetic variations, including SNPs, insertions, and deletions.

## Features

- Simulate reads from a given target region
- Generate SNPs, insertions, and deletions
- Simulate sequencing errors
- Customizable read parameters

## Installation

[Add installation instructions here]

## Usage

### Preparing the target region

Before running ReadSimulator, you need to prepare a FASTA file of all target regions. The header should contain the genomic location of the target region, for example:

```shell
samtools faidx hg19.fa chr13:32,908,718-32,916,660 > target.fa
samtools faidx target.fa
python ReadSimulator/ReadsSimulator.py -G target.fa -S simulated_variants.out -F simulate.reads -E simulated_sequencing_error.txt --snp --ins --del 
```

