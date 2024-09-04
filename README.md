# NGS-Reads-Simulator

NGS-Reads-Simulator is a tool for simulating artificial reads for Next-Generation Sequencing (NGS). It allows you to generate simulated reads with various genetic variations, including SNPs, insertions, and deletions.

## Features

- Simulate reads from a given target region
- Generate SNPs, insertions, and deletions
- Simulate sequencing errors
- Customizable read parameters


## Usage
```shell
optional arguments:
  -h, --help            show this help message and exit
  -G REF_GENOME, --ref_genome REF_GENOME
                        Reference genome file
  -S VAR_OUT_FILE, --var_out_file VAR_OUT_FILE
                        Variant output file
  -F READS_OUT_FILE, --reads_out_file READS_OUT_FILE
                        Reads output file prefix
  -E ERROR_OUT_FILE, --error_out_file ERROR_OUT_FILE
                        Error output file
  --coverage COVERAGE   Expected coverage depth
  --error_rate ERROR_RATE
                        Error rate
  --read_length READ_LENGTH
                        Estimated DNA fragment length
  --snp_percentage SNP_PERCENTAGE
                        SNP percentage regarding to the length of given sequence. 0.01 represents 1 SNP per 100 base pairs
  --deletion_length DELETION_LENGTH
                        Define the expected length of deletion and count of each length, e.g. 1*200,2*40 represents 200 1bp deletions and 40 2bp deletions.
  --insertion_length INSERTION_LENGTH
                        Define the expected length of insertion and count of each length, e.g. 1*200,2*40 represents 200 1bp insertions and 40 2bp insertions.
  --var_fraction VAR_FRACTION
                        Variant fraction
  --score_range SCORE_RANGE SCORE_RANGE
                        Sequencing score range
  --buffer_region BUFFER_REGION
                        Buffer region size
  -1 READ1_LEN, --read1_len READ1_LEN
                        Read 1 length
  -2 READ2_LEN, --read2_len READ2_LEN
                        Read 2 length
  --snp                 Enable SNP simulation
  --deletion            Enable deletion simulation
  --insertion           Enable insertion simulation
  --seed SEED           Random seed
  --log LOG             Log file name
```

## Preparing the target region

Before running ReadSimulator, you need to prepare a FASTA file of all target regions. The header should contain the genomic location of the target region, for example:

```shell
samtools faidx hg19.fa chr13:32,908,718-32,916,660 > target.fa
samtools faidx target.fa
python main.py -G target.fa -S simulated_variants.out -F simulate.reads -E simulated_sequencing_error.txt --snp --insertion --deletion
```

