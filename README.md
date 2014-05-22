ReadSimulator
=============

Simulate artificial reads for NGS sequencing.

1.0 
Generated single-end or paired end reads.
Add SNPs to simulated reads.
Add sequencing error ot simulated reads.

1.1
Add argument parser.
SNPs are tractable from reads name.
Reads sequencing quality score are dynammic.
The flow is wraped up to a main function.

1.2
SNP missing bug is fixed. Some SNPs are missing because the SNPs located into the gap within reads.
Coverage bug is fixed. The previous version can only have maximum 200X coverage.

1.3
The gab between paired reads are dynamic between 200bp and 300bp now.

1.4
Error generating fuction was optimized to have faster speed.

1.5
The final version for test SNP variant calling.

1.6
Add deletion into simulated reads.
Default parameters are wrapped into argument parser.
Disabled single-end support.

1.7
Add deletion with variable length at the sametime, eg 2bp,5bp,10bp...
