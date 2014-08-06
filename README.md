ReadSimulator
=============

Simulate artificial reads for NGS sequencing.

0.1 
Generated single-end or paired end reads.
Add SNPs to simulated reads.
Add sequencing error ot simulated reads.

0.2
Add argument parser.
SNPs are tractable from reads name.
Reads sequencing quality score are dynammic.
The flow is wraped up to a main function.

0.3
SNP missing bug is fixed. Some SNPs are missing because the SNPs located into the gap within reads.
Coverage bug is fixed. The previous version can only have maximum 200X coverage.

0.4
The gab between paired reads are dynamic between 200bp and 300bp now.

0.5
Error generating fuction was optimized to have faster speed.

0.6
The final version for test SNP variant calling.

0.7
The version for testing 1bp deletion.
Add deletion into simulated reads.
Default parameters are wrapped into argument parser.
Disabled single-end support.
Bug existed: score of error is shifted to either sides when there is deletion on the reads.

0.8
The version for testing 2bp,5bp,10bp,15bp,20bp,30bp,50bp deletions.
Add deletion with variable length at the sametime, eg 2bp,5bp,10bp...
Bug: score of error is shifted to either sides when there is deletion on the reads.
Bug: Variations are not written to files.

0.85
The version for testing 1bp,2bp,5bp,10bp,15bp,20bp,30bp and 50bp insertions.
Sequencing score of error shift to the left side of error.

0.9
The version for testing 1bp,2bp,5bp,10bp,15bp,20bp,30bp,50bp deletions.
The abosolute position of sequencing error are wirten to a file.
Bug in score shifting is fixed.


1.0
Combined simulating SNPs, deletions and insertions.
Fix the bug that deletions at the beginning are not added correctly. 

1.1
For continous alleles, SNPs and Indels are shift to left at the simulation for easy comprison with called variants.

