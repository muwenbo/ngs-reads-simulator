
import random
import math

class SNP:
    """
    Represents and manages Single Nucleotide Polymorphisms (SNPs) in a genetic sequence.

    This class is used to initialize SNPs before simulating reads, store SNP information
    including coverage, reference/alternate alleles, and locations. It provides methods
    for adding SNP information through reads, writing SNP data to files, and printing
    SNP information.

    Attributes:
        sequence (str): The genetic sequence without buffer regions.
        snp_count (int): The number of SNPs in the sequence.
        snp_index (list): Sorted list of SNP positions in the sequence.
        snp_buffer_index (list): SNP positions including the buffer region.
        snp_alleles (list): Alternative alleles for each SNP.
        snp_fraction (list): Fraction of reads that should contain each SNP.
        added_count (list): Count of reads where each SNP was added.
        record_id (str): Identifier for the record.
        strand (str): Strand information ('+' or '-').
        total_reads_count (list): Total number of reads covering each SNP position.
        absolute_position (list): Absolute positions of SNPs in the original sequence.
    """

    def __init__(self, snp_percentage, sequence, absolute_start, strand, record_id, snp_fraction, buffer_region, var_indexes):
        """
        Initialize the SNP object with given parameters.

        Args:
            snp_percentage (float): Percentage of positions that should be SNPs.
            sequence (str): The full genetic sequence.
            absolute_start (int): Start position of the sequence in the genome.
            strand (str): Strand information ('+' or '-').
            record_id (str): Identifier for the record.
            snp_fraction (float): Fraction of reads that should contain SNPs.
            buffer_region (int): Size of the buffer region to exclude from SNP generation.
            var_indexes (list): Indexes of existing variations to avoid.
        """
        self.sequence = sequence[buffer_region:-buffer_region]
        self.snp_count = int(math.floor(len(self.sequence) * snp_percentage))
        self.snp_index = sorted(random.sample(range(len(self.sequence)), self.snp_count))
        
        # Remove SNPs in homopolymer regions or existing variation positions
        self.snp_index = [index for index in self.snp_index if 
                          list(self.sequence[index-3:index+4]).count(self.sequence[index]) != 7 and 
                          index not in var_indexes]
        self.snp_count = len(self.snp_index)
        
        self.snp_buffer_index = [x + buffer_region for x in self.snp_index]
        self.snp_alleles = self.find_allele(self.snp_index, self.sequence)
        self.snp_fraction = [snp_fraction] * self.snp_count    
        self.added_count = [0] * self.snp_count
        self.record_id = record_id
        self.strand = strand
        self.total_reads_count = [0] * self.snp_count
        self.absolute_position = [x + absolute_start for x in self.snp_buffer_index]
 
    def find_allele(self, snp_index, sequence):
        """
        Assign alternative alleles to pre-defined SNPs.

        Args:
            snp_index (list): List of SNP positions.
            sequence (str): The genetic sequence.

        Returns:
            list: List of alternative alleles for each SNP.
        """
        snp_alleles = []
        cand_alleles = ['A', 'C', 'G', 'T']
        for x in snp_index:
            cand_alleles.remove(sequence[x])
            snp_alleles.append(random.choice(cand_alleles))
            cand_alleles = ['A', 'C', 'G', 'T']
        return snp_alleles
    
    def write_to_file(self, file_handle):
        """
        Write SNP information to a file.

        Args:
            file_handle (file): An open file handle to write the SNP data.
        """
        for i in range(self.snp_count):
            file_handle.write('\t'.join([
                self.record_id,
                str(self.absolute_position[i]),
                self.sequence[self.snp_index[i]],
                self.snp_alleles[i],
                str(self.added_count[i]),
                str(self.total_reads_count[i])
            ]) + '\n')

    def printer(self):
        """
        Print SNP information for debugging purposes.
        """
        print("SNP Indexes:", self.snp_index)
        print("SNP Alleles:", self.snp_alleles)
        print("SNP Fractions:", self.snp_fraction)
