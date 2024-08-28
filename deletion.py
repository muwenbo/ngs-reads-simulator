import random

class Deletion:
    """
    Represents and manages deletions in a genetic sequence.

    This class is used to initialize deletions before simulating reads, store deletion
    information including coverage, reference/alternate alleles, and locations. It provides
    methods for writing deletion data to files and printing deletion information.

    Attributes:
        sequence (str): The genetic sequence without buffer regions.
        del_len (list): Lengths of each deletion.
        del_index (list): Indices of deletions in the sequence.
        del_allele (list): Deleted alleles for each deletion.
        pre_allele (list): Alleles preceding each deletion.
        del_count (int): Total number of deletions.
        del_fraction (list): Fraction of reads that should contain each deletion.
        del_buffer_index (list): Deletion indices including the buffer region.
        added_count (list): Count of reads where each deletion was added.
        record_id (str): Identifier for the record.
        strand (str): Strand information ('+' or '-').
        total_reads_count (list): Total number of reads covering each deletion position.
        absolute_position (list): Absolute positions of deletions in the original sequence.
    """

    def __init__(self, sequence, absolute_start, strand, record_id, del_fraction, buffer_region, var_indexes):
        """
        Initialize the Deletion object with given parameters.

        Args:
            sequence (str): The full genetic sequence.
            absolute_start (int): Start position of the sequence in the genome.
            strand (str): Strand information ('+' or '-').
            record_id (str): Identifier for the record.
            del_fraction (float): Fraction of reads that should contain deletions.
            buffer_region (int): Size of the buffer region to exclude from deletion generation.
            var_indexes (list): Indexes of existing variations to avoid.
        """
        self.sequence = sequence[buffer_region:-buffer_region]
        self.del_len = []
        self.del_index = []
        self.del_allele = []
        self.pre_allele = []

        # Define possible deletion lengths and their frequencies
        len_set = [1]*200 + [2]*40 + [5]*20 + [10]*10 + [15]*10 + [20]*5 + [30]*5 + [50]*2
        
        # Start from a random index within the first 50 bases
        index = random.sample(range(50), 1)[0]

        while index < len(self.sequence):
            cur_del_len = random.sample(len_set, 1)[0]
            if cur_del_len + index >= len(self.sequence):
                break

            # Check for homopolymer regions
            if list(self.sequence[index:index+cur_del_len]).count(self.sequence[index]) == cur_del_len:
                while sequence[index+buffer_region-1] == sequence[index+buffer_region]:
                    index -= 1
                if index < 0 or list(self.sequence[index:index+7]).count(self.sequence[index]) == 7:
                    break

            self.del_len.append(cur_del_len)
            self.del_allele.append(self.sequence[index:index+cur_del_len])
            self.pre_allele.append(self.sequence[index-1])
            self.del_index.append(index)
            var_indexes.extend(range(index, index+cur_del_len))

            # Simulate next deletion 150-350bp away
            index = index + random.sample(range(150, 350), 1)[0]

        self.del_count = len(self.del_len)
        self.del_fraction = [del_fraction] * self.del_count
        self.del_buffer_index = [x + buffer_region for x in self.del_index]
        self.added_count = [0] * len(self.del_index)
        self.record_id = record_id
        self.strand = strand
        self.total_reads_count = [0] * len(self.del_index)
        self.absolute_position = [x + absolute_start - 1 for x in self.del_buffer_index]

    def write_to_file(self, file_h):
        """
        Write deletion information to a file.

        Args:
            file_h (file): An open file handle to write the deletion data.
        """
        for i in range(len(self.added_count)):
            file_h.write('\t'.join([
                self.record_id,
                str(self.absolute_position[i] - 1),
                self.pre_allele[i],
                str(self.pre_allele[i]) + '-' * self.del_len[i],
                str(self.added_count[i]),
                str(self.total_reads_count[i])
            ]) + '\n')

    def printer(self):
        """
        Print deletion information for debugging purposes.
        """
        print(self.del_index)
        print(self.del_allele)
        print(self.del_fraction)