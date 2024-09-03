import random

class Insertion:
    """
    Represents and manages insertions in a genetic sequence.

    This class is used to initialize insertions before simulating reads, store insertion
    information including coverage, reference/alternate alleles, and locations. It provides
    methods for writing insertion data to files and printing insertion information.

    Attributes:
        sequence (str): The genetic sequence without buffer regions.
        ins_len (list): Lengths of each insertion.
        ins_index (list): Indices of insertions in the sequence.
        ins_allele (list): Inserted alleles for each insertion.
        current_allele (list): Current alleles at insertion positions.
        ins_count (int): Total number of insertions.
        ins_fraction (list): Fraction of reads that should contain each insertion.
        ins_buffer_index (list): Insertion indices including the buffer region.
        added_count (list): Count of reads where each insertion was added.
        record_id (str): Identifier for the record.
        strand (str): Strand information ('+' or '-').
        total_reads_count (list): Total number of reads covering each insertion position.
        absolute_position (list): Absolute positions of insertions in the original sequence.
    """

    def __init__(self, len_set_string, sequence, absolute_start, strand, record_id, ins_fraction, buffer_region, var_indexes):
        """
        Initialize the Insertion object with given parameters.

        Args:
            sequence (str): The full genetic sequence.
            absolute_start (int): Start position of the sequence in the genome.
            strand (str): Strand information ('+' or '-').
            record_id (str): Identifier for the record.
            ins_fraction (float): Fraction of reads that should contain insertions.
            buffer_region (int): Size of the buffer region to exclude from insertion generation.
            var_indexes (list): Indexes of existing variations to avoid.
        """
        self.sequence = sequence[buffer_region:-buffer_region]
        self.ins_len = []
        self.ins_index = []
        self.ins_allele = []
        self.current_allele = []

        # len_set_string = '1*200,2*40,5*20,10*10,15*10,20*5'
        # Define possible deletion lengths and their frequencies
        # len_set = [1]*200 + [2]*40 + [5]*20 + [10]*10 + [15]*10 + [20]*5
        len_set = [int(item.split('*')[0]) for item in len_set_string.split(',') for _ in range(int(item.split('*')[1]))]

        
        # Start from a random index within the first 50 bases
        index = random.sample(range(50), 1)[0]

        while index < len(self.sequence):
            if index >= len(self.sequence):
                break

            # Avoid homopolymer regions
            while sequence[index+buffer_region] == sequence[index+buffer_region-1]:
                index -= 1
            if index < 0 or list(self.sequence[index:index+7]).count(self.sequence[index]) == 7 or index in var_indexes:
                break

            cur_ins_len = random.sample(len_set, 1)[0]
            in_allele = ''.join([random.sample(['A','T','C','G'], 1)[0] for _ in range(cur_ins_len)])
            
            self.ins_len.append(cur_ins_len)
            self.ins_allele.append(in_allele)
            self.current_allele.append(self.sequence[index])
            self.ins_index.append(index)
            var_indexes.append(index)

            # Insertions are simulated with 300-600 bp between each other
            index = index + random.sample(range(300, 600), 1)[0]

        self.ins_count = len(self.ins_len)
        self.ins_fraction = [ins_fraction] * self.ins_count
        self.ins_buffer_index = [x + buffer_region for x in self.ins_index]
        self.added_count = [0] * len(self.ins_index)
        self.record_id = record_id
        self.strand = strand
        self.total_reads_count = [0] * len(self.ins_index)
        self.absolute_position = [x + absolute_start for x in self.ins_buffer_index]

    def write_to_file(self, file_h):
        """
        Write insertion information to a file.

        Args:
            file_h (file): An open file handle to write the insertion data.
        """
        for i in range(len(self.added_count)):
            file_h.write('\t'.join([
                self.record_id,
                str(self.absolute_position[i]),
                self.current_allele[i],
                self.current_allele[i] + self.ins_allele[i],
                str(self.added_count[i]),
                str(self.total_reads_count[i])
            ]) + '\n')

    def printer(self):
        """
        Print insertion information for debugging purposes.
        """
        print(self.ins_index)
        print(self.ins_allele)
        print(self.ins_fraction)