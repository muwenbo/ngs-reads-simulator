
import random
import array
import math


class Read:
    """
    A class representing a simulated DNA read with various attributes and methods for read manipulation.

    This class handles the creation, modification, and output of simulated DNA reads,
    including the addition of errors, SNPs, insertions, and deletions.
    """

    def __init__(self, read_index, sequence, id_prefix, read1_len, read2_len, ref_chr, ref_start):
        """
        Initialize a Read object with given parameters.

        Args:
            read_index (int): The starting position of the read in the reference sequence.
            sequence (str): The reference sequence.
            id_prefix (str): Prefix for the read ID.
            read1_len (int): Length of the first read in a paired-end setup.
            read2_len (int): Length of the second read in a paired-end setup.
            ref_chr (str): Reference chromosome identifier.
            ref_start (int): Start position in the reference chromosome.
        """
        self.start_position = read_index
        self.id = '@' + id_prefix + '_SEQID' + str(read_index)
        self.read1_len = read1_len
        self.read2_len = read2_len
        self.ref_start = ref_start
        self.ref_chr = ref_chr
        
        # Calculate dynamic gap and new read length
        dynamic_gap = min(random.sample(range(200,300),1)[0], len(sequence)-read_index-read1_len-read2_len)
        new_read_len = read1_len + dynamic_gap + read2_len
        sub_sequence = sequence[read_index:read_index + new_read_len]
        
        self.seq = array.array('u', sub_sequence).tolist()
        self.seq_len = len(sub_sequence)
        self.suffix = ''
       
    def add_error(self, error, error_rate):
        """
        Add random errors to the read sequence based on the given error rate.

        Args:
            error (Error): An Error object to store error information.
            error_rate (float): The rate at which errors should be introduced.
        """
        cand_alleles = ['A','C','G','T']
        error_prob = [random.uniform(0,1) for x in range(self.seq_len)]
        error_index = [error_prob.index(x) for x in error_prob if x <= error_rate]    
        for i in error_index:
            raw_allele = self.seq[i]
            new_allele = random.sample([x for x in cand_alleles if x != raw_allele],1)[0]
            if raw_allele != new_allele:
                self.seq[i] = raw_allele
                error.add_item(self.ref_chr, self.ref_start, self.start_position, i, raw_allele, new_allele)
                self.score[i] = chr(random.sample(range(0,20),1)[0]+33)

    def score_generator(self, score_range):
        """
        Generate quality scores for the read.

        Args:
            score_range (list): A list of two integers representing the range of scores.
        """
        score_r = list(range(int(math.ceil(self.seq_len/2.0))))
        raw_score1 = [int(score_range[1] - random.gammavariate(1,0.5) - x*random.normalvariate(0.1,0.03)) for x in score_r]
        score_r.reverse()
        if (self.seq_len/2.0)%1 != 0:
            score_r = score_r[:-1]
        raw_score2 = [int(score_range[1] - random.gammavariate(1,0.5) - x*random.normalvariate(0.1,0.03)) for x in score_r]
        self.score = [chr(x+33) for x in raw_score1 + raw_score2]

    def add_snp(self, snp):
        """
        Add SNPs (Single Nucleotide Polymorphisms) to the read.

        Args:
            snp (SNP): An SNP object containing SNP information.
        """
        for index in range(len(snp.snp_buffer_index)):
            p = random.uniform(0,1)
            if((snp.snp_buffer_index[index] >= self.start_position \
                and snp.snp_buffer_index[index] < self.start_position + self.read1_len) \
                or (snp.snp_buffer_index[index] >= self.start_position + self.seq_len - self.read2_len \
                    and snp.snp_buffer_index[index] < self.start_position + self.seq_len)):
                snp.total_reads_count[index] = snp.total_reads_count[index] + 1
                if(p <= snp.snp_fraction[index]+0.01):
                    self.seq[snp.snp_buffer_index[index] - self.start_position] = snp.snp_alleles[index] 
                    snp.added_count[index] = snp.added_count[index] + 1
                    self.suffix = self.suffix + '_' + str(snp.absolute_position[index])
                              
    def add_deletion(self, dels):
        """
        Add deletions to the read.

        Args:
            dels (Deletion): A Deletion object containing deletion information.
        """
        for index in range(len(dels.del_buffer_index)):
            p = random.uniform(0,1)
            if((dels.del_buffer_index[index]+dels.del_len[index] >= self.start_position \
                and dels.del_buffer_index[index] < self.start_position + self.read1_len) \
                    or (dels.del_buffer_index[index]+dels.del_len[index] >= self.start_position + self.seq_len - self.read2_len \
                        and dels.del_buffer_index[index] < self.start_position + self.seq_len)):
                dels.total_reads_count[index] = dels.total_reads_count[index] + 1
                if(p < dels.del_fraction[index]+0.01):
                    del_index = dels.del_buffer_index[index]-self.start_position
                    if(del_index <0):
                        self.seq[:del_index+dels.del_len[index]] = '-'*(del_index+dels.del_len[index])
                    else:
                        self.seq[del_index:del_index+dels.del_len[index]] = '-' * dels.del_len[index]
                    dels.added_count[index] = dels.added_count[index]+1
                    self.suffix = self.suffix + '_' + str(dels.absolute_position[index])
        
    def add_insertion(self, ins):
        """
        Add insertions to the read.

        Args:
            ins (Insertion): An Insertion object containing insertion information.
        """
        for index in range(len(ins.ins_buffer_index)):
            p = random.uniform(0,1)
            if((ins.ins_buffer_index[index] >= self.start_position \
                and ins.ins_buffer_index[index] < self.start_position + self.read1_len) \
                    or (ins.ins_buffer_index[index] >= self.start_position + self.seq_len - self.read2_len \
                        and ins.ins_buffer_index[index] < self.start_position + self.seq_len)):
                ins.total_reads_count[index] = ins.total_reads_count[index] + 1
                if(p < ins.ins_fraction[index]+0.01):
                    ins_index = ins.ins_buffer_index[index]-self.start_position
                    if self.seq[ins_index] == '-':
                        print(self.seq)
                        print(ins.ins_buffer_index[index],self.start_position)
                        self.seq[ins_index] = str(self.seq[ins_index]) + ins.ins_allele[index]
                    ins.added_count[index] = ins.added_count[index]+1
                    self.suffix = self.suffix + '_' + str(ins.absolute_position[index])
    
    def paired_end_converter(self):
        """
        Convert the read to a paired-end format.
        """
        self.read1 = self.seq[:self.read1_len]
        base_complement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
        tmp_read2 = self.seq[-self.read2_len:]
        self.read2 = [base_complement[x] for x in tmp_read2]
        self.read2.reverse()
        self.read1_score = self.score[:self.read1_len]
        self.read2_score = self.score[-self.read2_len:]
        self.read2_score.reverse()        

    def reads_finalizer(self, error_rate, error):
        """
        Finalize the read by removing gaps, generating scores, and adding errors.

        Args:
            error_rate (float): The rate at which errors should be introduced.
            error (Error): An Error object to store error information.
        """
        index = 0
        tmp = self.seq
        while index < len(self.seq):
            if len(self.seq[index]) > 1:
                self.seq = self.seq[:index] + list(self.seq[index]) + self.seq[index+1:]
                index = index + len(self.seq[index])
            elif self.seq[index] != '-':
                index = index + 1
            else:
                self.seq.pop(index)
        self.seq_len = len(self.seq)
        self.score_generator([0,40])
        if(error_rate > 0):
            self.add_error(error, error_rate)
        self.paired_end_converter()

    def write_to_file(self, reads_file_handle):
        """
        Write the read to output files in FASTQ format.

        Args:
            reads_file_handle (object): An object with file handles for writing reads.
        """
        reads_file_handle.reads1_f.write(self.id + self.suffix + '/1\n')
        reads_file_handle.reads1_f.write(''.join(self.read1) + '\n')
        reads_file_handle.reads1_f.write('+' + '\n')
        reads_file_handle.reads1_f.write(''.join(self.read1_score) + '\n')
        reads_file_handle.reads2_f.write(self.id + self.suffix + '/2\n')
        reads_file_handle.reads2_f.write(''.join(self.read2) + '\n')
        reads_file_handle.reads2_f.write('+' + '\n')
        reads_file_handle.reads2_f.write(''.join(self.read2_score) + '\n')
    
    def printer(self):
        """
        Print read information for debugging purposes.
        """
        print(self.id)
        print(self.seq)
        print(self.score)

