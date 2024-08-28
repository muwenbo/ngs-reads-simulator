#!/bin/python

import sys
import random
import logging
from Bio import SeqIO
from sequencing_error import SequencingError
from deletion import Deletion
from insertion import Insertion
from read_file_handler import ReadsFileHandle
from snp import SNP
from read import Read
import argparse

def extract_regions_from_bed(bed_file):
    regions = []
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#'):  # Skip comment lines
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                chrom, start, end = fields[:3]
                regions.append((chrom, int(start), int(end)))
    
    return regions


def generate_reads(seq_biopython, reads_file_h, parameters):
    """
    Generate simulated reads from a given sequence with specified parameters.

    Args:
        seq_biopython (SeqRecord): A BioPython SeqRecord object containing the sequence.
        reads_file_h (ReadsFileHandle): A file handler for writing output.
        parameters (dict): A dictionary of simulation parameters.

    This function simulates reads from the given sequence, applying various genetic variations
    (SNPs, insertions, deletions) and sequencing errors based on the provided parameters.
    """
    sequence = seq_biopython.seq.upper()
    # Parse the sequence ID to extract chromosome and start position
    id_prefix = seq_biopython.id.split(':')
    chr_id = id_prefix[0]
    start_bp = int(id_prefix[1].split('-')[0].replace(',',''))

    strand = '+'
    var_indexes = []
    buffer_region = parameters['buffer_region']

    # Initialize variation objects if enabled in parameters
    if parameters['del']:
        dels = Deletion(sequence, start_bp, strand, chr_id, parameters['var_fraction'], buffer_region, var_indexes)
    if parameters['ins']:
        ins = Insertion(sequence, start_bp, strand, chr_id, parameters['var_fraction'], buffer_region, var_indexes)
    if parameters['snp']:
        snp = SNP(parameters['snp_percentage'], sequence, start_bp, strand, chr_id, parameters['var_fraction'], buffer_region, var_indexes)

    error = SequencingError()
    seq_len = len(sequence)

    # Generate reads for each position in the sequence
    for i in range(0, seq_len - buffer_region):
        read1_len = parameters['read1_len']
        read2_len = parameters['read2_len']
        multiple = parameters['coverage'] / float(read1_len + read2_len)

        while multiple > 0:
            p = random.uniform(0, 1)
            if p <= multiple:
                id_prefix_suffix = seq_biopython.id + '_' + str(multiple)            
                read = Read(i, sequence, id_prefix_suffix, read1_len, read2_len, chr_id, start_bp)

                # Apply variations if enabled
                if parameters['del']:
                    read.add_deletion(dels)
                if parameters['ins']:
                    read.add_insertion(ins)
                if parameters['snp']:
                    read.add_snp(snp)

                try:
                    read.reads_finalizer(parameters['error_rate'], error)
                except:
                    logging.error(f"Error in read finalization: ins_index={ins.ins_index}, ins_len={ins.ins_len}")
                    logging.error(f"del_index={dels.del_index}, del_len={dels.del_len}")
                    logging.error(f"var_indexes={var_indexes}")
                    raise
                read.write_to_file(reads_file_h)
            multiple = multiple - 1

    # Write variation and error information to output files
    if parameters['snp']:
        snp.write_to_file(reads_file_h.var_fh) 
    if parameters['del']:
        dels.write_to_file(reads_file_h.var_fh)
    if parameters['ins']:
        ins.write_to_file(reads_file_h.var_fh)
    error.write_to_file(reads_file_h.error_fh)


def argument_parser(argv):
    """
    Parse command-line arguments for the NGS Reads Simulator.

    Args:
        argv (list): List of command-line arguments.

    Returns:
        dict: A dictionary containing the parsed arguments and their values.
    """
    parser = argparse.ArgumentParser(description='NGS Reads Simulator')
    
    # Required arguments
    parser.add_argument('-G', '--ref_genome', required=True, help='Reference genome file')
    #parser.add_argument('-T', '--target_bed', required=True, help='Target region file')
    parser.add_argument('-S', '--var_out_file', required=True, help='Variant output file')
    parser.add_argument('-F', '--reads_out_file', required=True, help='Reads output file prefix')
    parser.add_argument('-E', '--error_out_file', required=True, help='Error output file')

    # Optional arguments with default values
    parser.add_argument('--coverage', type=int, default=50, help='Expected coverage depth')
    parser.add_argument('--error_rate', type=float, default=0.01, help='Error rate')
    parser.add_argument('--read_length', type=int, default=400, help='Estimated DNA fragment length')
    parser.add_argument('--snp_percentage', type=float, default=0.01, help='SNP percentage')
    parser.add_argument('--var_fraction', type=float, default=0.5, help='Variant fraction')
    parser.add_argument('--score_range', nargs=2, type=int, default=[0, 40], help='Sequencing score range')
    parser.add_argument('--buffer_region', type=int, default=100, help='Buffer region size')
    parser.add_argument('-1', '--read1_len', type=int, default=150, help='Read 1 length')
    parser.add_argument('-2', '--read2_len', type=int, default=150, help='Read 2 length')

    # Boolean flags for enabling different types of variations
    parser.add_argument('--snp', action='store_true', help='Enable SNP simulation')
    parser.add_argument('--deletion', action='store_true', help='Enable deletion simulation')
    parser.add_argument('--insertion', action='store_true', help='Enable insertion simulation')

    # Random seed for reproducibility
    parser.add_argument('--seed', type=int, help='Random seed')

    # Log file
    parser.add_argument('--log', type=str, default='simulation.log', help='Log file name')

    args = parser.parse_args(argv[1:])

    # Convert parsed arguments to a dictionary
    parameters = vars(args)
    # Rename 'deletion' to 'del' and 'insertion' to 'ins' for consistency
    parameters['del'] = parameters.pop('deletion')
    parameters['ins'] = parameters.pop('insertion')

    return parameters

def usage():
    """
    Print usage information for the script.
    """
    print("Usage: python script.py -G <ref_genome> -S <var_out_file> -F <reads_out_file> -E <error_out_file> [options]")
    print("For more information, use the -h or --help option.")


def main():
    """
    Main function to run the NGS Reads Simulator.
    """
    para = argument_parser(sys.argv)

    # Set up logging
    logging.basicConfig(filename=para['log'], level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')

    # Log the command-line arguments
    logging.info(f"Command-line arguments: {' '.join(sys.argv)}")

    if para['seed'] is not None:
        random.seed(para['seed'])
        logging.info(f"Random seed set to: {para['seed']}")
    else:
        random_seed = random.randint(1, 1000)
        logging.info(f"No random seed provided. Generated seed {random_seed} is used.")


    # Log important parameters
    logging.info(f"Reference genome: {para['ref_genome']}")
    logging.info(f"Coverage: {para['coverage']}")
    logging.info(f"Error rate: {para['error_rate']}")
    logging.info(f"SNP percentage: {para['snp_percentage']}")
    logging.info(f"Variant fraction: {para['var_fraction']}")

    reads_file_h = ReadsFileHandle(para['reads_out_file'], para['var_out_file'], para['error_out_file'])
    reads_file_h.file_delete()
    reads_file_h.file_open()

    logging.info('Start simulating.')
    print('Start simulating.')

    #target_regions = extract_regions_from_bed(para['target_bed'])
    #generate_variations(para['ref_genome'], target_regions, para)

    for seq_record in SeqIO.parse(para['ref_genome'], 'fasta'):
        generate_reads(seq_record, reads_file_h, para)
    reads_file_h.file_close()

    logging.info('Simulation completed.')

if __name__ == "__main__":
    main()
