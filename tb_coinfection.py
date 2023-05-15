import argparse, os, sys, utils, mixtas, execute
import logging
import datetime
import pandas as pd
import numpy as np

# parse arguments (-B)
parser = argparse.ArgumentParser(description="Script to infer a potential co-infection")

# INPUT parameters
# parser.add_argument("vcf_file", help="vcf file to extract htz positions")
# parser.add_argument("fastq_file", help="fastq file to extract htz positions")
parser.add_argument("inp_dir", help="input directory that contains the fastq files")
parser.add_argument("--out_dir", help="Output directory", default=".")
parser.add_argument("--file_sep", help="File separator", default="\t")
parser.add_argument("--cov_file", help='Cov file to get covered positions', 
                    default='', required=False)

# OPTIONAL
parser.add_argument("--min_DP", help="minimum frequency (depth) to accept a SNP", default=10,
                    type=int)
parser.add_argument("--min_HOM", help="minimum proportion for homocygosis",
                    default=0.85, type=float)
parser.add_argument("--ambiguity", help="min value to segregate", default=0.55, type=float)
parser.add_argument("--max_mean_htz", help="maximum htz proportion",
                    default=0.75, type=float)
parser.add_argument("--max_std_htz", help="maximum deviation of htz proportion",
                    default=0.08, type=float)
parser.add_argument("--max_extra_std_htz", help="maximum extra deviation of htz proportion",
                    default=0.015, type=float)
parser.add_argument("--SNPs_out", help="percentage of SNPs out mean htz +- (std htz + extra_std)",
                    default=0.3, type=float)
# parser.add_argument("--episode", help="fasta file to align our samples",
#                     required=False)
# parser.add_argument("--snipit", help="snipit analysis", 
#                     action="store_true")
parser.add_argument("--compare", help="directory with vcf files to validate the performance",
                    required=False, default='')

# logger

logger = logging.getLogger()

# COLORS AND AND FORMATTING

END_FORMATTING = '\033[0m'
WHITE_BG = '\033[0;30;47m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
BLUE = '\033[34m'
CYAN = '\033[36m'
YELLOW = '\033[93m'
DIM = '\033[2m'

# main
def main():

    # script directory
    abs_path = os.path.abspath(sys.argv[0])
    script_dir = os.path.dirname(abs_path)

    # arguments
    args = parser.parse_args()

    # Output directory

    output_dir = os.path.join(args.out_dir + '/results')
    utils.check_create_dir(output_dir)
    out_align_dir = os.path.join(output_dir + '/ALN/')
    utils.check_create_dir(out_align_dir)
    out_seq_dir = os.path.join(output_dir, 'fasta_sequences')
    utils.check_create_dir(out_seq_dir)
    out_annot_dir = os.path.join(output_dir, 'ANNOT/')
    utils.check_create_dir(out_annot_dir)
    out_stats_dir = os.path.join(output_dir, "Stats")
    utils.check_create_dir(out_stats_dir)

    # LOGGING
    # Create log file with date and time
    name = args.out_dir.split('/')[-1]
    right_now = str(datetime.datetime.now())
    right_now_full = "_".join(right_now.split(" "))
    log_filename = name + "_" + right_now_full + ".log"
    log_folder = os.path.join(output_dir, 'Logs')
    utils.check_create_dir(log_folder)
    log_full_path = os.path.join(log_folder, log_filename)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s:%(message)s')

    file_handler = logging.FileHandler(log_full_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    # stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    logger.info("\n\n" + BLUE + BOLD +
                "STARTING PIPELINE" + END_FORMATTING)

    today = str(datetime.date.today())

    logger.info("ARGUMENTS:")
    logger.info(str(args))

    # see filename
    sample_list = set([x.split('_R')[0] for x in os.listdir(args.out_dir) if 'fastq.gz' in x])
    
    # autosnippy and freebayes
    vcf_files = execute.call_variants(output_dir, args, sample_list)

    logger.info(CYAN + BOLD + "\nvcf parsing and strings segregation"+END_FORMATTING)
    prior = datetime.datetime.now()


    for file, vcf_file in zip(sample_list, vcf_files):

        logger.info(GREEN + "Sample: %s" %(file)+END_FORMATTING)

        # parse vcf
        df = utils.parse_vcf(args, vcf_file, output_dir, args.min_DP, file)
        if df.empty: 
            continue
        # df.to_csv(os.path.join(args.out_dir, '%s_df_to_process.csv' %(file)))

        # get alignment
        df_aln = mixtas.segregate_alignment(df, args, out_align_dir, file)

        # get sequences
        mixtas.store_sequences(df_aln, file, out_seq_dir)

        # annotate
        df_annotated = mixtas.annotation(df, df_aln, out_annot_dir, file)

        # stats
        mixtas.quality_control(df_annotated, args, file, out_stats_dir)

        # get low_pos_certain
        mixtas.low_certain_segregation(df, df_aln, out_stats_dir+'/'+file)
    
    after = datetime.datetime.now()
    logger.info(CYAN + "done with function in: %s" %(after - prior)+ END_FORMATTING)

    # if compare
    if args.compare:

        out_compare_dir = os.path.join(output_dir+'/compare')    
        utils.check_create_dir(out_compare_dir)

        logger.info(CYAN + BOLD + "\nVCFs for compare supplied"+END_FORMATTING)
        for file in sample_list:
            if os.path.isfile('%s/%s/sample_aln.csv' %(out_align_dir, file)):
                df_aln = pd.read_csv('%s/%s/sample_aln.csv' %(out_align_dir, file), header=None, index_col=0)
                df_aln = df_aln.T
                mixtas.compare_results(args, out_compare_dir, file, out_seq_dir)
                # mixtas.compare_vcf(args, out_compare_dir, df_aln, file, out_seq_dir)

    # If all OK
    logger.info(DIM + BOLD + "\n***END OF THE PIPELINE***"+END_FORMATTING)
    exit(0)



    exit(1)
#--------------------------------------------------
    # Output directory
    output_dir = os.path.join(args.out_dir + '/results_dp_%s' %(args.min_DP))
    utils.check_create_dir(output_dir)

    # autosnippy and freebayes
    vcf_file, filename = execute.call_variants(output_dir, args)

    # parse vcf
    logger.info(CYAN + BOLD + "\nvcf parsing and strings segregation"+END_FORMATTING)
    prior = datetime.datetime.now()

    df = utils.parse_vcf(args, vcf_file, output_dir, args.min_DP)

    # get alignment
    df_aln = mixtas.segregate_alignment(df, args, output_dir)

    # get sequences
    mixtas.store_sequences(df_aln, filename, output_dir)

    # get low_pos_certain
    mixtas.low_certain_segregation(df, df_aln, output_dir)

    # annotate
    df_annotated = mixtas.annotation(df, df_aln, output_dir)

    # stats
    mixtas.quality_control(df_annotated, args, filename, output_dir)

    # # if episode
    # if args.episode:
    #     mixtas.compare_episode(args, output_dir)
    
    after = datetime.datetime.now()
    logger.info(CYAN + "\ndone with function in: %s" %(after - prior)+ END_FORMATTING)

    # if compare
    if args.compare:
        logger.info(CYAN + BOLD + "\nVCF for compare supplied"+END_FORMATTING)
        mixtas.compare_vcf(args, output_dir, df_aln)

    # If all OK
    logger.info(DIM + BOLD + "\n***END OF THE PIPELINE***"+END_FORMATTING)
    exit(0)

if __name__ == "__main__":
    main()
    



    

