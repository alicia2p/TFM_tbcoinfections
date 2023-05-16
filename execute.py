import utils
import os
import logging
import datetime

logger = logging.getLogger()


# COLORS AND AND FORMATTING
"""
http://ozzmaker.com/add-colour-to-text-in-python/
The above ANSI escape code will set the text colour to bright green. The format is;
\033[  Escape code, this is always the same
1 = Style, 1 for normal.
32 = Text colour, 32 for bright green.
40m = Background colour, 40 is for black.
"""
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

def call_variants(output_dir, args, sample_list):
    '''Uses autosnippy to obtain BAM file 
        and freebayes to obtain VCF file'''

    # Set the input directory and output directory for autosnippy
    input_dir = args.inp_dir
    output_results = os.path.join(output_dir + '/autosnippy')
    utils.check_create_dir(output_results)

    # AUTOSNIPPY
    logger.info(CYAN + BOLD + "\nget .bam file and quality information with autosnippy"+ END_FORMATTING)
    # Run autosnippy to obtain .bam files and quality information
    cmd_autosnippy = "python ./autosnippy_mixtas/autosnippy.py -i %s -r ./data/MTB_ancestorII_reference.fa -o %s" %(input_dir, output_results)
    logger.debug(cmd_autosnippy)
    os.system(cmd_autosnippy)

    # Obtain the paths of .bam files and corresponding VCF files
    bam_files = [os.path.join('%s/Variants/%s/%s.bam' %(output_results, filename, filename))
                for filename in sample_list]
    vcf_files = [os.path.join('%s/%s.vcf' %(output_dir, filename))
                for filename in sample_list]

    # FREEBAYES
    logger.info(CYAN + BOLD + "\nget .vcf file with freebayes"+ END_FORMATTING)
    # Iterate over each VCF and corresponding BAM file
    for vcf, bam in zip(vcf_files, bam_files):
        if os.path.isfile(vcf):
            # If the VCF file already exists, skip the processing
            logger.info(YELLOW + "vcf file: %s EXISTS" %(vcf.split('/')[-1])+ END_FORMATTING)
            continue
        prior = datetime.datetime.now()
        logger.info(GREEN + BOLD + "run freebayes %s" %(bam.split('/')[-1])+ END_FORMATTING)
        # Run freebayes to obtain the VCF file
        cmd = 'freebayes-parallel <(fasta_generate_regions.py ./data/MTB_ancestorII_reference.fa.fai 100000) 36 -f ./data/MTB_ancestorII_reference.fa -p 1 -C 0 %s > %s' %(bam, vcf)
        logger.debug(cmd)
        os.system('/bin/bash -c "%s"' %(cmd))
        after = datetime.datetime.now()
        logger.info(CYAN + "done with function in: %s" %(after - prior)+ END_FORMATTING)

    # Return the list of VCF files
    return vcf_files

