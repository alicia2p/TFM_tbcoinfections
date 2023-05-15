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
    '''Use autosnippy to obtain .bam file and do variant calling with freebayes'''

    input_dir = args.inp_dir
    output_results = os.path.join(output_dir + '/autosnippy')

    utils.check_create_dir(output_results)
    logger.info(CYAN + BOLD + "\nget .bam file and quality information with autosnippy"+ END_FORMATTING)
    cmd_autosnippy = "python ./autosnippy_mixtas/autosnippy.py -i %s -r ./data/MTB_ancestorII_reference.fa -o %s" %(input_dir, output_results)
    # cmd_autosnippy = "python ./autosnippy_mixtas/autosnippy.py -i %s -r ./data/Mycobacterium_abscessus_ATCC_19977_complete_genome.fasta -o %s" %(input_dir, output_results)
    logger.debug(cmd_autosnippy)
    os.system(cmd_autosnippy)


    bam_files = [os.path.join('%s/Variants/%s/%s.bam' %(output_results, filename, filename))
                for filename in sample_list]
    vcf_files = [os.path.join('%s/%s.vcf' %(output_dir, filename))
                for filename in sample_list]
    print('bam_files', bam_files)

    logger.info(CYAN + BOLD + "\nget .vcf file with freebayes"+ END_FORMATTING)

    for vcf, bam in zip(vcf_files, bam_files):
        if os.path.isfile(vcf):
            logger.info(YELLOW + "vcf file: %s EXISTS" %(vcf.split('/')[-1])+ END_FORMATTING)
            continue
        prior = datetime.datetime.now()
        logger.info(GREEN + BOLD + "run freebayes %s" %(bam.split('/')[-1])+ END_FORMATTING)
        cmd = 'freebayes-parallel <(fasta_generate_regions.py ./data/MTB_ancestorII_reference.fa.fai 100000) 36 -f ./data/MTB_ancestorII_reference.fa -p 1 -C 0 %s > %s' %(bam, vcf)
        # cmd = 'freebayes-parallel <(fasta_generate_regions.py ./data/Mycobacterium_abscessus_ATCC_19977_complete_genome.fasta.fai 100000) 36 -f ./data/Mycobacterium_abscessus_ATCC_19977_complete_genome.fasta -p 1 -C 0 %s > %s' %(bam, vcf)
        logger.debug(cmd)
        os.system('/bin/bash -c "%s"' %(cmd))
        after = datetime.datetime.now()
        logger.info(CYAN + "done with function in: %s" %(after - prior)+ END_FORMATTING)

    return vcf_files
    
    # input_dir = args.inp_dir
    # output_results = os.path.join(output_dir + '/autosnippy')
    
    # files = [x for x in os.listdir(input_dir) if 'fastq' in x]
    # filename = files[0][0:files[0].find('_R')]

    # bam_file = os.path.join('%s/Variants/%s/%s.bam' %(output_results, filename, filename))
    # vcf_file = os.path.join('%s/%s.vcf' %(output_dir, filename))

    # if os.path.isfile(bam_file):
    #     logger.info(GREEN + "\nbam file already exists, continue with the pipe"+ END_FORMATTING)
    # else:
        # utils.check_create_dir(output_results)
        # logger.info(CYAN + BOLD + "\nget .bam file and quality information with autosnippy"+ END_FORMATTING)
        # cmd_autosnippy = "python ./autosnippy_mixtas/autosnippy.py -i %s -r ./data/MTB_ancestorII_reference.fa -o %s" %(input_dir, output_results)
        # logger.debug(cmd_autosnippy)
        # os.system(cmd_autosnippy)
    

    # if os.path.isfile(vcf_file):
    #     logger.info(GREEN + "\nvcf file already exists, continue with the pipe"+ END_FORMATTING) 
    #     return vcf_file, filename
    # else:
    #     prior = datetime.datetime.now()
    #     logger.info(CYAN + BOLD + "\nrun freebayes"+ END_FORMATTING)
    #     cmd = 'freebayes-parallel <(fasta_generate_regions.py ./data/MTB_ancestorII_reference.fa.fai 100000) 36 -f ./data/MTB_ancestorII_reference.fa -p 1 -C 0 %s > %s' %(bam_file, vcf_file)
    #     logger.debug(cmd)
    #     os.system('/bin/bash -c "%s"' %(cmd))
    #     after = datetime.datetime.now()
    #     logger.info(CYAN + "\ndone with function in: %s" %(after - prior)+ END_FORMATTING)

    # return vcf_file, filename

