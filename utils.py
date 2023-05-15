import logging
import re
import sys
import os
import numpy as np
import os
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
import logging

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


#################################################################
########################## create dir ###########################
#################################################################


def check_create_dir(path):
    if os.path.exists(path): pass
    else: os.mkdir(path)

#################################################################
##################### extract sample name #######################
#################################################################

def extract_read_list_legacy(input_dir):
    """
    Search files in a directory sort by name and extract comon name of R1 and R2
    with extract_sample() function
    190615 - Limit only parent folder, not subdirectories
    """
    input_dir = os.path.abspath(input_dir)
    r1_list = []
    r2_list = []
    for root, _, files in os.walk(input_dir):
        if root == input_dir:  # This only apply to parent folder, not subdirectories
            for name in files:
                filename = os.path.join(root, name)
                is_fasta = re.match(r'.*\.f(ast)*[aq](\.gz)*', name)
                r1 = re.match(
                    r'.*(_R1_|_1|_1_|_R1).*\.f(ast)*[aq](\.gz)*$', name)
                r2 = re.match(
                    r'.*(_R2_|_2|_2_|_R2).*\.f(ast)*[aq](\.gz)*$', name)
                if is_fasta:
                    if r1:
                        r1_list.append(filename)
                    elif r2:
                        r2_list.append(filename)
                    else:
                        logger.info(
                            RED + "ERROR, file is not R1 nor R2" + END_FORMATTING)
                        sys.exit(1)
    r1_list = sorted(r1_list)
    r2_list = sorted(r2_list)
    return r1_list, r2_list


def extract_sample(R1_file, R2_file):
    """
    Extract sample from R1, R2 files.
    """
    basename_R1 = os.path.basename(R1_file)
    basename_R2 = os.path.basename(R2_file)

    sample_name_R = os.path.commonprefix([basename_R1, basename_R2])

    long_suffix = re.search('_S.*', sample_name_R)
    short_suffix = re.search('_R.*', sample_name_R)
    bar_suffix = re.search('_$', sample_name_R)
    dot_suffix = re.search('.R$', sample_name_R)

    if long_suffix:
        match = long_suffix.group()
        sample_name = sample_name_R.split(match)[0]
    elif short_suffix:
        match = short_suffix.group()
        sample_name = sample_name_R.split(match)[0]
    elif bar_suffix:
        match = bar_suffix.group()
        sample_name = sample_name_R.rstrip("_")
    elif dot_suffix:
        match = dot_suffix.group()
        sample_name = sample_name_R.rstrip(".R")
    else:
        sample_name = sample_name_R

    return sample_name



#################################################################
########################## VCF -> TSV ###########################
#################################################################

def parse_fasta(file):
    with(open(file, 'r')) as f:
        f = f.readlines()

    gen = [y for x in f if x[0] != '>' for y in x if y != '\n']
    header = f[0].split()

    return header, gen

def discard_SNP_in_DEL(df):
    '''discards those SNP that are in delected regions'''
    indel_pos = []
    pos_remove = []
    for index, row in df.iterrows():
        if row['TYPE'] == 'del': 
            indel_pos += [x for x in range(row['POS'], row['POS']+len(row['REF']))]   
        elif row['TYPE'] == 'snp' and row['POS'] in indel_pos: 
            pos_remove.append(index)
    df.drop(index = pos_remove, inplace=True)


def import_VCF42_freebayes_to_tsv(vcf_file, output_dir, sep='\t'):
    
    vcf_file = os.path.abspath(vcf_file)
    name_vcf = os.path.basename(vcf_file)
    name_tsv = name_vcf.split('.')[0]

    tsv_file = output_dir + '/' + name_tsv + ".tsv"

    if os.path.isfile(tsv_file):
        logger.info(YELLOW + "tsv file: %s EXISTS" %(name_tsv + ".tsv")+ END_FORMATTING) 
        return tsv_file

    headers = []
    extra_fields = ['TYPE', 'DP', 'RO', 'AO']
    with open(tsv_file, 'w+') as fout:
        with open(vcf_file, 'r') as f:
            next_line = f.readline().strip()
            while next_line.startswith("#"):
                next_line = f.readline().strip()
                if next_line.startswith('#CHROM'):
                    headers = next_line.split('\t')

        headers = headers[:7] + extra_fields + ['OLDVAR']
        fout.write(("\t").join(headers) + "\n")

        with open(vcf_file, 'r') as f:
            for line in f:
                extra_field_list = []
                # and not 'complex' in line and not 'mnp' in line
                if not line.startswith("#"):
                    line_split = line.split(sep)[:8]
                    info = line_split[-1].split(";")
                    for field in extra_fields:
                        extra_field_list.append(
                            [x.split("=")[-1] for x in info if field in x][0])
                    if 'OLDVAR' in line:
                        extra_field_list.append([x.split("=")[-1]
                                                 for x in info if 'OLDVAR' in x][0].split(',')[0])
                    output_line = ("\t").join(
                        line_split[:7] + extra_field_list)
                    fout.write(output_line + "\n")
    
    return tsv_file

def import_tsv_freebayes_to_df(tsv_file, dp, sep='\t'):

    df = pd.read_csv(tsv_file, sep=sep)

    if df.empty:
        logger.info(RED + BOLD + 'There is not information about this sample: %s' %(tsv_file.split('/')[-1]) + END_FORMATTING)
        return pd.DataFrame()

    df.rename(columns={'#CHROM': 'REGION', 'RO': 'REF_DP',
                       'DP': 'TOTAL_DP', 'AO': 'ALT_DP', 'QUAL': 'ALT_QUAL'}, inplace=True)

    df = df[df['TOTAL_DP'] > dp]

# after modify snippy, it started to raise a new error in which we had a line in df['ALT_DP'] with a ','
# this code splits this row with sep=',', adding a new row for each comma to the df

    try:
        df['REF_FREQ'] = df['REF_DP']/df['TOTAL_DP']
        df['ALT_FREQ'] = df['ALT_DP']/df['TOTAL_DP']

    except:
        # to split lines with ',' in ALT_DP column
        l = -1
        for index, row in df.iterrows():
            if ',' in str(row['ALT_DP']):
                
                l+=1
                
                for i in range(len(row['ALT_DP'].split(','))):
                    if i == 0: df.loc[index] = pd.Series([row[x] if ',' not in str(row[x]) else row[x].split(',')[i] for x in df.columns], index=df.columns)
                    else: 
                        new_row = pd.DataFrame([row[x] if ',' not in str(row[x]) else row[x].split(',')[i] for x in df.columns], index=df.columns)
                        df = pd.concat([df, new_row.transpose()])

        to_int = ['POS', 'REF_DP', 'ALT_DP', 'len_AD']

        for column in df.columns:
            if column in to_int:
                df[column] = df[column].astype(int)

        df['REF_FREQ'] = df['REF_DP']/df['TOTAL_DP']
        df['ALT_FREQ'] = df['ALT_DP']/df['TOTAL_DP']

    df = df.sort_values(by=['POS']).reset_index(drop=True)

    return df[['REGION', 'POS', 'ID', 'REF', 'ALT', 'ALT_QUAL', 'FILTER', 'TOTAL_DP', 'TYPE', 'REF_DP', 'ALT_DP', 'REF_FREQ', 'ALT_FREQ', 'OLDVAR']]

def different_nt(string1, string2):
    '''finds the different str between two strings
    used in remove_extra_nt'''
    pos=0
    for x, y in zip(string1, string2):
        if x != y: return [x, y, pos]
        pos+=1

def remove_extra_nt(df):
    '''detects the SNP in rows with lenght > 1'''
    for index, x in df.iterrows():
        if len(x['REF']) > 1:
            ref, alt, pos = different_nt(x['REF'], x['ALT'])
            df.at[index,'REF']= ref
            df.at[index, 'ALT']= alt
            df.at[index, 'POS']= x['POS']+pos

def cov_selection(df, args, sample_file):

    cov_files = [os.path.join(args.cov_file,x) for x in os.listdir(args.cov_file) 
                 if x.split('.')[0] in sample_file]

    if len(cov_files) != 2:
        logger.info(RED + BOLD + "there is a problem with arg cov_file: 2 cov files are required"+ END_FORMATTING)

        exit(1)

    df_cov = pd.DataFrame()

    for file in cov_files:
        filename = file.split('/')[-1]
        df_cov_int = pd.read_csv(file, sep='\t', names=['region', 'POS', 'COV'])
        df_cov_int = df_cov_int[df_cov_int['COV'] >= 40] # minimum coverage that original reads must have

        df_cov['POS_%s' %(filename)] = df_cov_int['POS']
        df_cov['COV_%s' %(filename)] = df_cov_int['COV']
    
    df_cov = df_cov.fillna(0)
    df_cov['cov_difference']= abs(df_cov.iloc[:, 1] - df_cov.iloc[:, 3])
    df_cov = df_cov[df_cov['cov_difference'].isin(np.arange(11))] # maximum reads difference between both samples

    df = df[df['POS'].isin(list(df_cov.iloc[:, 0]))]

    return df

def parse_vcf(args, vcf_file, output_dir, dp, filename, compare=False):

    '''Obtain df format to do the alignment, annotation and its analysis'''
    
    # transform from VCF to TSV and from TSV to df
    tsv_file = import_VCF42_freebayes_to_tsv(vcf_file, output_dir)
    df = import_tsv_freebayes_to_df(tsv_file, dp=dp)

    if df.empty:
        return df    

    # remove those SNP that are in delections
    discard_SNP_in_DEL(df)

    # # get positions with a dp value => dp variable
    # df_int = df[['POS', 'TOTAL_DP']]
    # df_int = df_int.groupby(['POS']).sum()
    # df_int = df_int[df_int['TOTAL_DP'] >= dp]
    # df_int = df_int.reset_index()
    # pos = list(df_int['POS'])
    # df = df[df['POS'].isin(pos)]


    df.drop_duplicates(subset=['POS', 'REF', 'ALT'], keep='first', inplace = True)

    if args.cov_file: # just with in vitro samples
        df = cov_selection(df, args, filename)

    # get the SNPs
    df = df[df.TYPE == 'snp']


    # Drop non relevant features
    df.drop(columns=["REGION", 'ID', "FILTER", 'TYPE', 'OLDVAR'], 
                        inplace=True)
    
    
    # Order columns
    df = df[["POS", "REF", "ALT", "TOTAL_DP", "REF_DP", "REF_FREQ",
            "ALT_DP", "ALT_FREQ"]]

    # Round ALT_FREQ
    # df['ALT_FREQ'] = df['ALT_FREQ'].astype(float).round(3)

    # Fill the columns
    df['REF_FREQ'] = 1 - df['ALT_FREQ']
    df['REF_DP'] = df['TOTAL_DP'] - df['ALT_DP']

    # Remove repetitive regions and extra nt when len(row) > 1
    df = remove_pe(df)

    if compare == False and len(df) < 10:
        logger.info(YELLOW + BOLD + "There is less than 10 snps in this sample. Analysis will not be performed" + END_FORMATTING)
        return []

    # Adapt format of 'REF' and 'ALT' values with length > 1
    remove_extra_nt(df)

    # Round ALT_FREQ and REF_FREQ
    df['ALT_FREQ'] = df['ALT_FREQ'].astype(float).round(3)
    df['REF_FREQ'] = df['REF_FREQ'].astype(float).round(3)

    return df


#################################################################
######################### REMOVE PE/PEE #########################
#################################################################


def bed2tsv_ranges(out_dir, bedfile):

    '''Generates a tsv file with the bed file ranges (not used)'''

    with open(bedfile, 'r') as file:
        ranges = file.read().split('\n')

    ranges = [x.split('\t') for x in ranges ]
    highs = np.array([x[2] for x in ranges if len(x)>1])
    lows = np.array([x[1] for x in ranges if len(x)>1])

    df = pd.DataFrame(data=(lows, highs))
    df = df.T

    df.to_csv("%srepeats_phage_coord.tsv" %(out_dir + "/"), sep="\t", header=False, index=False)


# function used by remove_pe
def in_range(x, lows, highs):
    return np.any((lows <= x) & (x <= highs))

# def remove_pe(df):

#     '''Remove PE/PEE regions'''

#     with open('./data/repeats_phage_coord.tsv', 'r') as file:
#         ranges = file.read().split('\n')

#     ranges = [x for x in ranges if len(x) != 0]
#     ranges = [y.split() for y in ranges]
#     highs = np.array([int(y[1]) for y in ranges])
#     lows = np.array([int(y[0]) for y in ranges])

#     for index, x in df.iterrows():
#         if in_range(x['POS'], lows, highs):
#             df = df.drop(index=index)
    
#     return df

def remove_pe(df):

    with open('./data/repeats_phage_coord.tsv', 'r') as file:
        ranges = file.read().split('\n')

    pos_range = [rang for x in ranges if len(x)>1 for rang in range(int(x.split('\t')[0]), int(x.split('\t')[1])+1)]

    df = df[~df['POS'].isin(pos_range)]

    return df
