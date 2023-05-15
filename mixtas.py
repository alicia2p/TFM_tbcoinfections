import os
import math
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import seaborn as sns
import utils
import warnings
warnings.filterwarnings("ignore")
import logging

mlogger = logging.getLogger('matplotlib')
mlogger.setLevel(logging.WARNING)

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
######################### GET ALIGNMENT #########################
#################################################################

def segregate_alignment(df, args, output_dir, filename):

    # homocygotic positions

    df_hom = df[df['ALT_FREQ'] >= args.min_HOM]
    df_hom.iloc[:,[2]]

    df_hom['sample_1'] = list(df_hom['ALT'])
    df_hom['sample_1+2'] = list(df_hom['ALT'])
    df_hom['sample_2'] = list(df_hom['ALT'])
    df_hom = df_hom.rename(columns={'REF': 'Reference'})

    df_hom = df_hom[['POS', 'Reference', 'sample_1', 'sample_1+2', 'sample_2']]

    # ambiguous positions

    df_htz = df[(df['ALT_FREQ'] <= args.ambiguity) & (df['ALT_FREQ'] > (1-args.ambiguity))].round(3)

    iupac =     {"AG": "R", "GA":"R",
                "CT":"Y", "TC":"Y",
                "GC":"S", "CG":"S",
                "AT":"W", "TA":"W",
                "GT":"K", "TG":"K",
                "AC":"M", "CA":"M"}

    sample_1 = ['%s (%s)' %(iupac[ref+alt], min([alt_freq, ref_freq]))  
                for alt, ref, alt_freq, ref_freq 
                in zip(df_htz['ALT'], df_htz['REF'], 
                    df_htz['ALT_FREQ'], df_htz['REF_FREQ'])]
    
    sample_2 = ['%s (%s)' %(iupac[ref+alt], max([alt_freq, ref_freq]))  
                for alt, ref, alt_freq, ref_freq 
                in zip(df_htz['ALT'], df_htz['REF'], 
                    df_htz['ALT_FREQ'], df_htz['REF_FREQ'])]
    
    sample_1_2 = ['%s/%s' %(ref, alt)
                for ref, alt in zip (df_htz['REF'], df_htz['ALT'])]
    
    df_htz['sample_1'] = sample_1
    df_htz['sample_1+2'] = sample_1_2
    df_htz['sample_2'] = sample_2
    df_htz = df_htz.rename(columns={'REF': 'Reference'})

    df_htz = df_htz[['POS', 'Reference', 'sample_1', 'sample_1+2', 'sample_2']]

    # HTZ positions

    df_segregate = df[ ((df.ALT_FREQ < args.min_HOM) & 
                        (df.ALT_FREQ > args.ambiguity)) |
                        ((df.ALT_FREQ <= (1 - args.ambiguity)) &
                        (df.ALT_FREQ >= (1 - args.min_HOM)))]

    sample_1 = ['%s (%s)' %(alt, min([alt_freq, ref_freq]))
                if min([alt_freq, ref_freq]) == alt_freq
                else '%s (%s)' %(ref, min([alt_freq, ref_freq])) 
                for alt, ref, ref_freq, alt_freq 
                in zip(df_segregate['ALT'], df_segregate['REF'], 
                    df_segregate['REF_FREQ'], df_segregate['ALT_FREQ'])]

    sample_2 = ['%s (%s)' %(alt, max([alt_freq, ref_freq]))
                if max([alt_freq, ref_freq]) == alt_freq
                else '%s (%s)' %(ref, max([alt_freq, ref_freq])) 
                for alt, ref, ref_freq, alt_freq 
                in zip(df_segregate['ALT'], df_segregate['REF'], 
                    df_segregate['REF_FREQ'], df_segregate['ALT_FREQ'])]

    sample_1_2 = ['%s/%s' %(ref, alt)
                for ref, alt in zip (df_segregate['REF'], df_segregate['ALT'])]

    df_segregate['sample_1'] = sample_1
    df_segregate['sample_2'] = sample_2
    df_segregate['sample_1+2'] = sample_1_2
    df_segregate = df_segregate.rename(columns={'REF': 'Reference'})

    df_segregate = df_segregate[["POS","Reference","sample_1","sample_1+2", 'sample_2', 'TOTAL_DP']]


    # create df_aln

    df_aln = pd.concat([df_hom, df_htz, df_segregate])
    df_aln = df_aln.sort_values(by='POS')

    # save df_aln and df_aln_HTZ
    df_aln_t = df_aln.T
    df_aln_HTZ_t = df_aln[df_aln['sample_1+2'].str.contains('/')].T

    out_dir = os.path.join(output_dir, filename)
    utils.check_create_dir(out_dir)

    df_aln_t.to_csv(out_dir+'/sample_aln.csv', header=False)
    df_aln_HTZ_t.to_csv(out_dir+'/sample_HTZ_aln.csv', header=False)

    return df_aln



# def segregate_alignment(df, args, output_dir, filename):
    
#     # if ALT_FREQ >= min_HOM

#     hom = df[df['ALT_FREQ'] >= args.min_HOM]

#     hom = hom[['POS', 'REF','ALT','TOTAL_DP']]

#     hom['sample_1'] = hom[['ALT']]
#     hom['sample_1+2'] = hom[['ALT']]
#     hom = hom.rename(columns={'ALT':'sample_2', 'REF':'Reference'})

#     hom = hom[["POS","Reference","sample_1","sample_1+2", 'sample_2', 'TOTAL_DP']]

#     # if ALT_FREQ between ambiguity and 1-ambiguity 

#     htz = df[df['ALT_FREQ'] <= args.ambiguity].round(3)
#     htz = htz[htz['ALT_FREQ'] > 1 - args.ambiguity].round(3)


#     iupac =     {"AG": "R", "GA":"R",
#                 "CT":"Y", "TC":"Y",
#                 "GC":"S", "CG":"S",
#                 "AT":"W", "TA":"W",
#                 "GT":"K", "TG":"K",
#                 "AC":"M", "CA":"M"}

#     sample_1 = []
#     sample_2 = []
#     sample_1_2 = []

#     for index, x in htz.iterrows():
#         nt = iupac[x['ALT']+x['REF']]
#         sample_2.append(nt +' ('+ str(max(x[['REF_FREQ', 'ALT_FREQ']]))+')')
#         sample_1.append(nt +' ('+ str(min(x[['REF_FREQ', 'ALT_FREQ']]))+')')
#         sample_1_2.append(x['REF']+'/'+x['ALT'])

#     htz = htz[['POS','REF','TOTAL_DP']]

#     htz['sample_1'] = sample_1
#     htz['sample_1+2'] = sample_1_2
#     htz['sample_2'] = sample_2
#     htz = htz.rename(columns={'REF':'Reference'})

#     htz = htz[["POS","Reference","sample_1","sample_1+2", 'sample_2', 'TOTAL_DP']]

#     # if ALT_FREQ between min_HOM and ambiguity

#     segregate = df[(df['ALT_FREQ'] < args.min_HOM) & (df['ALT_FREQ'] > args.ambiguity)] # & (result['ALT_FREQ'] <= 0.45) & (result['ALT_FREQ'] >= 0.15)]
#     segregate2 = df[(df['ALT_FREQ'] <= 1- args.ambiguity) & (df['ALT_FREQ'] > 1 - args.min_HOM)] # dani quita el igual en >= 1-min_hom ya que no considera los que tengan 0.15

#     sample_1 = []
#     sample_2 = []
#     sample_1_2 = []

#     for index, x in segregate.iterrows():
#         sample_2.append(x['ALT'] + ' (' + str(x['ALT_FREQ']) + ')')
#         sample_1.append(x['REF'] + ' (' + str(x['REF_FREQ']) + ')')
#         sample_1_2.append(x['ALT'] + '/' + x['REF'])

#     segregate = segregate[['POS','REF','TOTAL_DP']]

#     segregate['sample_1'] = sample_1
#     segregate['sample_1+2'] = sample_1_2
#     segregate['sample_2'] = sample_2
#     segregate = segregate.rename(columns={'REF':'Reference'})

#     segregate = segregate[["POS","Reference","sample_1","sample_1+2", 'sample_2', 'TOTAL_DP']]

#     sample_1 = []
#     sample_2 = []
#     sample_1_2 = []

#     for index, x in segregate2.iterrows():
#         sample_2.append(x['REF'] + ' (' + str(x['REF_FREQ']) + ')')
#         sample_1.append(x['ALT'] + ' (' + str(x['ALT_FREQ']) + ')')
#         sample_1_2.append(x['REF'] + '/' + x['ALT'])

#     segregate2 = segregate2[['POS','REF','TOTAL_DP']]

#     segregate2['sample_1'] = sample_1
#     segregate2['sample_1+2'] = sample_1_2
#     segregate2['sample_2'] = sample_2
#     segregate2 = segregate2.rename(columns={'REF':'Reference'})

#     segregate2 = segregate2[["POS","Reference","sample_1","sample_1+2", 'sample_2', 'TOTAL_DP']]

#     segregate = pd.concat([segregate, segregate2])

#     # create df_aln

#     df_aln = pd.concat([hom, htz, segregate])
#     df_aln = df_aln.sort_values(by='POS')

#     # save df_aln and df_aln_HTZ
#     df_aln_t = df_aln.T
#     df_aln_HTZ_t = df_aln[df_aln['sample_1+2'].str.contains('/')].T

#     out_dir = os.path.join(output_dir, filename)
#     utils.check_create_dir(out_dir)

#     df_aln_t.to_csv(out_dir+'/sample_aln.csv', header=False)
#     df_aln_HTZ_t.to_csv(out_dir+'/sample_HTZ_aln.csv', header=False)

#     return df_aln

def store_sequences(df_aln, filename, output_dir):

    # obtain ref_gen, sample_1 and sample_2

    with(open('./data/MTB_ancestorII_reference.fa', 'r')) as f:
        f = f.readlines()

    ref_gen = [y for x in f if x[0] != '>' for y in x if y != '\n']
    # header = f[0].split()


    d_sample_1 = {int(y):x[0] for y,x in zip(list(df_aln['POS']), df_aln['sample_1'])}
    d_sample_2 = {int(y):x[0] for y,x in zip(list(df_aln['POS']), df_aln['sample_2'])}

    

    sample_1 = [x if int(y) not in d_sample_1.keys() else d_sample_1[y] 
                for x, y in zip(ref_gen, np.arange(1, len(ref_gen)+1))]
    
    sample_2 = [x if int(y) not in d_sample_2.keys() else d_sample_2[y] 
                for x, y in zip(ref_gen, np.arange(1, len(ref_gen)+1))]

    df_store = pd.DataFrame({'sample_1':sample_1, 'sample_2':sample_2})

    # store sequences to fasta files
    out_seq_dir = os.path.join(output_dir, filename)
    utils.check_create_dir(out_seq_dir)

    for sample in ['sample_1', 'sample_2']:
        seq = open(out_seq_dir + '/' + '%s.fasta' %(sample), "w")
        to_write = '>%s_%s\n' %(filename, sample) + ''.join(list(df_store[sample]))
        seq.write(to_write)
        seq.close()


def low_certain_segregation(df, df_aln, output_dir):

    # Compute mean htz proportion and std of majority allel
    HTZ_SNVs = df[(df.ALT_FREQ <= 0.85) & (df.ALT_FREQ >= (1 - 0.85))]
    # number htz SNPs
    n_HTZ_SNPs = HTZ_SNVs.shape[0]
    if n_HTZ_SNPs:
        # select upper proportion htz
        upper_HTZ_prop_l = HTZ_SNVs[["ALT_FREQ", "REF_FREQ"]].max(axis=1).to_list()
        # mean_htz and std
        mean_ALT_HTZ_prop = round(np.mean(upper_HTZ_prop_l), 2)
        std_ALT_HTZ_prop = round(np.std(upper_HTZ_prop_l), 2)
    else:
        mean_ALT_HTZ_prop = 0
        std_ALT_HTZ_prop = 0

    max_HTZ = df_aln[df_aln['sample_2'].str.contains('\(')]

    # Check if SNP is of low certain segregation

    freq = [float(x[x.find('(')+1:x.find(')')]) for x in list(max_HTZ['sample_2'])]
    max_HTZ.insert(6, 'FREQ', freq)
    max_HTZ = max_HTZ[['POS', 'FREQ']]

    pos_low_certain = max_HTZ[max_HTZ['FREQ'] > mean_ALT_HTZ_prop + (std_ALT_HTZ_prop + 0.015)]
    pos_low_certain_2 = max_HTZ[max_HTZ['FREQ'] < mean_ALT_HTZ_prop - (std_ALT_HTZ_prop + 0.015)]

    pos_low_certain = pd.concat([pos_low_certain, pos_low_certain_2])

    # save pos_low_certain postiions and its frequency

    moderate_conf_segregat_snps = output_dir + '/moderate_confidence_segregated_SNPs.csv'

    pos_low_certain.to_csv(moderate_conf_segregat_snps, index=False)


#################################################################
########################### ANNOTATION ##########################
#################################################################

def make_autopct(values):
    '''function used in pie_plot_annotation autopct'''
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100))
        return '{p:.2f}%  ({v:d})'.format(p=pct,v=val)
    return my_autopct

def pie_plot_annotation(samples, output_dir):
    '''Plots the fraction of each lineage in each sample'''
    
    try:
        
        fig, ax = plt.subplots(1, 2, figsize=(11, 6), )
        sns.set_theme(palette="Set2", font_scale= 1)

        n=0
        for sample in list(samples.columns):
            
            ax[n].pie(  samples[sample], 
                        radius=1,
                        wedgeprops=dict(width=0.4, 
                                        edgecolor='white'), 
                        autopct=make_autopct(list(samples[sample])),
                        labels=list(samples.index),  
                        rotatelabels=False, 
                        labeldistance=None, 
                        frame=False, 
                        pctdistance=1.3, 
                        textprops=dict( color="black", 
                                        size=11, 
                                        weight='normal' ))
            
            ax[n].legend(bbox_to_anchor=(1, 0))
            
            ax[n].set_title(sample, y=1, pad=10)
            n+=1

        plt.subplots_adjust(bottom=0.1, right=1, top=0.7)
        plt.savefig(output_dir+'/annotation_pie_plot.png',
                        bbox_inches='tight')
    except:
        logger.info(YELLOW + BOLD + "unable to construct pie plot: check matplotlib version"+ END_FORMATTING)

def bar_plot_annotation(samples, output_dir):
    
    try:
        samples_T = samples.T

        ax = samples_T.plot(kind = 'bar',
                            stacked = True,
                            title = 'number of marker SNPs in each sample',
                            rot = 0)
        
        show_values=[]
        for column in samples_T.columns:
            show_values = show_values + list(samples_T[column])

        for rect, value in zip(ax.patches, show_values):
            if value != 0:
                h = rect.get_height() /2.
                w = rect.get_width() /2.
                x, y = rect.get_xy()
                ax.text(x+w, y+h,value,horizontalalignment='center',verticalalignment='center')

        ax.legend(bbox_to_anchor=(1.65, 1))

        plt.savefig("%s/bar_plot_annotation.png" %(output_dir), bbox_inches='tight')

    except:
        logger.info('')
        logger.info(YELLOW + BOLD + "unable to construct bar plot: check matplotlib version"+ END_FORMATTING)


def annotation(df, df_aln, output_dir, filename):

    '''Annotate snps markers using ./data/annotation.tsv'''

    annot = pd.read_csv('./data/annotation.tsv', sep='\t')
    annot = annot.fillna('')
    annot = annot.drop(columns='Unnamed: 0')
    annot['ANNOT'] = annot['species'] + '_' + annot['lineage'] + '_' + annot['sublineage']
    annot['ANNOT'] = ['_'.join(x.split(' ')) for x in list(annot['ANNOT'])]
    annot = annot.rename(columns={'SNP':'ALT'})

    annotpie = annot[['POS', 'ALT', 'ANNOT']]
    annotated = annotpie[annotpie.POS.isin(df['POS'])]
    annotated = annotated.sort_values(by=['POS'])

    out_dir = os.path.join(output_dir + filename)
    utils.check_create_dir(out_dir)

    samples = ['sample_1', 'sample_2']
    c = 1

    for sample in samples:
        
        snp = [x[0] for x in df_aln[sample]]
        pos = list(df_aln['POS'])

        colname = {'POS':pos, 'ALT':snp}

        sample = pd.DataFrame(colname)
        sample = pd.merge(sample, annotated, how="left", on=["POS", "ALT"])
        
        sample_save = sample
        sample_save.to_csv(out_dir+'/sample_'+ str(c) +'_annot.csv', index=False)
        
        sample = sample.dropna()
        sample = sample.sort_values(by='ANNOT')
        sample.to_csv(out_dir+'/sample_'+ str(c) +'_annot_only.csv', index=False)

        sample_list = list(sample['ANNOT'])

        if c ==1 : 
            d = {x:[sample_list.count(x)] for x in set(sample_list)}

        else: 
            d2 = {x:[sample_list.count(x)] for x in set(sample_list)}

        c+=1

    sample_1 = pd.DataFrame(d, index=['sample_1']).T
    sample_2 = pd.DataFrame(d2, index=['sample_2']).T

    samples = sample_1.join(sample_2, how='outer')
 
    samples.to_csv(out_dir+'/amount_snp_annot.csv')

    samples = samples[['sample_1', 'sample_2']].fillna(0)

    pie_plot_annotation(samples, out_dir)
    bar_plot_annotation(samples, out_dir)


    df_annotated = pd.merge(df, annot, how="left", on=["POS", "ALT"])
    df_annotated_dropna = df_annotated.dropna()

    df_annotated.to_csv(out_dir+'/snps_annotated.csv', index=False)
    df_annotated_dropna.to_csv(out_dir+'/snps_annotated_only.csv', index=False)

    return df_annotated


#################################################################
############################# STATS #############################
#################################################################

def split(a, n):
    """ Splits a list into N parts of approximately equal length"""
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))

def plot_proportions(HTZ_SNVs, name_stats_file):

    '''Stacked bar plots of the proportions of the segregated samples'''

    variant_name = os.path.basename(name_stats_file)

    sample_1 = list(HTZ_SNVs[["ALT_FREQ", "REF_FREQ"]].min(axis=1).to_list())
    sample_2 = list(HTZ_SNVs[["ALT_FREQ", "REF_FREQ"]].max(axis=1).to_list())
    positions = list(HTZ_SNVs['POS'].astype(str))

    # number of graphs
    n_graphs = math.ceil(len(positions)/100)

    sample_1_s = list(split(sample_1, n_graphs))
    sample_2_s = list(split(sample_2, n_graphs))
    positions_s = list(split(positions, n_graphs))

    width = 0.50      # the width of the bars: can also be len(x) sequence

    # mean low horizontal line
    high_mean = round(np.mean(sample_2), 2)
    high_std = round(np.std(sample_2), 2)

    # start plottingreturn
    for n in range(n_graphs):

        fig, ax = plt.subplots(figsize=(70, 15))

        ax.bar(positions_s[n], sample_2_s[n], width, label='sample 2', color="lightblue")
        ax.bar(positions_s[n], sample_1_s[n], width, bottom=sample_2_s[n],
            label='sample 1', color="darkgreen")

        ax.set_ylabel('Frequency', fontsize=20, labelpad=30)
        ax.set_xlabel('Positions', fontsize=20, labelpad=30)
        ax.set_title("HTZ positions %s" %(variant_name), fontsize=30)
        ax.legend()

        plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5,
            0.6, 0.7, 0.8, 0.9, 1], fontsize=15)
        plt.xticks(rotation = 90, fontsize=15)

        # 0.5 horizontal line
        plt.axhline(y=0.5, color='r', linestyle='-')
        # mean low horizontal line
        plt.axhline(y=high_mean + high_std, color='black', linestyle='--')
        plt.axhline(y=high_mean, color='black', linestyle='-')
        plt.axhline(y=high_mean - high_std, color='black', linestyle='--')

        plt.savefig("%s_%s.png" %(name_stats_file, n))
    
    x_ticks = [x[0] for x in positions_s]

    if n_graphs <= 20:
        size = (70, 15)
        title_size = 35
        font_size = 20
        legend_size = 20
    else:
        size = (140, 30)
        title_size = 70
        font_size = 40
        legend_size = 40

    fig, ax = plt.subplots(figsize=size)

    ax.bar(positions, sample_2, width, label='sample 2', color="lightblue")
    ax.bar(positions, sample_1, width, bottom=sample_2,
        label='sample 1', color="darkgreen")

    ax.set_ylabel('Frequency', fontsize=font_size, labelpad=40)
    ax.set_xlabel('Positions', fontsize=font_size, labelpad=40)
    ax.set_title("HTZ positions", fontsize=title_size)
    ax.legend(prop={'size': legend_size})

    plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5,
        0.6, 0.7, 0.8, 0.9, 1], fontsize=font_size)
    plt.xticks(x_ticks, fontsize=font_size, rotation=90)

    # 0.5 horizontal line
    plt.axhline(y=0.5, color='r', linestyle='-', linewidth=4)
    # mean low horizontal line
    plt.axhline(y=high_mean + high_std, color='black', linestyle='--', linewidth=5)
    plt.axhline(y=high_mean, color='black', linestyle='-', linewidth=5)
    plt.axhline(y=high_mean - high_std, color='black', linestyle='--', linewidth=5)

    plt.savefig("%s_overall.png" %(name_stats_file))

def quality_control(df, args, filename, output_dir):

    '''Results stats'''

    # out_dir_stats
    dir_name_tsv_stats = os.path.join(output_dir, filename)
    utils.check_create_dir(dir_name_tsv_stats)

    # stats file
    name_stats_file = os.path.join(dir_name_tsv_stats, filename)

    if os.path.isfile("%s_stats.csv" %(name_stats_file)):
        logger.info(YELLOW + 'stats file: %s_stats.csv EXISTS' %(filename) + END_FORMATTING)
        return 
    
    stats_file = open("%s_stats.csv" %(name_stats_file), "w")

    ##### HTZ PROPORTION
    # std = 0.08
    upper_std = 0.015
    # max_prop = 0.75
    # points = 0

    # fields
    fields = ["Muestra", "total_HTZ", "mean_htz_proportion",
                "std_htz_proportion", "N_SNPs_between_std", "%_SNPs_between_std", 
                "N_SNPs_higher_mean_std", "mean_dist_higher_mean_std", "std_dist_higher_mean_std",
                "N_SNPs_lower_mean_std", "mean_dist_lower_mean_std", "std_dist_lower_mean_std", 
                "Confidence_segregation"]
    row = [filename]

    # sample_1 = [float(x[x.find('(')+1:x.find(')')]) for x in df_aln_HTZ['Sample_1']]
    HTZ_SNVs = df[(df.ALT_FREQ <= args.min_HOM) & (df.ALT_FREQ >= (1 - args.min_HOM))]
    
    #INSPECT HTZ_SNVs
    # HTZ_SNVs.to_csv(os.path.join(output_dir, 'df_HTZ_snps_analysis.csv'))

    upper_HTZ_prop_l = HTZ_SNVs[["ALT_FREQ", "REF_FREQ"]].max(axis=1).to_list()

    n_HTZ_SNPs = HTZ_SNVs.shape[0]

    if n_HTZ_SNPs == 0: # when using cov files this variable might be 0
        fields +=  ["N_SNPs_not_lineage", "%_SNPs_not_lineage", "N_SNPs_min", "N_SNPs_max"]
        row+=[str(n_HTZ_SNPs)]
        row += ['nan']*15
    
        to_write = ",".join(fields) + "\n" + ",".join(row) + "\n"

        stats_file.write(to_write)
        stats_file.close()

    else:
        # Statistics
        mean_ALT_HTZ_prop = round(np.mean(upper_HTZ_prop_l), 2)
        std_ALT_HTZ_prop = round(np.std(upper_HTZ_prop_l), 2)
        N_SNPs_between_std = len([ p for p in upper_HTZ_prop_l 
                                                    if p <= mean_ALT_HTZ_prop + (std_ALT_HTZ_prop + upper_std) and
                                                    p >= mean_ALT_HTZ_prop - (std_ALT_HTZ_prop + upper_std)])
        SNPs_in_mean_limits = round(N_SNPs_between_std / len(upper_HTZ_prop_l), 2)
        # if len(upper_HTZ_prop_l) >= 1: 
        #     SNPs_in_mean_limits = round(N_SNPs_between_std / len(upper_HTZ_prop_l), 2)
        # else: SNPs_in_mean_limits = np.nan

        # Info SNPs out interval mean +- std
        ## HIGHER
        N_SNPs_higher_mean_std = [p - (mean_ALT_HTZ_prop + std_ALT_HTZ_prop) for p in upper_HTZ_prop_l if p > mean_ALT_HTZ_prop + (std_ALT_HTZ_prop + upper_std)]
        mean_dist_higher_mean_std = round(np.mean(N_SNPs_higher_mean_std), 2)
        std_dist_higher_mean_std = round(np.std(N_SNPs_higher_mean_std), 2)

        ## LOWER
        N_SNPs_lower_mean_std = [(mean_ALT_HTZ_prop - std_ALT_HTZ_prop) - p for p in upper_HTZ_prop_l if p < mean_ALT_HTZ_prop - (std_ALT_HTZ_prop + upper_std)]
        mean_dist_lower_mean_std = round(np.mean(N_SNPs_lower_mean_std), 2)
        std_dist_lower_mean_std = round(np.std(N_SNPs_lower_mean_std), 2)

        row += [str(n_HTZ_SNPs), str(mean_ALT_HTZ_prop), str(std_ALT_HTZ_prop), str(N_SNPs_between_std), str(SNPs_in_mean_limits),
                str(len(N_SNPs_higher_mean_std)), str(mean_dist_higher_mean_std), str(std_dist_higher_mean_std),
                str(len(N_SNPs_lower_mean_std)), str(mean_dist_lower_mean_std), str(std_dist_lower_mean_std)]

        # If possible segregation
        if mean_ALT_HTZ_prop > 0.6:
            row += ["1"]
        else:
            row += ["0"]

        ##### HTZ DISTRIBUTION
        # number htz not related to lineage
        not_lineage_HTZ_SNVs = HTZ_SNVs[HTZ_SNVs['species'].isna()]
        not_lineage_HTZ_SNVs = not_lineage_HTZ_SNVs[((not_lineage_HTZ_SNVs["REF_FREQ"] > args.ambiguity) |
                                    ((not_lineage_HTZ_SNVs["ALT_FREQ"] > args.ambiguity)))]

        max_index = not_lineage_HTZ_SNVs[["ALT_FREQ", "REF_FREQ"]].idxmax(axis=1).to_list()

        fields += ["N_SNPs_not_lineage", "%_SNPs_not_lineage", "N_SNPs_min", "N_SNPs_max"]

        row += [str(len(max_index)), str(round(len(max_index)/n_HTZ_SNPs, 2)),
                str(round(max_index.count("REF_FREQ") / len(max_index), 2)),
                str(round(max_index.count("ALT_FREQ") / len(max_index), 2))]
                
        #### WRITE FILE
        
        to_write = ",".join(fields) + "\n" + ",".join(row) + "\n"
        
        stats_file.write(to_write)
        stats_file.close()

        # plot_proportions

        plot_proportions(HTZ_SNVs, name_stats_file)

        # proportion in genome's ranges
        
        part_1 = HTZ_SNVs[HTZ_SNVs['POS'] <= 1000000]
        part_2 = HTZ_SNVs[(HTZ_SNVs['POS'] > 1000000) & (HTZ_SNVs['POS'] <= 2000000)]
        part_3 = HTZ_SNVs[(HTZ_SNVs['POS'] > 2000000) & (HTZ_SNVs['POS'] <= 3000000)]
        part_4 = HTZ_SNVs[HTZ_SNVs['POS'] > 4000000]

        df_htz_div = pd.DataFrame()
        mean_div = [round(np.mean(x[["ALT_FREQ", "REF_FREQ"]].max(axis=1).to_list()), 2) for x in [part_1, part_2, part_3, part_4]]
        std_div = [round(np.std(x[["ALT_FREQ", "REF_FREQ"]].max(axis=1).to_list()), 2) for x in [part_1, part_2, part_3, part_4]]
        range_div = ['0 - 10⁶', '10⁶ - 2*10⁶', 
                    '2*10⁶ - 3*10⁶', '3*10⁶ - end']
        n_HTZ = [len(part_1), len(part_2), len(part_3), len(part_4)]

        df_htz_div['genome_range'] = range_div
        df_htz_div['number_htz'] = n_HTZ
        df_htz_div['mean_htz_proportion'] = mean_div
        df_htz_div['std_htz_proportion'] = std_div

        df_htz_div.to_csv(os.path.join(dir_name_tsv_stats, 'stats_div_genome.csv'))


#################################################################
######################## COMPARE EPISODE ########################
#################################################################



### COMPARISON DESDE CERO PORQUE YA NO PUEDO MÁS


def add_sequence(args, df_all, path_vcf, output_dir):
    
    name_vcf = path_vcf.split('/')[-1]
    name_vcf = name_vcf.split('.')[0]
    # name_vcf = name_vcf.split('_')[0]

    df_new_seq = utils.parse_vcf(args, path_vcf, output_dir, 20,
                                 name_vcf, compare=True)

    df_new_seq = df_new_seq[['POS', 'REF', 'ALT', 'TOTAL_DP', 'REF_DP', 'REF_FREQ',
       'ALT_DP', 'ALT_FREQ']]
    df_new_seq = df_new_seq[df_new_seq['ALT_FREQ'] >= 0.15]
    
    d = {pos:nt for nt, pos in zip(df_new_seq.ALT, df_new_seq.POS)}
    
    new_seq = [nt if pos not in d.keys() else d[pos] for nt, pos in zip(df_all.MTB_anc, df_all.POS)]

    df_all[name_vcf] = new_seq

    return df_new_seq, name_vcf

def fill_results(results, df_all, df_original, name_original, sample):

    # snps in original
    # df_total_snps_original = df_all[['POS', name_original]][df_all.MTB_anc != df_all[name_original]]
    df_total_snps_original = df_all[df_all.MTB_anc != df_all[name_original]]
    snps_original = len(df_total_snps_original)

    # snps in sample
    df_total_snps_sample = df_all[df_all.MTB_anc != df_all[sample]]
    snps_sample = len(df_total_snps_sample)

    # diferent snps original-sample
    df_snps_differents = df_all[['POS', sample, name_original, 'MTB_anc']][df_all[sample] != df_all[name_original]]

    df_snps_differents = pd.merge( df_snps_differents,
                                    df_original[['POS', 'ALT_FREQ', 'TOTAL_DP']],
                                    on='POS', how='left')
    snps_differents = len(df_snps_differents)

    # total positions match original-sample
    df_positions_match = df_all[df_all[sample] == df_all[name_original]]
    positions_match = len(df_positions_match)

    # snps original that match with sample
    # df_snps_original = df_all[df_all.POS.isin(df_total_snps_original.POS)]
    # df_snps_original_match_sample = df_snps_original[df_snps_original[sample] == df_snps_original[name_original]]
    # snps_match = len(df_snps_original_match_sample)

    df_snps_original_match_sample = df_total_snps_original[df_total_snps_original[sample] == df_total_snps_original[name_original]]
    snps_match = len(df_snps_original_match_sample)

    # ambiguous nt sample
    ambiguous_nt = len([x for x in df_all[sample] if x not in ['A', 'C', 'G', 'T']])

    # total snps
    total_snps = snps_match + snps_differents

    name_original_sample = '%s_%s' %(name_original, sample)

    values = [  name_original_sample,   # original-sample
                snps_original,          # snps_original
                snps_sample,            # snps_sample
                snps_match+snps_differents, # total_snps
                snps_match,             # snps_match_original
                snps_differents,        # different_alleles
                ambiguous_nt,           # ambiguous
                round(snps_match*100/snps_original, 2), # snps_original_match_per
                round(snps_match*100/(total_snps), 2)] # snps_match_percentage
    
    # add row to results
    results.loc[len(results)] = values

    return df_snps_differents, df_snps_original_match_sample

def htz_positions(df_original, name_original):
  
    # HTZ positions in original
    df_HTZ_original = df_original[['POS']][(df_original.ALT_FREQ <= 0.85) & (df_original.ALT_FREQ >= (1 - 0.85))]
    df_HTZ_original[name_original] = ['HTZ']*len(df_HTZ_original)

    return(df_HTZ_original)

def results_sample_original_htz(out_dir, results, name_original, sample, df_dif_htz, df_match_htz):
    
    results_sample = results[results['original-sample'] == name_original+'_'+sample]
    
    htz_dif_sample = len(df_dif_htz[~df_dif_htz['HTZ'].isna()])
    htz_match_sample_1 = len(df_match_htz[~df_match_htz['HTZ'].isna()])
    
    results_sample['HTZ'] = htz_dif_sample + htz_match_sample_1

    snps_match_no_htz = len(df_match_htz[df_match_htz['HTZ'].isna()])
    snps_dif_no_htz =  len(df_dif_htz[df_dif_htz['HTZ'].isna()])
    snps_no_htz_match_perc = (snps_match_no_htz*100/(snps_match_no_htz + snps_dif_no_htz))

    results_sample['snps_no_htz_match_percentage'] = snps_no_htz_match_perc

    results_sample.to_csv(os.path.join(out_dir, 'results_%s_%s.csv' %(name_original, sample)))


def compare_results(args, output_dir, name_mix, out_seq_dir):
    
    logger.info(GREEN + "Sample: %s" %(name_mix)+ END_FORMATTING)

    out_dir = os.path.join(output_dir, name_mix)
    utils.check_create_dir(out_dir)

    # read segregated and reference sequences
    header_1, sample_1 = utils.parse_fasta('%s/%s/sample_1.fasta' %(out_seq_dir, name_mix))
    header_2, sample_2 = utils.parse_fasta('%s/%s/sample_2.fasta' %(out_seq_dir, name_mix))
    header_ref, reference = utils.parse_fasta('./data/MTB_ancestorII_reference.fa')   

    df_all = pd.DataFrame()
    df_all['POS'] = [x for x in range(1, len(reference)+1)]
    df_all['MTB_anc'] = reference
    df_all['sample_1'] = sample_1
    df_all['sample_2'] = sample_2


    # labels of results
    results = pd.DataFrame(columns=['original-sample',
                                'snps_original', 
                                'snps_sample', 
                                'total_snps', 
                                'snps_match_original', 
                                'different_alleles', 
                                'ambiguous',
                                'snps_original_match_per', 
                                'snps_match_percentage'])
    
    c = 0
    d = 0 # flags to calculate dif_htz

    vcf_to_compare = os.listdir(args.compare)

    for vcf in vcf_to_compare:

        vcf_path = os.path.join(args.compare, vcf)

        original_df_all = df_all

        df_original, name_original = add_sequence(args, original_df_all, vcf_path, output_dir)

        df_snps_differents_1, df_snps_match_1 = fill_results(results, original_df_all, df_original, name_original, 'sample_1')
        df_snps_differents_2, df_snps_match_2 = fill_results(results, original_df_all, df_original, name_original, 'sample_2')

        if name_original == name_mix.split('_')[1]:
            df_HTZ_original_1 = htz_positions(df_original, name_original)
            name_original_1 = name_original
            df_dif_1 = df_snps_differents_1
            df_match_1 = df_snps_match_1
            c += 1
        
        if name_original == name_mix.split('_')[0]:
            df_HTZ_original_2 = htz_positions(df_original, name_original)
            name_original_2 = name_original
            df_dif_2 = df_snps_differents_2
            df_match_2 = df_snps_match_2
            c += 1
        
        if c == 2 and d == 0: # see htz

            d += 1

            df_HTZ = pd.merge(df_HTZ_original_1,
                            df_HTZ_original_2,
                            on=['POS'], how='outer')

            df_both = df_HTZ[['POS']][df_HTZ[name_original_1] == df_HTZ[name_original_2]]
            df_both['HTZ'] = ['both'] * len(df_both) 

            df_HTZ_original_1 = df_HTZ_original_1[['POS']][~df_HTZ_original_1.POS.isin(df_both['POS'])]
            df_HTZ_original_1['HTZ'] = [name_original_1] * len(df_HTZ_original_1)

            df_HTZ_original_2 = df_HTZ_original_2[['POS']][~df_HTZ_original_2.POS.isin(df_both['POS'])]
            df_HTZ_original_2['HTZ'] = [name_original_2] * len(df_HTZ_original_2)

            df_HTZ = pd.concat([df_both, df_HTZ_original_1, df_HTZ_original_2])

            df_dif_1 = pd.merge(df_dif_1,
                                    df_HTZ,
                                    on='POS', how='left')
            print(len(df_match_1))
            df_match_1 = pd.merge(df_match_1,
                                    df_HTZ,
                                    on='POS', how='left')
            print(len(df_match_1))
            
            df_dif_1.to_csv(os.path.join(out_dir, 'diff_%s_sample_1.csv' %(name_original_1)))

            df_dif_2 = pd.merge(df_dif_2,
                                    df_HTZ,
                                    on='POS', how='left')
            
            df_match_2 = pd.merge(df_match_2,
                                    df_HTZ,
                                    on='POS', how='left')

            df_dif_2.to_csv(os.path.join(out_dir, 'diff_%s_sample_2.csv' %(name_original_2)))

    results = results.sort_values(by='snps_match_percentage', ascending=False)
    results.to_csv(os.path.join(out_dir, 'total_compare_%s.csv' %(name_mix)))

    results_sample_original_htz(out_dir, results, name_original_1, 'sample_1', df_dif_1, df_match_1)
    results_sample_original_htz(out_dir, results, name_original_2, 'sample_2', df_dif_2, df_match_2)
