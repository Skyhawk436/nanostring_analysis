#! python3

""" compile raw RCC files, normalize with HK genes, log2 tranform and QC """

import pandas as pd
import glob
import numpy as np
import os
from scipy.stats.mstats import gmean
import seaborn as sns
import matplotlib.pyplot as plt



def load_rcc(directory):

    os.chdir(directory)
    temp_folder = os.mkdir(directory + '\\temp_files')

    sample_ids = []

    for file in os.listdir():
        if 'RCC' in file:

            rcc = file
            sample_id = file.split('_')[2]
            sample_ids.append(sample_id)
            # data begins on line 23 of an RCC file
            df = pd.read_csv(file, sep=',', names=['CodeClass','Gene','Accession','Count'], header=23)
            df['SampleID'] = sample_id
            df['RCC'] = rcc
            df.to_csv(os.path.join('.\\temp_files', '{}.csv'.format(str(file))))
            

    data = pd.concat([pd.read_csv(f) for f in glob.glob('.\\temp_files\\*.csv')]).drop(columns={'Unnamed: 0'})

    print('\nRaw dataframe created:')
    print(data.head(3))

    return data


def get_annotations(directory):# create sample annotation dataframe if 'annotations' file is in the directory

    os.chdir(directory)

    for file in os.listdir():
        if 'annotations' in file:
            annotations = pd.read_csv('{}'.format(file))
            print('\nSample annotations created from {} file'.format(file))
            return annotations
        
        else:
            print('\n    Note: No annotation file was found in this directory')
    
            return None
    
def pos_qc(data):

    control_data = data[(data['CodeClass'] == 'Positive') | (data['CodeClass']=='Negative')].reset_index(drop=True)
    control_data.replace({'Count':0},1,inplace=True)
    pos_geomean = gmean(control_data[control_data.CodeClass == 'Positive']['Count'])

    background = np.mean(control_data[control_data.CodeClass == 'Negative']['Count']) + (2 * np.std(control_data[control_data.CodeClass == 'Negative']['Count']))
    print('\nNEG control background is {:.2f}\n'.format(background))

    pos_titration = {'Pos_A': 128, 'Pos_B': 32, 'Pos_C': 8, 'Pos_D': 2, 'Pos_E': 0.5, 'Pos_F': 0.125}
    
    control_data_piv = control_data.pivot(index='RCC', columns='Gene', values='Count')
    print('Control data table created..\n')

    fig, ax = plt.subplots()
    ax = sns.boxplot(x=control_data.Gene, y=np.log2(control_data.Count),order=['POS_A(128)','POS_B(32)','POS_C(8)','POS_D(2)','POS_E(0.5)','POS_F(0.125)','NEG_A(0)','NEG_B(0)','NEG_C(0)',
                                                                              'NEG_D(0)','NEG_E(0)','NEG_F(0)','NEG_G(0)'])
    #ax= sns.swarmplot(x=control_data.Gene, y=np.log2(control_data.Count),order=['POS_A(128)','POS_B(32)','POS_C(8)','POS_D(2)','POS_E(0.5)','POS_F(0.125)','NEG_A(0)','NEG_B(0)','NEG_C(0)',
    #                                                                        'NEG_D(0)','NEG_E(0)','NEG_F(0)','NEG_G(0)'], color='b')
    ax.axhline(y=np.log2(background), linestyle='--',color='r', label='background')
    ax.set_title('Control counts')
    ax.legend()
    

def endog_qc(data):

    # drop the 'message' rows at the bottom, and keep only Endog. and HK genes for analysis
    data2 = data[(data.CodeClass == 'Endogenous') | (data.CodeClass == 'Housekeeping')].reset_index()

    #replace 0 counts with 1 (to avoid annoying divide by 0 issues with log-transformed counts)
    data2.replace({'Count':0},1,inplace=True)

    #plot histogram of raw gene counts
    fig2, ax2 = plt.subplots()
    ax2 = sns.distplot(np.log2(data2.Count), label='raw count distribution')
    ax2.legend()
    ax2.set_xlabel('log2 counts')
    
    
    hk_genes = list(set(data2.query("CodeClass=='Housekeeping'")['Gene']))
    
    #transform the data table to 'tidy' form for analysis
    df_piv = data2.pivot(index='RCC',columns='Gene',values='Count')

    # add columns for HK geomean, HK normalization factors and then plot normalization factors
    df_piv['hk_gmean'] = df_piv[hk_genes].apply(lambda x: gmean(x), axis=1)
    avg_hk = np.mean(df_piv['hk_gmean'])
    df_piv['hk_normfactor'] = df_piv['hk_gmean'] / avg_hk

    #line plot of HK normalization factors. Should be near 1 if equal loading of total RNA in hybridizations

    fig3, ax3 = plt.subplots()
    ax3 = sns.lineplot(data=df_piv['hk_normfactor'], label='ref gene normalization factor')
    ax3.set_ylim([0,3])
    ax3.axhline(y=1, linestyle='--',color='r', label='reference line=1')
    ax3.set_xlabel('')
    ax3.set_xticks('')
    ax3.set_ylabel('normalization factor')
    plt.legend()

    return df_piv
    

def hk_normalize(data, annotations, save_directory):

    #multiply raw gene counts by the housekeeper normalization factor, then log2 transform

    normdata = data.loc[:,'ABCC8':'ZNRF2'].apply(lambda x: x * data['hk_normfactor'],axis=0)

    log2_norm = np.log2(normdata).reset_index()

    #merge dataframe with sample annotation file if it exists:
    if annotations:

        log2_norm_annot = log2_norm.merge(annotations, left_on="RCC", right_on = "RCC",how='left')

        log2_norm_annot.to_csv(save_directory + '/log2_normalized_data.csv')
        print(log2_norm_annot.head(3))
        
    else:
        log2_norm.to_csv(save_directory + '/log2_normalized_data.csv')
        print(log2_norm.head(3))

    print('\nLog2-HK-normalized data saved to: {}'.format(save_directory))
    print('\n')

##    groups = list(log2_norm_annot['sample_id'].unique())
##    colors = sns.xkcd_palette(["windows blue", "amber", "greyish", "faded green", "dusty purple"])
##    group_colors = colors[:len(groups)]
##    sns.clustermap(log2_norm_annot.set_index('RCC').T, z_score=0, cmap='bwr',col_colors = group_colors, figsize=(10,8))


def process_rcc_data():
    
    print('\n\n')
    print("--------------------- Load and analyse Nanostring RCC data ------------------------")
    print('\n')

    directory = input(r'Enter the path to RCC directory: ')

    annotations = get_annotations(directory)
    data = load_rcc(directory)
    pos_qc(data)
    data2 = endog_qc(data)
    hk_normalize(data2, annotations, directory)

    plt.show()
    
    
    
          
if __name__ == '__main__':

    process_rcc_data()
        
