"""
create DE data and plot as a volcano plot

"""
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt

def get_groups(data,column):
        
    # returns the list of annotation groups from the annotation column of dataframe
    groups = list(data[column].unique())

    return groups


def make_volcano(data, gene_list, annotation_column, baseline_group, treatment_group, label_alpha=.05):
    """create a volcano plot from two sample groups.
    Parameters
    ------------
    data: pd dataframe. Must be log2-transformed counts, in tidy format (colunms->genes, rows->samples).
    gene_list: list, gene names/columns to use for DE
    annotation_column: str, the label of column in dataframe used for sample group information
    baseline_group: str, the annotation used for baseline in DE
    treatment_group: str, the annotation used for text group in DE
    label_alpha: float, the p-value threshold below which data points will be labeled with 'gene name'.
    """
    
    fold_change_dict = {}
    ttest_dict = {}

    group1 = data[data[annotation_column] == baseline_group]
    group2 = data[data[annotation_column] == treatment_group]

    genes = gene_list

    for gene in genes:

        baseline = group1[gene] 
        test = group2[gene]

        fold_change_dict[gene] = test.mean() - baseline.mean()
        p_val = -np.log10(ttest_ind(baseline, test)[1]) #get neg log of p value
        ttest_dict[gene] = p_val

        volcano_df = pd.DataFrame.from_dict(fold_change_dict, orient="Index")
        volcano_df.columns = ['fold_change']
        volcano_df['pval'] = np.fromiter(ttest_dict.values(),dtype=float)
        volcano_df.reset_index(inplace=True)

        
    ax = sns.scatterplot(x="fold_change", y="pval", data= volcano_df, s=volcano_df["pval"]*10,label="")
    for i,row in volcano_df.iterrows():
        if row['pval'] > -np.log10(label_alpha):
            ax.text(row['fold_change'],row['pval'],str(row['index']))
    ax.set(title="DE genes: {} vs. baseline {}".format(treatment_group,baseline_group), xlabel="log2 fold-change", ylabel="-log p-value")


    return ax
