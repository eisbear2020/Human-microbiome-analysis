################################################################################
#   Contains all helper functions for analyse_data.py
#
#
#
#
#
#
################################################################################
from scipy.stats import mannwhitneyu # for Mann-Whiteny U test
from scipy.stats import wilcoxon # For one sampke Wilcoxon signed rank test
import matplotlib.pyplot as plt
from math import log as ln
import pandas as pd
import numpy as np

def plot_hist(d1,d2,d3,col_name):
    plt.subplot(1,3,1)
    plt.hist(d1[col_name],bins = 30)
    plt.title("UC")
    plt.subplot(1,3,2)
    plt.hist(d2[col_name],bins = 30)
    plt.title("CD")
    plt.subplot(1,3,3)
    plt.hist(d3[col_name],bins = 30)
    plt.title("control")
    plt.suptitle(col_name)
    plt.show()
    pass

def save_subsets(comb_data):
    # separate UC,CD and nonIBD and save as pickle
    data_UC = comb_data.loc[comb_data["diagnosis"] == "UC"]
    data_CD = comb_data.loc[comb_data["diagnosis"] == "CD"]
    data_cont = comb_data.loc[comb_data["diagnosis"] == "nonIBD"]

    data_UC.to_pickle("01 temp data/data_UC")
    data_CD.to_pickle("01 temp data/data_CD")
    data_cont.to_pickle("01 temp data/data_cont")

def return_subsets(comb_data):
    # separate UC,CD and nonIBD and save as pickle
    data_UC = comb_data.loc[comb_data["diagnosis"] == "UC"]
    data_CD = comb_data.loc[comb_data["diagnosis"] == "CD"]
    data_cont = comb_data.loc[comb_data["diagnosis"] == "nonIBD"]
    return data_UC, data_CD, data_cont


def load_subsets():
    # load UC,CD and nonIBD data
    data_UC = pd.read_pickle("01 temp data/data_UC")
    data_CD = pd.read_pickle("01 temp data/data_CD")
    data_cont = pd.read_pickle("01 temp data/data_cont")
    return data_UC, data_CD, data_cont

def nr_zero_entries(comb_data,c_start,c_end):
    # non zero values in each column of dataframe between column c_start and
    # c_end
    temp = comb_data.copy()
    nr_nz = temp.fillna(0).astype(bool).sum(axis=0)

    return nr_nz[c_start:(c_end+1)]

def normalize_data_samples(data, c_start):
    # normalize values in columns between c_start and c_end by maximum value
    # normalization is done for each row separately
    non_bac_col = data.iloc[:,:c_start]
    col_to_norm = data.iloc[:,c_start:]

    col_to_norm= col_to_norm.div(col_to_norm.max(axis=1), axis=0)

    return pd.concat([non_bac_col,col_to_norm],axis = 1)



def collapse_data(sel_level,data,dic):
    # collapse data:
    # "life","phylum","class","order","family","genus"

    # split taxonomy column in several columns (corresponding to bact. taxonomy)
    temp_dic = dic.copy()
    split = temp_dic["taxonomy"].str.split(";",expand=True)
    # remove taxonomy column
    temp_dic.drop(columns = "taxonomy", inplace = True)
    dic_tax = pd.concat([temp_dic,split], axis = 1)
    dic_tax.rename(index=str, columns={0:"life", 1:"phylum",2:"class",3:"order",4:"family",5:"genus"}, inplace = True)

    # find unique elements on selected level
    unique_el = dic_tax[sel_level].unique()

    # create data frame for merged data
    merged_data = pd.DataFrame(columns = unique_el)

    # go through unique elements
    for id in unique_el:
        indices_matched = dic_tax.loc[dic_tax[sel_level] == id,"#OTU ID"]
        merged_data[id] = data[indices_matched].sum(axis="columns")

    # merge metadata and collapsed data
    data_col_other = data.iloc[:,:data.columns.get_loc("IP8BSoli")]
    #print(data_col_other)

    return pd.concat([data_col_other,merged_data],axis = 1), unique_el

def data_per_phylum(sel_level,data,dic):
    # split taxonomy column in several columns (corresponding to bact. taxonomy)
    temp_dic = dic.copy()
    split = temp_dic["taxonomy"].str.split(";",expand=True)
    # remove taxonomy column
    temp_dic.drop(columns = "taxonomy", inplace = True)
    dic_tax = pd.concat([temp_dic,split], axis = 1)
    dic_tax.rename(index=str, columns={0:"life", 1:"phylum",2:"class",3:"order",4:"family",5:"genus"}, inplace = True)

    # find unique elements on selected level
    unique_el = dic_tax[sel_level].unique()

    # OTU array
    otu_arr = []

    # go through unique elements
    for id in unique_el:
        indices_matched = dic_tax.loc[dic_tax[sel_level] == id,"#OTU ID"]
        #print("ID:",id, "OTU ID: ",indices_matched.iloc[:])
        otu_arr.append(np.array(indices_matched))

    return unique_el, otu_arr


    # merge metadata and collapsed data
    #data_col_other = data.iloc[:,:data.columns.get_loc("IP8BSoli")]
    #print(data_col_other)

    #return pd.concat([data_col_other,merged_data],axis = 1), unique_el






def split_for_stat_test(data_UC, data_CD, data_cont):
    # splits data into two subsets using control subset:
    # 1) Mann-whitney U test: variables for which all values are zero in healthy
    #    subjects
    # 2) Wilcoxon test: all other variables
    # returns
    #---------------------------------------------------------------------------

    # only use taxonomic data and exclude metadata
    data_UC = data_UC.iloc[:,data_UC.columns.get_loc("Tube A and B received at Broad:")+1:]
    data_CD = data_CD.iloc[:,data_CD.columns.get_loc("Tube A and B received at Broad:")+1:]
    data_cont = data_cont.iloc[:,data_cont.columns.get_loc("Tube A and B received at Broad:")+1:]

    # dataset containing variables having only zero values
    data_cont_WC = data_cont.loc[:, (data_cont == 0).all(axis=0)]
    # data for variables which are all not zero
    data_cont_MW = data_cont.loc[:, (data_cont != 0).any(axis=0)]

    #Filtering UC and CD based on the above phylums:
    #variables in which "healthy" values were 0
    data_UC_WC = data_UC[data_cont_WC.columns.values]
    #variables in which all "healthy" values were NOT 0
    data_UC_MW = data_UC.drop(columns= data_cont_WC.columns.values)

    data_CD_WC = data_CD[data_cont_WC.columns.values]
    data_CD_MW = data_CD.drop(columns= data_cont_WC.columns.values)

    return data_cont_MW, data_UC_MW, data_UC_WC, data_CD_MW, data_CD_WC

def multi_mann_whitney(alpha, d1, d2):
    # computes mann-whitney u test for two samples d1 and d2
    stat_list=[]
    #p values
    p_list=[]
    # gives the index of columns for the significantly different variables
    ind_diff=[]

    if d1.shape[1] != d2.shape[1]:
        raise ValueError("data sets of unequal length")

    for i in range(d1.shape[1]):
        s1 =d1.iloc[:,i]
        s2 = d2.iloc[:,i]

        stat,p = mannwhitneyu(s1,s2)
        stat_list.append(stat)
        p_list.append(p)
        if p < alpha:
            ind_diff.append(i)
    return stat_list, p_list, ind_diff

def multi_wilcoxon_test(alpha, d1):
    # computes wilcoxon one sample test
    stat_list=[]
    #p values
    p_list=[]
    # gives the index of columns for the significantly different variables
    ind_diff=[]

    # delete all variables with all zero values
    d1_temp = d1.copy()
    d1_temp = d1_temp.loc[:, (d1_temp != 0).any(axis=0)]

    for i in range(d1_temp.shape[1]):
        s1 =d1_temp.iloc[:,i]
        stat,p = wilcoxon(s1)
        stat_list.append(stat)
        p_list.append(p)
        if p < alpha:
            ind_diff.append(i)
    return stat_list, p_list, ind_diff

def cal_fold_change(d1_ref,d2_in):
    # compute fold change between to data sets d1_ref and d2
    # outputs log and non-transformed values
    # d1_ref is the reference data and d2_in is related to the reference
    d1 = d1_ref.copy()
    d2 = d2_in.copy()

    if "Tube A and B received at Broad:" in d1.columns:
        # only use taxonomic data
        d1 = d1.iloc[:,d1.columns.get_loc("Tube A and B received at Broad:")+1:]

    if "Tube A and B received at Broad:" in d2.columns:
        # only use taxonomic data
        d2 = d2.iloc[:,d2.columns.get_loc("Tube A and B received at Broad:")+1:]

    fold_change_list =[]
    fold_change_list_log =[]

    d1_mean = d1.mean().to_frame().T
    d2_mean = d2.mean().to_frame().T
    # add small value to elements for log
    d1_mean_log = d1_mean #+ 0.0000001
    d2_mean_log = d2_mean #+ 0.0000001

    for i in range(d1.shape[1]):
        # for log transform
        if d1_mean.iloc[0,i] == 0 or d2_mean.iloc[0,i] == 0:
            # append 0 if one of the values is zero
            fold_change_list_log.append(0)
        else:
            fold_change_list_log.append(ln(d2_mean_log.iloc[0,i],2.0)-ln(d1_mean_log.iloc[0,i],2.0))
        # for non log transform
        if d1_mean.iloc[0,i] == 0 or d2_mean.iloc[0,i] == 0:
            fold_change_list.append(0)
        elif d1_mean.iloc[0,i]/d2_mean.iloc[0,i] >= 1:
            # ratio
            fold_change_list.append(-(d1_mean.iloc[0,i]/d2_mean.iloc[0,i]))
        else:
            fold_change_list.append((d2_mean.iloc[0,i]/d1_mean.iloc[0,i]))

    return fold_change_list, fold_change_list_log, d1.columns.values

def calc_shannon_div(d):
    # compute shannon diversity of data set d
    shannon_div = []
    if "Tube A and B received at Broad:" in d.columns:
        # only use taxonomic data
        d_tax = d.iloc[:,d.columns.get_loc("Tube A and B received at Broad:")+1:]
    else:
        d_tax = d
    # add small value for log transformation
    d_tax = d_tax + 0.0000001

    #Finding proportion:
    for i in range(d_tax.shape[0]):
        d_tax.iloc[i,:]/= d_tax.iloc[i,:].sum()

    # calculate sum of terms
    for i in range(d_tax.shape[0]):
        a=0
        for j in range(d_tax.shape[1]):
            a += (d_tax.iloc[i,j]*(ln(d_tax.iloc[i,j])))
        shannon_div.append(a)
    shannon_div = np.absolute(shannon_div)

    return shannon_div

def ben_hoch_corr(p_list, fdr):
    # benjamin hochberg correction
    # fdr: false discovery rate

    # remember indices
    sort_ind = np.argsort(p_list)
    p_list_sorted = sorted(p_list)
    c_list=[]
    m = len(p_list_sorted)
    for i in range(len(p_list_sorted)):
        c_list.append(((i+1)/m)*fdr)
    # list with p-values and corrected values
    sig_list =list(zip(p_list_sorted,c_list))
    # only rows where the corrected value is smaller than the initial p-value the
    # result is significant
    sig_ind = []
    for i in range(len(p_list_sorted)):
        if p_list_sorted[i] < c_list[i]:
            sig_ind.append(sort_ind[i])
    return sig_ind
