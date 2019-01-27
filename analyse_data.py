################################################################################
#   Script to analyse microbiome data
#
#   input: combined_data and OTUID_taxonomy_dic
#
#   structure:
#       1.filtering & preprocessing data
#       2.abundance analysis
#           - statistical tests (wilcoxon & mann-whitney u)
#           - fold change
#       3.shannon (alpha) diversity
#
################################################################################

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu # for Mann-Whiteny U test
from scipy.stats import wilcoxon # For one sampke Wilcoxon signed rank test
# shouldnt import like this
from helper_functions import *

# params
#-------------------------------------------------------------------------------
PLOT = False
# fold change, influence parameters, div for all data, div for each phylum,
# zero entries analysis, volcano plot
plot_list = [0,0,0,0,0,1]
save_plots = False
alpha = 0.05
# false discovery rate
fdr = 0.05

if __name__ == '__main__':
    # load data
    data = pd.read_pickle("01 temp data/combined_data")
    dic = pd.read_pickle("01 temp data/OTUID_taxonomy_dic")
    #---------------------------------------------------------------------------
    # 1.filtering & preprocessing data
    #---------------------------------------------------------------------------

    # delete rows that dont have consent_age value, location
    data.dropna(subset = ["consent_age"],inplace = True)
    data.dropna(subset = ["biopsy_location"],inplace = True)

    # Drop columns which have all values as zero
    #data= data.loc[:, (data != 0).any(axis=0)]

    # normalize data
    norm_data = normalize_data_samples(data,data.columns.get_loc("IP8BSoli"))
    # for OTU ID: ["#OTU ID"]
    header_norm_data = dic["#OTU ID"]
    #norm_data.to_csv("01 temp data/nomalized_data.csv")

    # collapsing the data
    data_phylum, header_phylum = collapse_data("phylum",norm_data,dic)

    # get UC,CD and nonIBD data
    data_UC, data_CD, data_cont = return_subsets(data_phylum)
    data_UC_uncol, data_CD_uncol, data_cont_uncol = return_subsets(norm_data)

    #---------------------------------------------------------------------------
    # 2.abundance analysis
    #---------------------------------------------------------------------------

    # phylum level
    #---------------------------------------------------------------------------

    # split data for wilcoxon and mann-whitney
    data_cont_MW, data_UC_MW, data_UC_WC, data_CD_MW, data_CD_WC = \
    split_for_stat_test(data_UC, data_CD, data_cont)

    # perform mann-whitney u test for UC <=> control
    stat_list_uc_MW, p_list_uc_MW, ind_diff_uc_MW = multi_mann_whitney(alpha,data_cont_MW, \
    data_UC_MW )

    # benjamin hochberg correction
    sig_ind_uc_MW = ben_hoch_corr(p_list_uc_MW, fdr)
    # nothing is significant!

    # perform mann-whitney u test for CD <=> control
    stat_list_cd_MW, p_list_cd_MW, ind_diff_cd_MW = multi_mann_whitney(alpha,data_cont_MW, \
    data_CD_MW )

    # benjamin hochberg correction
    sig_ind_cd_MW = ben_hoch_corr(p_list_cd_MW, fdr)
    p_list_cd_MW = np.array(p_list_cd_MW)
    # significant differences for phyla cont vs. cd
    print("Significantly different abundance after correction phylum level (CD vs. control):")
    print(header_phylum[sig_ind_cd_MW])
    print(p_list_cd_MW[sig_ind_cd_MW])

    # perform wilcoxon test for UC
    stat_list_uc_WC, p_list_uc_WC, ind_diff_uc_WC = multi_wilcoxon_test(alpha, \
    data_UC_WC)

    # perform wilcoxon test for CD
    stat_list_cd_WC, p_list_cd_WC, ind_diff_cd_WC = multi_wilcoxon_test(alpha, \
    data_CD_WC)

    # calculate fold change using data for mann-whitney test (values are not
    # all zero in the healthy subjects) --> only consider variables that showed
    # a significant difference
    fold_change_list_UC, fold_change_list_log_UC, h = \
    cal_fold_change(data_cont_MW, data_UC_MW)

    fold_change_list_CD, fold_change_list_log_CD, h = \
    cal_fold_change(data_cont_MW, data_CD_MW)

    # fold change on phylum level
    fold_change_list_phyl_CD, fold_change_list_phyl_log_CD, header_phyl_CD = \
    cal_fold_change(data_cont, data_CD)

    fold_change_list_phyl_UC, fold_change_list_phyl_log_UC, header_phyl_UC = \
    cal_fold_change(data_cont, data_UC)

    #print(np.array(fold_change_list_phyl_log_CD).shape)

    # species level
    #---------------------------------------------------------------------------

    # split data for wilcoxon and mann-whitney
    data_cont_unc_MW, data_UC_unc_MW, data_UC_unc_WC, data_CD_unc_MW, data_CD_unc_WC = \
    split_for_stat_test(data_UC_uncol, data_CD_uncol, data_cont_uncol)

    # UC vs. control

    # perform mann-whitney u test for UC <=> control
    stat_list_uc_unc_MW, p_list_uc_unc_MW, ind_diff_uc_unc_MW = multi_mann_whitney(alpha,data_cont_unc_MW, \
    data_UC_unc_MW)

    # benjamin hochberg correction
    sig_ind_uc_MW = ben_hoch_corr(p_list_uc_unc_MW, fdr)
    p_list_uc_MW = np.array(p_list_uc_unc_MW)
    # significant differences for phyla cont vs. cd
    print("Significantly different abundance after correction (UC vs. control):")
    print(header_norm_data[sig_ind_uc_MW])
    print(p_list_uc_MW[sig_ind_uc_MW])

    fold_change_list_unc_uc, fold_change_list_unc_uc_log, h = \
    cal_fold_change(data_cont_uncol[header_norm_data[sig_ind_uc_MW]], data_CD_uncol[header_norm_data[sig_ind_uc_MW]])
    print(fold_change_list_unc_uc_log)

    # CD vs. control

    # perform mann-whitney u test for UC <=> control
    stat_list_cd_unc_MW, p_list_cd_unc_MW, ind_diff_cd_unc_MW = multi_mann_whitney(alpha,data_cont_unc_MW, \
    data_CD_unc_MW)

    # benjamin hochberg correction
    sig_ind_cd_MW = ben_hoch_corr(p_list_cd_unc_MW, fdr)
    p_list_cd_MW = np.array(p_list_cd_unc_MW)
    # significant differences for phyla cont vs. cd
    print("Significantly different abundance after correction (CD vs. control):")
    print(header_norm_data[sig_ind_cd_MW])
    print(p_list_cd_MW[sig_ind_cd_MW])

    fold_change_list_unc_cd, fold_change_list_unc_cd_log, h = \
    cal_fold_change(data_cont_uncol[header_norm_data[sig_ind_cd_MW]], data_CD_uncol[header_norm_data[sig_ind_cd_MW]])
    print(fold_change_list_unc_cd_log)

    #---------------------------------------------------------------------------
    # 3.alpha diversity
    #---------------------------------------------------------------------------

    # species level
    #---------------------------------------------------------------------------
    # use uncollapsed data (species level)
    sd_cont_spec = calc_shannon_div(data_cont_uncol)
    sd_UC_spec = calc_shannon_div(data_UC_uncol)
    sd_CD_spec = calc_shannon_div(data_CD_uncol)

    # statistical significance using mann-whitney
    s_u, p_u = mannwhitneyu(sd_cont_spec,sd_UC_spec)
    s_c, p_c = mannwhitneyu(sd_cont_spec, sd_CD_spec)
    s, p = mannwhitneyu(sd_UC_spec, sd_CD_spec)

    #print(p_u, p_c, p)

    # species level for each phylum
    #---------------------------------------------------------------------------
    unique_el, otu_arr = data_per_phylum("phylum",norm_data,dic)
    div_phyl = []

    for i,el in enumerate(unique_el):
        sd_cont = calc_shannon_div(data_cont_uncol[otu_arr[i]])
        sd_UC = calc_shannon_div(data_UC_uncol[otu_arr[i]])
        sd_CD = calc_shannon_div(data_CD_uncol[otu_arr[i]])
        # statistical significance using mann-whitney

        labels = ["control","UC","CD"]
        fig, ax = plt.subplots()
        plt.boxplot([sd_cont, sd_UC\
        ,sd_CD], labels = labels)
        plt.title(unique_el[i])
        plt.savefig("02 Plots/Phylum_level_alpha/"+"plot"+unique_el[i]+".svg")
        try:
            s_u, p_u = mannwhitneyu(sd_cont,sd_UC)
            s_c, p_c = mannwhitneyu(sd_cont, sd_CD)
            s, p = mannwhitneyu(sd_UC, sd_CD)
            to_add = [unique_el[i],p_u,p_c,p]
            div_phyl.append(to_add)
        except:
            continue
    div_phyl = np.array(div_phyl)

    # non zero elements analysis
    #---------------------------------------------------------------------------
    nr_nz_bac = nr_zero_entries(data,data.columns.get_loc("IP8BSoli"),data.shape[1])

    if PLOT:
        if plot_list[0]:
            # plot fold change
            #-------------------------------------------------------------------
            fold_change_list_phyl_UC = np.array(fold_change_list_phyl_UC)

            bars = (header_phyl_UC)
            y_pos = np.arange(fold_change_list_phyl_UC.shape[0])

            plt.barh(y_pos, fold_change_list_phyl_UC)
            plt.yticks(y_pos, bars)
            plt.title("Fold change (phylum level): UC vs. control")
            plt.xlabel("Fold change")
            plt.ylabel("Phylum name")
            #plt.xlim(np.min(fold_change_list_phyl_UC)-20,\
            #np.max(fold_change_list_phyl_UC)+20)
            if save_plots:
                plt.savefig("02 Plots/fold_change_phylum_UC_vs_control.svg")
            plt.show()



            fold_change_list_phyl_CD = np.array(fold_change_list_phyl_CD)
            bars = (header_phyl_CD)
            y_pos = np.arange(fold_change_list_phyl_CD.shape[0])
            plt.barh(y_pos, fold_change_list_phyl_CD, color= "orange")
            plt.yticks(y_pos, bars)
            plt.title("Fold change (phylum level): CD vs. control")
            plt.xlabel("Fold change")
            plt.ylabel("Phylum name")
            plt.xlim(-10,10)
            #plt.xlim(np.min(fold_change_list_phyl_CD)-200,\
            #np.max(fold_change_list_phyl_CD)+50)
            if save_plots:
                plt.savefig("02 Plots/fold_change_phylum_CD_vs_control.svg")
            plt.show()


        if plot_list[1]:
            # differences between uc, healthy and cd for age, gender, biopsy location
            #-----------------------------------------------------------------------
            plot_hist(data_UC,data_CD,data_cont,"consent_age")
            #plot_hist(data_UC,data_CD,data_cont,"biopsy_location")
            labels = ["control","UC","CD"]
            fig, ax = plt.subplots()
            plt.boxplot([data_UC["consent_age"], data_CD["consent_age"]\
            ,data_cont["consent_age"]], labels = labels)
            ax.set_ylabel("consent age")
            ax.set_title("Age")
            if save_plots:
                plt.savefig("02 Plots/age_boxplot.svg")
            plt.show()

        if plot_list[2]:
            # alpha diversity species level for all data
            #-------------------------------------------------------------------
            labels = ["control","UC","CD"]
            fig, ax = plt.subplots()
            plt.boxplot([sd_cont_spec, sd_UC_spec, sd_CD_spec], labels = labels)
            ax.set_ylabel("Shannon index of biodiversity")
            ax.set_title("Alpha diversity for entire data")
            # add statistical annotation
            x1, x2 = 1, 2   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
            y, h, col = np.max(sd_cont_spec) + 0.1, 0.1, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y+h, "ns", ha='center', va='bottom', color=col)
            x1, x2 = 1, 3   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
            y, h, col = np.max(sd_cont_spec) + 0.3, 0.1, 'k'
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y+h, "*", ha='center', va='bottom', color=col)
            if save_plots:
                plt.savefig("02 Plots/alpha_diversity.svg")
            plt.show()

        if plot_list[3]:
            # alpha diversity for each phylum
            #-------------------------------------------------------------------
            # plot heatmap of p-values for comparisons of alpha diversity for each
            # phylum
            phyla_list = div_phyl[:,0]
            comp_list = ["cont vs. UC", "cont vs. CD", "UC vs. CD"]
            dat = div_phyl[:,1:].astype(float)
            ax = plt.axes()
            sns.heatmap(dat, annot=True,  \
            ax = ax, linewidths=.5,xticklabels=comp_list,yticklabels = phyla_list,cbar_kws={'label': 'p-value'})
            ax.set_title('Comparison of alpha diversity for each phylum')
            if save_plots:
                plt.savefig("02 Plots/alpha_diversity_phylum.svg")
            plt.show()
            #plt.imshow(div_phyl[:,1:], cmap='hot', interpolation='nearest')
            #plt.show()

        if plot_list[4]:
            # non zero elements analysis
            #-------------------------------------------------------------------
            plt.subplot(1,2,1)
            plt.hist(nr_nz_bac, bins=range(min(nr_nz_bac), max(nr_nz_bac) + 1, 1))
            plt.xlabel("#samples")
            plt.ylabel("#present bacteria")
            plt.subplot(1,2,2)
            plt.scatter(np.arange(nr_nz_bac.shape[0]),nr_nz_bac)
            plt.xlabel("bacteria ID")
            plt.ylabel("#samples where bacterium is present")
            if save_plots:
                plt.savefig("02 Plots/zero_element_analysis.svg")
            plt.show()

        if plot_list[5]:
            # volcano plot
            #-------------------------------------------------------------------
            # cd
            p_values_cd = np.absolute(np.log10(p_list_cd_MW[sig_ind_cd_MW]))
            fold_change_cd = fold_change_list_unc_cd_log
            print(fold_change_cd)
            plt.scatter(fold_change_cd,p_values_cd)
            plt.xlabel("log2(fold change)")
            plt.ylabel("-log10(p-value)")
            plt.title("Volcano plot (species level): CD vs. control")
            for i, txt in enumerate(header_norm_data[sig_ind_cd_MW]):
                if np.absolute(fold_change_cd[i]) > 3.5 or p_values_cd[i] > 6:
                    plt.annotate(txt, (fold_change_cd[i], p_values_cd[i]+0.05))
            if save_plots:
                plt.savefig("02 Plots/volcano_CD_vs_control.svg")
            plt.show()


            p_values_uc = np.absolute(np.log10(p_list_uc_MW[sig_ind_uc_MW]))
            fold_change_uc = fold_change_list_unc_uc_log
            plt.scatter(fold_change_uc,p_values_uc)
            plt.xlabel("log2(fold change)")
            plt.ylabel("-log10(p-value)")
            plt.title("Volcano plot (species level): UC vs. control")
            for i, txt in enumerate(header_norm_data[sig_ind_uc_MW]):
                plt.annotate(txt, (fold_change_uc[i], p_values_uc[i]+0.005))
            if save_plots:
                plt.savefig("02 Plots/volcano_UC_vs_control.svg")
            plt.show()
