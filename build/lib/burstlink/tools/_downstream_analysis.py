# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import zscore, ttest_ind
from sklearn.linear_model import BayesianRidge
from scipy.stats import t
from sklearn.utils import resample
import matplotlib.pyplot as plt
import seaborn as sns

from _utils import _plotting as _plotting
from  plotting import _plotting as _Plotting


def tf_tg_analysis(burst_info_file, degree_file):
    """
    Downstream analysis for gene roles in network.
    """
    burst_info_ = pd.read_csv(burst_info_file)
    burst_info = np.matrix(burst_info_)[:, 1::]
    degree_ = pd.read_csv(degree_file)
    degree = np.matrix(degree_)
    indices_tf = np.where(degree[:, 0] - degree[:, 1] < -3)[0]
    indices_tg = np.where(degree[:, 0] - degree[:, 1] > 3)[0]
    expressioninfo_tf = burst_info[indices_tf, :]
    expressioninfo_tg = burst_info[indices_tg, :]
    tf_tg_violin_plot(np.log10(expressioninfo_tf[:, 2]), np.log10(expressioninfo_tg[:, 2]))
    tf_tg_violin_plot(np.log10(expressioninfo_tf[:, 0]), np.log10(expressioninfo_tg[:, 0]))
    tf_tg_violin_plot(np.log10(expressioninfo_tf[:, 1]), np.log10(expressioninfo_tg[:, 1]))
    return

def tf_tg_violin_plot(data1, data2):
    t_stat, p_value = stats.ttest_ind(data1, data2)
    fig, ax = plt.subplots(dpi=900)
    data = [np.array(data1).flatten(), np.array(data2).flatten()]
    ax.violinplot(data, showmedians=True)
    ax.set_xticks(np.arange(1, len(data) + 1))
    ax.set_xticklabels(['tf gene', 'tg gene'])
    ax.set_ylabel('Value')
    ax.set_title('Violin Plot')  
    plt.show()
    return

def interaction_burst_regression(infer_result_filename):
    """
    Bayesian regression for gene regulatory loops and transcrptional bursting.
    """
    filtered_inference_result, unique_genename, result_unique_ = unique_infer_info(infer_result_filename, [-3.5, 2.5, -2.0, 3.0])
    # result_unique_ = [bf, bs, cv2, expression_level, overall_interaction]
    expression_level = result_unique_[:, 3]
    indices_low = np.where(np.log10(expression_level) < np.mean(np.log10(expression_level)))[0]
    indices_high = np.where(np.log10(expression_level) > np.mean(np.log10(expression_level)))[0]
    interaction_matrix_low, cmi_low = interaction_info_matrix(filtered_inference_result, unique_genename, indices_low)
    interaction_matrix_high, cmi_high = interaction_info_matrix(filtered_inference_result, unique_genename, indices_high)

    bf_results_low_ = bayes_ridge_prediction(np.log10(result_unique_[indices_low, 0]), interaction_matrix_low, cmi_low, unique_genename)
    cv2_results_low_ = bayes_ridge_prediction(np.log10(result_unique_[indices_low, 2]), interaction_matrix_low, cmi_low, unique_genename)
    bf_results_high_ = bayes_ridge_prediction(np.log10(result_unique_[indices_high, 0]), interaction_matrix_high, cmi_high, unique_genename)
    cv2_results_high_ = bayes_ridge_prediction(np.log10(result_unique_[indices_high, 2]), interaction_matrix_high, cmi_high, unique_genename)
  
    bf_results_low = bayes_ridge_prediction(result_unique_[indices_low, 0], interaction_matrix_low, cmi_low, unique_genename)
    cv2_results_low = bayes_ridge_prediction(result_unique_[indices_low, 2], interaction_matrix_low, cmi_low, unique_genename)
    bf_results_high = bayes_ridge_prediction(result_unique_[indices_high, 0], interaction_matrix_high, cmi_high, unique_genename)
    cv2_results_high = bayes_ridge_prediction(result_unique_[indices_high, 2], interaction_matrix_high, cmi_high, unique_genename)
    
    bf = np.array([bf_results_low[0][0], bf_results_low[2][0], bf_results_low[4][0], bf_results_high[0][0], bf_results_high[2][0], bf_results_high[4][0]])
    cv2 = np.array([cv2_results_low[0][0], cv2_results_low[2][0], cv2_results_low[4][0], cv2_results_high[0][0], cv2_results_high[2][0], cv2_results_high[4][0]])
    bf_sd = np.array([bf_results_low_[1][0], bf_results_low_[3][0], bf_results_low_[5][0], bf_results_high_[1][0], bf_results_high_[3][0], bf_results_high_[5][0]])
    cv2_sd = np.array([cv2_results_low_[1][0], cv2_results_low_[3][0], cv2_results_low_[5][0], cv2_results_high_[1][0], cv2_results_high_[3][0], cv2_results_high_[5][0]])
    
    _plotting.scatter_plot(np.log10(bf), np.log10(cv2), bf_sd, cv2_sd)
    return

def interaction_info_matrix(filtered_inference_result, unique_genename, indices):
    unique_genename_ = np.asarray(unique_genename)[indices]
    interaction_info_ = np.vstack([filtered_inference_result[:, np.array([0, 1, 21, 23])], filtered_inference_result[:, np.array([1, 0, 20, 22])]])
    indices = np.where(np.isin(interaction_info_[: ,0], unique_genename_))[0]
    interaction_info_ = interaction_info_[indices, :]
    interaction_info = np.matrix(sorted(interaction_info_.tolist(), key = lambda x: (x[0], x[1])))
    interaction_matrix = np.zeros([len(unique_genename_), len(unique_genename)])
    index = 0
    for genename in unique_genename_:
        indices_ = np.where(interaction_info[:, 0] == genename)[0]
        indices__ = np.where(np.isin(unique_genename, interaction_info[indices_ , 1]))[0]
        interaction_matrix[index, indices__] = np.array(interaction_info[indices_, 2].astype(float)).flatten() *  np.array(interaction_info[indices_, 3].astype(float)).flatten()
        index = index + 1
    interaction_matrix_ = np.where(interaction_matrix == 0, np.nan, interaction_matrix)
    column_cmi = np.abs(np.nanmean(interaction_matrix_, axis=0))
    return(interaction_matrix, column_cmi)
   
def bayes_ridge_prediction(y, x, w, unique_genename):
    br = BayesianRidge()
    br.fit(x, y)
    data_none = np.zeros([1, len(unique_genename)])
    data_positive = w.reshape([1, len(unique_genename)])
    data_negative = -w.reshape([1, len(unique_genename)])
    predictions1, std_dev1 = br.predict(data_positive, return_std=True)
    predictions2, std_dev2 = br.predict(data_negative, return_std=True)
    predictions3, std_dev3 = br.predict(data_none, return_std=True)
    return(predictions1, np.sqrt(std_dev1), predictions2, np.sqrt(std_dev2), predictions3, np.sqrt(std_dev3))

def unique_infer_info(infer_result_filename, thresholds):
    filtered_inference_result = _Plotting.filter_inference_results(infer_result_filename, thresholds)
    result_ = np.vstack([filtered_inference_result[:, np.array([0, 11, 13, 15, 21, 23])], filtered_inference_result[:, np.array([1, 12, 14, 16, 20, 22])]])
    sorted_result = np.matrix(sorted(result_.tolist(), key = lambda x: (x[0])))
    unique_genename = sorted(list(set(np.array(sorted_result[:, 0]).flatten().tolist())))
    result_unique_ = np.empty([len(unique_genename), 5])
    index = 0
    for genename in unique_genename:
        indices_ = np.where(sorted_result[:, 0] == genename)[0]
        if (len(indices_) == 1):
            for n in np.arange(3):
                result_unique_[index, n] = sorted_result[indices_, n+1].astype(float)
            result_unique_[index, 3] = result_unique_[index, 0] * result_unique_[index, 1]
            result_unique_[index, 4] = sorted_result[indices_, 4].astype(float) * sorted_result[indices_, 5].astype(float)
        elif (len(indices_) > 1):
            for n in np.arange(3):
                data = np.array(sorted_result[indices_, n+1].astype(float)).flatten()
                z_scores = np.abs(stats.zscore(data))
                filtered_data = np.array(data)[z_scores < 2]
                result_unique_[index, n] = np.mean(filtered_data)
            result_unique_[index, 3] = result_unique_[index, 0] * result_unique_[index, 1]
            result_unique_[index, 4] = np.sum(np.array(sorted_result[indices_, 4].astype(float)).flatten() * np.array(sorted_result[indices_, 5].astype(float)).flatten())
        index = index + 1  
    return(filtered_inference_result, unique_genename, result_unique_)

def burst_interaction_overall(burst_info_file):
    """
    Downstream analysis for gene overall regulation and transcrptional bursting.
    """
    burst_info_ = pd.read_csv(burst_info_file)
    burst_info = np.matrix(burst_info_)[:, 1::]
    expression_level = burst_info[:, 3]
    indices_low = np.where(np.log10(expression_level) < np.mean(np.log10(expression_level)))[0]
    indices_high = np.where(np.log10(expression_level) > np.mean(np.log10(expression_level)))[0]
    burst_info_low = burst_info[indices_low, :]
    burst_info_high = burst_info[indices_high, :]
    pvalue_bf_low, pvalue_bs_low, pvalue_cv2_low, pvalue_bf_high, pvalue_bs_high, pvalue_cv2_high = burst_interaction_overall_plot(burst_info_low, burst_info_high)
    return(pvalue_bf_low, pvalue_bs_low, pvalue_cv2_low, pvalue_bf_high, pvalue_bs_high, pvalue_cv2_high)

def burst_interaction_overall_plot(burstinfo1, burstinfo2):
    interaction_level1 = (burstinfo1[:, 4] - np.mean(burstinfo1[:, 4])) / np.std(burstinfo1[:, 4])
    indices_negative1 = np.where(interaction_level1 < 0)[0]
    indices_positive1 = np.where(interaction_level1 > 0)[0]
    interaction_level2 = (burstinfo2[:, 4] - np.mean(burstinfo2[:, 4])) / np.std(burstinfo2[:, 4])
    indices_negative2 = np.where(interaction_level2 < 0)[0]
    indices_positive2 = np.where(interaction_level2 > 0)[0]
    p_value1_bf, p_value1_bs, p_value1_cv2, p_value2_bf, p_value2_bs, p_value2_cv2 = _plotting.box_plot3_2(np.log10(burstinfo1[indices_positive1, 0]), np.log10(burstinfo1[indices_negative1, 0]), np.log10(burstinfo1[indices_positive1, 1]), np.log10(burstinfo1[indices_negative1, 1]), np.log10(burstinfo1[indices_positive1, 2]), np.log10(burstinfo1[indices_negative1, 2]), np.log10(burstinfo2[indices_positive2, 0]), np.log10(burstinfo2[indices_negative2, 0]), np.log10(burstinfo2[indices_positive2, 1]), np.log10(burstinfo2[indices_negative2, 1]), np.log10(burstinfo2[indices_positive2, 2]), np.log10(burstinfo2[indices_negative2, 2]))
    return(p_value1_bf, p_value1_bs, p_value1_cv2, p_value2_bf, p_value2_bs, p_value2_cv2)

def burst_info_summarize(infer_result_filename):
    inference_result__ = os.path.abspath(infer_result_filename)
    inference_result_ = pd.read_csv(inference_result__)
    inference_result = np.matrix(inference_result_)[:, 1::]
    # filter out
    indices_infer = np.where(inference_result[:, 17].astype(float) == 0)[0]
    for index in np.arange(4):
        indices_ = np.where(inference_result[:, index+4].astype(float) == 1e-4)[0]
        indices__ = np.where(inference_result[:, index+4].astype(float) == 1e4)[0]
        indices___ = list(set(indices_) | set(indices__))
        indices_infer = list(set(indices_infer) | set(indices___))
    for index in np.arange(2):
        indices_ = np.where(np.log10(inference_result[:, index+11].astype(float)) < -3.5)[0]
        indices__ = np.where(np.log10(inference_result[:, index+11].astype(float)) > 2.5)[0]
        indices___ = list(set(indices_) | set(indices__))
        indices_infer = list(set(indices_infer) | set(indices___))
    for index in np.arange(2):
        indices_ = np.where(np.log10(inference_result[:, index+13].astype(float)) < -2.0)[0]
        indices__ = np.where(np.log10(inference_result[:, index+13].astype(float)) > 3.0)[0]
        indices___ = list(set(indices_) | set(indices__))
        indices_infer = list(set(indices_infer) | set(indices___))
    filtered_inference_result = np.delete(inference_result, indices_infer, axis = 0)
    result_ = np.vstack([filtered_inference_result[:, np.array([0, 4, 6, 11, 13, 15, 21, 23])], filtered_inference_result[:, np.array([1, 5, 7, 12, 14, 16, 20, 22])]])
    sorted_result = np.matrix(sorted(result_.tolist(), key = lambda x: (x[0])))
    unique_genename = sorted(list(set(np.array(sorted_result[:, 0]).flatten().tolist())))
    result_unique_ = np.empty([len(unique_genename), 7])
    index = 0
    for genename in unique_genename:
        indices_ = np.where(sorted_result[:, 0] == genename)[0]
        if (len(indices_) == 1):
            for n in np.arange(5):
                result_unique_[index, n] = sorted_result[indices_, n+1].astype(float)
            result_unique_[index, 5] = result_unique_[index, 2] * result_unique_[index, 3]
            result_unique_[index, 6] = sorted_result[indices_, 6].astype(float) * sorted_result[indices_, 7].astype(float)
        elif (len(indices_) > 1):
            for n in np.arange(5):
                data = np.array(sorted_result[indices_, n+1].astype(float)).flatten()
                z_scores = np.abs(stats.zscore(data))
                filtered_data = np.array(data)[z_scores < 2]
                result_unique_[index, n] = np.mean(filtered_data)
            result_unique_[index, 5] = result_unique_[index, 2] * result_unique_[index, 3]
            result_unique_[index, 6] = np.sum(np.array(sorted_result[indices_, 6].astype(float)).flatten() * np.array(sorted_result[indices_, 7].astype(float)).flatten())
        index = index + 1 
    return(result_unique_)  

def affinity_burst(burst_info):
    """
    Downstream analysis for TF binding affinity and transcrptional bursting.
    """
    r = burst_info[:, 0] / burst_info[:, 1]
    bf = burst_info[:, 2]
    bs = burst_info[:, 3]
    expression_level = burst_info[:, 5]
    e, a, b, p_values_e, p_values_a, p_values_b = _plotting.fit_3d(np.log10(bf), np.log10(bs), np.log10(r))
    fig = plt.figure(figsize=(2, 1.5), dpi=300)
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(np.log10(bf), np.log10(bs), np.log10(r), c=np.log10(expression_level), s=3.5, cmap='viridis', marker='o', alpha=1)
    ax.set_xlabel('log10(bf)', fontsize=4, labelpad=-10)
    ax.set_ylabel('log10(bs)', fontsize=4, labelpad=-10)
    ax.set_zlabel('log10(r)', fontsize=4, labelpad=-10)
    ax.xaxis.label.set_position((0.5, -0.15))  
    ax.yaxis.label.set_position((-0.2, 0.5)) 
    ax.zaxis.label.set_position((-0.15, 0.5))  
    ax.tick_params(axis='x', which='major', labelsize=3, width=0.15, pad=-2) 
    ax.tick_params(axis='y', which='major', labelsize=3, width=0.15, pad=-2)  
    ax.tick_params(axis='z', which='major', labelsize=3, width=0.15, pad=-2)
    ax.w_xaxis.line.set_linewidth(0.2)
    ax.w_yaxis.line.set_linewidth(0.2)
    ax.w_zaxis.line.set_linewidth(0.2)
    plt.show()
    p_value_bf = _plotting.scatter_2d(np.log10(bf), np.log10(r), np.log10(expression_level), 'log10(bf)', 'log10(r)')
    p_value_bs = _plotting.scatter_2d(np.log10(bs), np.log10(r), np.log10(expression_level), 'log10(bs)', 'log10(r)')
    return(a, b, p_value_bf, p_value_bs)

def burst_interactionlevel_positive(readfile_name, burst_info, cv2_info):
    """
    Downstream analysis for TF binding numbers and transcrptional bursting.
    """
    degree_result_ = pd.read_csv(readfile_name)
    degree_result = np.matrix(degree_result_)
    interaction_level = zscore(burst_info[:, 6])
    indices_positive = np.where(interaction_level > 0)[0]
    indegree = np.asarray(degree_result[indices_positive, 0])
    bf = burst_info[indices_positive, 2]
    bs = burst_info[indices_positive, 3]
    expression_level = burst_info[indices_positive, 5]
    inter, a, b, p_values_inter, p_values_a, p_values_b = _plotting.fit_3d(np.log10(bf), np.log10(bs), np.log10(indegree))
    p_value_bf = _plotting.scatter_2d(np.log10(bf), np.log10(indegree), np.log10(expression_level), 'log10(bf)', 'log10(indegree)')
    p_value_bs = _plotting.scatter_2d(np.log10(bs), np.log10(indegree), np.log10(expression_level), 'log10(bs)', 'log10(indegree)')

    sorted_indegree_cv2_ = pd.read_csv(cv2_info)
    sorted_indegree_cv2 = np.array(sorted_indegree_cv2_)[:, 1::]
    x = sorted_indegree_cv2[:, 0]
    y = sorted_indegree_cv2[:, 1]
    degree = 4
    coefficients = np.polyfit(x, y, degree) 
    polynomial = np.poly1d(coefficients) 
    x_fit = np.linspace(min(x), max(x), 100)
    y_fit = polynomial(x_fit)
    n = len(y) 
    y_hat = polynomial(x) 
    residuals = y - y_hat
    mean_x = np.mean(x)
    df = n - degree - 1  
    residual_std_error = np.sqrt(np.sum(residuals**2) / df)
    t_value = t.ppf(0.95, df)
    ci = t_value * residual_std_error * np.sqrt(1/n + (x_fit - mean_x)**2 / np.sum((x - mean_x)**2))

    fig = plt.figure(figsize=(3, 1.2), dpi=300)
    ax1 = fig.add_subplot(121, projection='3d')
    scatter = ax1.scatter(np.log10(bf), np.log10(bs), np.log10(indegree), c=np.log10(expression_level), s=3, cmap='viridis', marker='o', alpha=1)
    ax1.set_xlabel('log10(bf)', fontsize=4, labelpad=-10)
    ax1.set_ylabel('log10(bs)', fontsize=4, labelpad=-10)
    ax1.set_zlabel('log10(indegree)', fontsize=4, labelpad=-10)
    ax1.xaxis.label.set_position((0.5, -0.15))  
    ax1.yaxis.label.set_position((-0.2, 0.5)) 
    ax1.zaxis.label.set_position((-0.15, 0.5))  
    ax1.tick_params(axis='x', which='major', labelsize=3, width=0.15, pad=-2) 
    ax1.tick_params(axis='y', which='major', labelsize=3, width=0.15, pad=-2)  
    ax1.tick_params(axis='z', which='major', labelsize=3, width=0.15, pad=-2)
    ax1.w_xaxis.line.set_linewidth(0.2)
    ax1.w_yaxis.line.set_linewidth(0.2)
    ax1.w_zaxis.line.set_linewidth(0.2)
    ax2 = fig.add_subplot(122) 
    ax2.scatter(x, y, s = 4)
    ax2.plot(x_fit, y_fit, linewidth=0.7)
    ax2.fill_between(x_fit, y_fit - ci, y_fit + ci, alpha=0.2)
    ax2.set_xlabel('log10(indegree)', fontsize=4)
    ax2.set_ylabel('log10(cv2)', fontsize=4)
    ax2.tick_params(axis='x', which='major', labelsize=3, width=0.08)  
    ax2.tick_params(axis='y', which='major', labelsize=3, width=0.08) 
    for spine in ax2.spines.values():
        spine.set_linewidth(0.3)
    plt.tight_layout()
    plt.show()
    t_stat_cv1, p_value_cv1 = stats.ttest_ind(sorted_indegree_cv2[0:15, 0], sorted_indegree_cv2[0:15, 1])
    t_stat_cv2, p_value_cv2 = stats.ttest_ind(sorted_indegree_cv2[14::, 0], sorted_indegree_cv2[14::, 1])
    return (a, b, p_value_bf, p_value_bs, p_value_cv1, p_value_cv2)

def dge_analysis(counts_data_file_Fibr, counts_data_file_ESC, inference_result_Fibr_, inference_result_ESC_):
    """ Function for differential gene analysis. """
    counts_data_Fibr_ = pd.read_csv(counts_data_file_Fibr)
    counts_data_ESC_ = pd.read_csv(counts_data_file_ESC)
    counts_data_Fibr = np.delete(np.matrix(counts_data_Fibr_.values.tolist()), np.s_[:1], axis=1)
    counts_data_ESC = np.delete(np.matrix(counts_data_ESC_.values.tolist()), np.s_[:1], axis=1)
    selected_countsdata_Fibr = np.zeros([1, counts_data_Fibr.shape[1]-1])
    selected_countsdata_ESC = np.zeros([1, counts_data_ESC.shape[1]-1])
    filtered_inference_result_Fibr, unique_genename_Fibr, result_unique_Fibr = unique_infer_info(inference_result_Fibr_, [-3.5, 2.5, -2.0, 3.0])
    filtered_inference_result_ESC, unique_genename_ESC, result_unique_ESC = unique_infer_info(inference_result_ESC_, [-3.5, 2.8, -2.6, 3.2])
    for n in np.arange(len(unique_genename_Fibr)):
        genename = unique_genename_Fibr[n]
        indice = np.where(np.array(unique_genename_ESC) == genename)[0]
        if (len(indice) == 1):
            index_Fibr = np.where(counts_data_Fibr[:, 0] == genename)[0][0]
            counts_data_Fibr__ = counts_data_Fibr[index_Fibr, 1::].astype(float)
            if np.isnan(counts_data_Fibr__).any():
                row_mean = np.nanmean(counts_data_Fibr__)
                counts_data_Fibr__ = np.where(np.isnan(counts_data_Fibr__), row_mean, counts_data_Fibr__)
            selected_countsdata_Fibr = np.vstack([selected_countsdata_Fibr, counts_data_Fibr__])
            index_ESC = np.where(counts_data_ESC[:, 0] == genename)[0][0]
            counts_data_ESC__ = counts_data_ESC[index_ESC, 1::].astype(float)
            if np.isnan(counts_data_ESC__).any():
                row_mean = np.nanmean(counts_data_ESC__)
                counts_data_ESC__ = np.where(np.isnan(counts_data_ESC__), row_mean, counts_data_ESC__)
            selected_countsdata_ESC = np.vstack([selected_countsdata_ESC, counts_data_ESC__])   
    group_Fibr = pd.DataFrame(selected_countsdata_Fibr, index=[f"Gene_{i+1}" for i in range(selected_countsdata_Fibr.shape[0])])
    group_ESC = pd.DataFrame(selected_countsdata_ESC, index=[f"Gene_{i+1}" for i in range(selected_countsdata_Fibr.shape[0])])
    group_ESC_upsampled = resample(group_ESC.T, replace=True, n_samples=counts_data_Fibr.shape[1]-1, random_state=123).T
    group_ESC = group_ESC_upsampled
    
    results = []
    for gene in group_Fibr.index:
        group_Fibr_values = group_Fibr.loc[gene].values
        group_ESC_values = group_ESC.loc[gene].values
        mean_group_Fibr = np.mean(group_Fibr_values)
        mean_group_ESC = np.mean(group_ESC_values)
        fc = mean_group_ESC / mean_group_Fibr
        log2_fc = np.log2(fc) if fc > 0 else 0
        t_stat, p_value = ttest_ind(group_Fibr_values, group_ESC_values, equal_var=False)
        results.append({'Gene': gene, 'log2FC': log2_fc, 'p_value': p_value})
        results_df = pd.DataFrame(results)
        results_df['-log10(p_value)'] = -np.log10(results_df['p_value'])
    log2_fc_threshold = 0.5
    p_value_threshold = 0.05
    results_df['Significance'] = 'Not Significant'
    results_df.loc[(results_df['log2FC'] > log2_fc_threshold) & (results_df['p_value'] < p_value_threshold), 'Significance'] = 'Upregulated'
    results_df.loc[(results_df['log2FC'] < -log2_fc_threshold) & (results_df['p_value'] < p_value_threshold),'Significance'] = 'Downregulated'

    plt.figure(dpi=300, figsize=(2.5, 2))
    sns.scatterplot(data=results_df, x='log2FC', y='-log10(p_value)', hue='Significance', palette={'Upregulated': '#FF6666', 'Downregulated': '#6699FF', 'Not Significant': 'grey'}, s=18, alpha=0.85)
    plt.axhline(-np.log10(p_value_threshold), color='black', linestyle='--', linewidth=0.4, label=f'p = {p_value_threshold}')
    plt.axvline(log2_fc_threshold, color='black', linestyle='--', linewidth=0.4, label=f'log2FC = ±{log2_fc_threshold}')
    plt.axvline(-log2_fc_threshold, color='black', linestyle='--', linewidth=0.4)
    plt.title('Volcano Plot', fontsize=5)
    plt.xlabel('log2(Fold Change)', fontsize=3.5)
    plt.ylabel('-log10(p-value)', fontsize=3.5)
    plt.xticks(fontsize=2.5) 
    plt.yticks(fontsize=2.5)  
    plt.tick_params(axis='both', which='both', width=0.3) 
    ax = plt.gca()
    for spine in ax.spines.values(): 
        spine.set_linewidth(0.4)
    plt.legend(title='Significance', fontsize=2, title_fontsize=2)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()
    return

def correlation_burst(inference_result_Fibr_, inference_result_ESC_):
    filtered_inference_result_Fibr, unique_genename_Fibr, result_unique_Fibr = unique_infer_info(inference_result_Fibr_, [-3.5, 2.5, -2.0, 3.0])
    filtered_inference_result_ESC, unique_genename_ESC, result_unique_ESC = unique_infer_info(inference_result_ESC_, [-3.5, 2.8, -2.6, 3.2])
    unique_results = (np.array(['bf_Fiber', 'bs_Fibr', 'bf_ES', 'bs_ES'])).reshape([1, 4])
    for n in np.arange(134):
        genename = unique_genename_Fibr[n]
        indice = np.where(np.array(unique_genename_ESC) == genename)[0]
        if (len(indice) == 1):
            result = (np.array([result_unique_Fibr[n, 0], result_unique_Fibr[n, 1], result_unique_ESC[indice, 0][0], result_unique_ESC[indice, 1][0]])).reshape([1, 4])
            unique_results = np.vstack([unique_results, result])
    unique_results = np.delete(unique_results, 62, axis=0)
    p_bf, p_bs = _plotting.scatter_interval(unique_results)
    return(p_bf, p_bs)