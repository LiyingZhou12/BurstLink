# -*- coding: utf-8 -*-
import os, math
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import pearsonr
from scipy.ndimage import gaussian_filter1d, gaussian_filter
import statsmodels.api as sm
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

from _utils import _synthetic_data
from tools import _burst_interactions as bi


def landscape(params, fig_index):
    # params = [kon1, kon2, koff1, koff2, ksyn1, ksyn2, kdeg, rK, h1, h2, eps_, T]
    S_stable, mkonx, mkony, mbf1, mbf2, mbs1, mbs2, mbf1_, mbf2_, mbs1_, mbs2_ = _synthetic_data.SSA_coexpression(params, verbose = 'burst', fig = False)
    simul_data = S_stable.astype(int)
    pxy = bi.density_estimation(simul_data)
    geneinfo = np.array(['X1', 'X2', -math.copysign(1, params[8]), -math.copysign(1, params[9])]).reshape([1, 4])
    infer_result = bi.genepair_inference(simul_data, geneinfo, 0, verbose1 = False, verbose2 = False, verbose3 = 0, test = False)
    stat_params = infer_result[0: 7]
    est_pxy, est_px, est_py = bi.Binary_PoBe(stat_params, simul_data)
    fig, ax = plt.subplots(dpi=900)
    g = sns.jointplot(x=simul_data[0, :], y=simul_data[1, :], kind="kde", marginal_kws={'fill': False, 'alpha': 0.0}, fill=True, cut=0)
    px = ax.hist(simul_data[0, :], bins = np.max(simul_data[0, :]), density=True, alpha=0)[0]
    py = ax.hist(simul_data[1, :], bins = np.max(simul_data[1, :]), density=True, alpha=0)[0]
    px = gaussian_filter1d(px, sigma=0.5)
    py = gaussian_filter1d(py, sigma=0.5)
    x_range = np.arange(len(px))
    y_range = np.arange(len(py))
    g.ax_marg_x.plot(x_range, px)
    g.ax_marg_y.plot(py, y_range)
    g.ax_marg_x.fill_between(x_range, px, color='#257AB6', alpha=0.15)  
    g.ax_marg_y.fill_betweenx(y_range, py, color='#257AB6', alpha=0.15)
    g.ax_marg_x.scatter(np.arange(np.max(simul_data[0, :])), est_px, facecolors='none', edgecolors='#257AB6', s=80, linewidths=2, alpha=0.65)
    g.ax_marg_y.scatter(est_py, np.arange(np.max(simul_data[1, :])), facecolors='none', edgecolors='#257AB6', s=80, linewidths=2, alpha=0.65)
    g.ax_joint.set_xlim([0, np.max(simul_data[0, :])])
    g.ax_joint.set_ylim([0, np.max(simul_data[1, :])])
    g.ax_marg_x.set_xlim([0, np.max(simul_data[0, :])])
    g.ax_marg_y.set_ylim([0, np.max(simul_data[1, :])])
    plt.show()
    fig = plt.figure(dpi=900)
    ax = fig.add_subplot(projection='3d')
    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    xx = np.arange(np.max(simul_data[0, :]))
    yy = np.arange(np.max(simul_data[1, :]))
    x, y = np.meshgrid(xx, yy)
    z = pxy
    ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=plt.cm.coolwarm)
    plt.xticks([])
    plt.yticks([])
    ax.set_zticks([])
    plt.show()
    fig = plt.figure(dpi=900)
    ax = fig.add_subplot(projection='3d')
    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    xx = np.arange(np.max(simul_data[0, :]))
    yy = np.arange(np.max(simul_data[1, :]))
    x, y = np.meshgrid(xx, yy)
    z = gaussian_filter(pxy, sigma=2)
    ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=plt.cm.plasma, alpha=0.95)
    plt.xticks([])
    plt.yticks([])
    ax.set_zticks([])  
    plt.show()
    return

def ROC_plot(readfile_name):
    filtered_inference_result = filter_inference_results(readfile_name, [-0.5, 1.5, 0.0, 1.0])
    interaction_grn = np.asarray(np.vstack([filtered_inference_result[:, 2].astype(float), filtered_inference_result[:, 3].astype(float)]))
    interaction_infer = np.asarray(np.vstack([np.abs(filtered_inference_result[:, 22].astype(float)), np.abs(filtered_inference_result[:, 23].astype(float))]))
    fpr, tpr, thresholds = roc_curve(interaction_grn, interaction_infer)
    auc_score = roc_auc_score(interaction_grn, interaction_infer)
    # Create ROC plot
    fig, ax = plt.subplots(dpi=900)
    plt.plot(fpr, tpr, label='ROC Curve (AUC = {:.2f})'.format(auc_score), c = '#4292c5')
    plt.plot([0, 1], [0, 1], 'k--', label='Random', c = '#9dc4e6')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend()
    plt.show()
    return

def filter_inference_results(readfile_name, thresholds):
    inference_result_ = pd.read_csv(readfile_name)
    inference_result = np.matrix(inference_result_)[:, 1::]
    # filter out
    indices_infer = np.where(inference_result[:, 17].astype(float) == 0)[0]
    for index in np.arange(4):
        indices_ = np.where(inference_result[:, index+4].astype(float) == 1e-4)[0]
        indices__ = np.where(inference_result[:, index+4].astype(float) == 1e4)[0]
        indices___ = list(set(indices_) | set(indices__))
        indices_infer = list(set(indices_infer) | set(indices___))
    for index in np.arange(2):
        indices_ = np.where(np.log10(inference_result[:, index+11].astype(float)) < thresholds[0])[0]
        indices__ = np.where(np.log10(inference_result[:, index+11].astype(float)) > thresholds[1])[0]
        indices___ = list(set(indices_) | set(indices__))
        indices_infer = list(set(indices_infer) | set(indices___))
    for index in np.arange(2):
        indices_ = np.where(np.log10(inference_result[:, index+13].astype(float)) < thresholds[2])[0]
        indices__ = np.where(np.log10(inference_result[:, index+13].astype(float)) > thresholds[3])[0]
        indices___ = list(set(indices_) | set(indices__))
        indices_infer = list(set(indices_infer) | set(indices___))
    filtered_inference_result = np.delete(inference_result, indices_infer, axis = 0)
    return(filtered_inference_result)

def network_umap_pre(read_filename, save_filename):
    filtered_inference_result = filter_inference_results(read_filename, [-2.5, 2.0, -1.5, 2.0])
    result_ = np.vstack([filtered_inference_result[:, np.array([0, 1, 11, 13, 15, 21, 23])], filtered_inference_result[:, np.array([1, 0, 12, 14, 16, 20, 22])]])
    sorted_result_ = np.matrix(sorted(result_.tolist(), key = lambda x: (x[0], x[1])))
    # sorted_result_ = [genename_target, genename_tf, bf, bs, cv2, cmi, feedback_infer]
    index_interaction = np.where(sorted_result_[:, 6].astype(float) == 0)[0]
    sorted_result = np.delete(sorted_result_, index_interaction, axis = 0)
    unique_genename = sorted(list(set(np.vstack([np.array(sorted_result[:, 0]), np.array(sorted_result[:, 1])]).flatten().tolist())))
    counts_data_ = pd.read_csv(os.path.abspath('readfiles/scRNA-seq/countsdata matrix.csv'))
    counts_data = np.delete(np.matrix(counts_data_.values.tolist()), np.s_[:1], axis=1)
    counts_index = np.where(np.isin(counts_data[:, 0], unique_genename))[0]
    counts = counts_data[counts_index, 1::].astype(float)
    row_mean = np.nanmean(counts, axis=1)
    inds = np.where(np.isnan(counts))
    counts[inds] = np.take(row_mean, inds[0])
    counts = np.hstack([np.asarray(unique_genename).reshape([len(unique_genename), 1]), counts])
    df = pd.DataFrame(counts)
    df.to_csv(save_filename) 
    return (counts)

def network_pre(readfile_name, savefile_name):
    filtered_inference_result = filter_inference_results(readfile_name, [-3.5, 2.5, -2.0, 3.0])
    result_ = np.vstack([np.array(filtered_inference_result[:, np.array([1, 0, 23, 21])]), np.array(filtered_inference_result[:, np.array([0, 1, 22, 20])])])
    sorted_result = np.matrix(sorted(result_.tolist(), key = lambda x: (x[0], x[1])))
    unique_genename = sorted(list(set(np.array(sorted_result[:, 0]).flatten().tolist())))
    indices = np.where(sorted_result[:, 2].astype(float) == 0)[0]
    sorted_result = np.delete(sorted_result, indices, axis = 0)
    feedback = sorted_result.copy()
    index = 1
    for genename in unique_genename:
        indices_ = np.where(sorted_result[:, 0] == genename)[0]
        feedback[indices_, 0] = index
        indices__ = np.where(np.isin(unique_genename, sorted_result[indices_, 1]))[0]
        feedback[indices_, 1] = indices__ + 1
        index = index + 1
    df = pd.DataFrame(feedback)
    df.to_csv(savefile_name) 
    return(feedback)

def burst_pre(readfile_name, savefile_name):
    filtered_inference_result = filter_inference_results(readfile_name, [-2.5, 2.0, -1.5, 2.0])
    result_ = np.vstack([filtered_inference_result[:, np.array([0, 11, 13, 15, 21, 23])], filtered_inference_result[:, np.array([1, 12, 14, 16, 20, 22])]])
    sorted_result = np.matrix(sorted(result_.tolist(), key = lambda x: (x[0])))
    # sorted_result = [genename, bf, bs, cv2, cmi, feedback_infer]
    unique_genename = sorted(list(set(np.array(sorted_result[:, 0]).flatten().tolist())))
    result_unique_ = np.empty([len(unique_genename), 5])
    # result_unique_ = [bf, bs, cv2, expression_level, feedback]
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
    df = pd.DataFrame(result_unique_)
    df.to_csv(savefile_name) 
    return(result_unique_)

def barplot(vals, ylim, savefile_name):
    fig, ax = plt.subplots(dpi=900)
    categories = ['A', 'B']
    plt.bar(categories, vals)
    plt.title('Sample Bar Plot')
    plt.xlabel('Categories')
    plt.ylabel('Values')
    plt.ylim(ylim[0], ylim[1])
    plt.savefig(savefile_name)     
    plt.show()
    return

def gene_roles_heatmap(vals):
    data = np.asarray(np.vstack([vals[:, 0].astype(float), vals[:, 1].astype(float)]))
    scaler = MinMaxScaler()
    normalized_data = scaler.fit_transform(data.reshape(-1, 1))
    scaled_transformed = np.power(normalized_data, 2) 
    degree_data = normalized_data.reshape([vals.shape[1], vals.shape[0]])
    plt.figure(dpi=900)
    # Create a heatmap
    sns.heatmap(degree_data, annot=False, annot_kws={'size': 100}, cmap=plt.cm.Greens, square=True, cbar=False, xticklabels=False, yticklabels=False,)
    plt.xlabel('Gene roles', fontsize=3) 
    plt.title('Identification of gene roles in gene regulatory intercations network', fontsize=5, fontweight='bold')  
    plt.show()
    return

def subnetwork_pre(readfile_name, genename_, savefile_name):
    filtered_inference_result = filter_inference_results(readfile_name, [-3.5, 2.5, -2.0, 3.0])
    result_all = np.vstack([np.array(filtered_inference_result[:, np.array([1, 0, 23, 21])]), np.array(filtered_inference_result[:, np.array([0, 1, 22, 20])])])
    unique_genename_all = sorted(list(set(np.array(result_all[:, 0]).flatten().tolist())))
    indices1 = np.where(filtered_inference_result[:, 0] == genename_)[0]
    indices2 = np.where(filtered_inference_result[:, 1] == genename_)[0]
    indices = list(set(indices1) | set(indices2))
    filtered_inference_result = filtered_inference_result[indices, :]
    result_ = np.vstack([np.array(filtered_inference_result[:, np.array([1, 0, 23, 21])]), np.array(filtered_inference_result[:, np.array([0, 1, 22, 20])])])
    sorted_result = np.matrix(sorted(result_.tolist(), key = lambda x: (x[0], x[1])))
    indices = np.where(sorted_result[:, 2].astype(float) == 0)[0]
    sorted_result = np.delete(sorted_result, indices, axis = 0)
    unique_genename = sorted(list(set(np.array(sorted_result[:, 0]).flatten().tolist())))
    genename_index = np.where(np.isin(unique_genename_all, unique_genename))[0] + 1
    # df = pd.DataFrame(genename_index)
    # df.to_csv(savefile_name) 
    return(genename_index)

def scatter_plot(x, y, x_sd, y_sd):
    x_err = x_sd
    y_err = y_sd
    colors = ['blue', 'orange', 'grey', 'blue', 'orange', 'grey']
    categories = ['negative for low-expression', 'positive for low-expression', 'none for low-expression', 'negative for high-expression', 'positive for high-expression', 'none for high-expression']
    plt.rcParams.update({
        'font.size': 3, 'axes.labelsize': 3, 'axes.titlesize': 4,       
        'xtick.labelsize': 4,'ytick.labelsize': 4,         #
        'legend.fontsize': 2,       
        'lines.linewidth': 0.4,      
        'lines.markersize': 2, 'errorbar.capsize': 0.9,       
        'grid.linewidth': 0.4,       
        'figure.dpi': 300, 'figure.figsize': (2, 1.5)})
    fig, ax = plt.subplots()
    for i in range(3):
        plt.errorbar(x[i], y[i], xerr=x_err[i], yerr=y_err[i], fmt='o',  color=colors[i], markersize=1.8, capsize=1, capthick=1, label=categories[i])
    for i in range(3, len(x)):
        plt.errorbar(x[i], y[i], xerr=x_err[i], yerr=y_err[i], fmt='s',  color=colors[i], markersize=1.8, capsize=1, capthick=1, label=categories[i])
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)  
    ax.tick_params(axis='both', which='major', width=0.4) 
    ax.tick_params(axis='both', which='minor', width=0.4) 
    plt.xlabel('log10(bf)')
    plt.ylabel('log10(cv2)')
    plt.grid(True) 
    legend_handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=1.5, label='negative for low-expression'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', markersize=1.5, label='positive for low-expression'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='grey', markersize=1.5, label='none for low-expression'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='blue', markersize=1.5, label='negative for high-expression'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='orange', markersize=1.5, label='positive for high-expression'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='grey', markersize=1.5, label='none for high-expression')]
    plt.legend(handles=legend_handles, loc='best', handletextpad=0.5, fontsize=2)
    plt.show()
    return

def box_plot3(data_bf1, data_bf2, data_bs1, data_bs2, data_cv21, data_cv22):
    data_bf = [data_bf1, data_bf2]
    data_bs = [data_bs1, data_bs2]
    data_cv2 = [data_cv21, data_cv22]
    data = [data_bf, data_bs, data_cv2]
    fig, ax = plt.subplots(dpi=900)
    positions = [[1, 2], [4, 5], [7, 8]] 
    colors = ['skyblue', 'peachpuff']  
    for i, group_data in enumerate(data):
        for j, group in enumerate(group_data):
            ax.boxplot(group, positions=[positions[i][j]], patch_artist=True,
                       boxprops=dict(facecolor=colors[j], linewidth=3),
                       capprops=dict(linewidth=1.0), whiskerprops=dict(linewidth=1.0),  
                       medianprops=dict(linewidth=0.9, color='red'),  
                       showfliers=False)
    ax.set_xticks([1.5, 4.5, 7.5])
    ax.set_xticklabels(['BF', 'BS', 'CV2'], fontsize=3)
    ax.set_ylabel('log10(Value)', fontsize=3)
    plt.show()
    t_stat_bf, p_value_bf = stats.ttest_ind(data_bf1, data_bf2)
    t_stat_bs, p_value_bs = stats.ttest_ind(data_bs1, data_bs2)
    t_stat_cv2, p_value_cv2 = stats.ttest_ind(data_cv21, data_cv22)
    return(p_value_bf, p_value_bs, p_value_cv2)

def box_plot3_2(data1_bf1, data1_bf2, data1_bs1, data1_bs2, data1_cv21, data1_cv22, data2_bf1, data2_bf2, data2_bs1, data2_bs2, data2_cv21, data2_cv22):
    data1_bf = [data1_bf1, data1_bf2]
    data1_bs = [data1_bs1, data1_bs2]
    data1_cv2 = [data1_cv21, data1_cv22]
    data1 = [data1_bf, data1_bs, data1_cv2]
    data2_bf = [data2_bf1, data2_bf2]
    data2_bs = [data2_bs1, data2_bs2]
    data2_cv2 = [data2_cv21, data2_cv22]
    data2 = [data2_bf, data2_bs, data2_cv2]
    fig, axes = plt.subplots(1, 2, figsize=(4, 1.5), dpi=300) 
    positions = [[1, 2], [4, 5], [7, 8]] 
    colors = ['skyblue', 'peachpuff']  
    for i, group_data in enumerate(data1):
        for j, group in enumerate(group_data):
            axes[0].boxplot(group, positions=[positions[i][j]], patch_artist=True,
                       widths=0.6,
                       boxprops=dict(facecolor=colors[j], linewidth=0.5),
                       capprops=dict(linewidth=0.5), whiskerprops=dict(linewidth=0.3),  
                       medianprops=dict(linewidth=0.5, color='red'),  
                       showfliers=False)
    axes[0].set_xticks([1.5, 4.5, 7.5])
    axes[0].set_xticklabels(['BF', 'BS', 'CV2'], fontsize=6)
    axes[0].set_ylabel('log10(Value)', fontsize=5.5)
    axes[0].tick_params(axis='y', labelsize=5)
    for spine in axes[0].spines.values():
        spine.set_linewidth(0.5) 
    axes[0].tick_params(axis='both', width=0.5) 
    for i, group_data in enumerate(data2):
        for j, group in enumerate(group_data):
            axes[1].boxplot(group, positions=[positions[i][j]], patch_artist=True,
                       widths=0.6,
                       boxprops=dict(facecolor=colors[j], linewidth=0.5),
                       capprops=dict(linewidth=0.5), whiskerprops=dict(linewidth=0.3),  
                       medianprops=dict(linewidth=0.5, color='red'),  
                       showfliers=False)
    axes[1].set_xticks([1.5, 4.5, 7.5])
    axes[1].set_xticklabels(['BF', 'BS', 'CV2'], fontsize=6)
    axes[1].set_ylabel('log10(Value)', fontsize=5.5)
    axes[1].tick_params(axis='y', labelsize=5)
    for spine in axes[1].spines.values():
        spine.set_linewidth(0.5) 
    axes[1].tick_params(axis='both', width=0.5) 
    plt.tight_layout()
    plt.show()
    t_stat1_bf, p_value1_bf = stats.ttest_ind(data1_bf1, data1_bf2)
    t_stat1_bs, p_value1_bs = stats.ttest_ind(data1_bs1, data1_bs2)
    t_stat1_cv2, p_value1_cv2 = stats.ttest_ind(data1_cv21, data1_cv22)
    t_stat2_bf, p_value2_bf = stats.ttest_ind(data2_bf1, data2_bf2)
    t_stat2_bs, p_value2_bs = stats.ttest_ind(data2_bs1, data2_bs2)
    t_stat2_cv2, p_value2_cv2 = stats.ttest_ind(data2_cv21, data2_cv22)
    return(p_value1_bf, p_value1_bs, p_value1_cv2, p_value2_bf, p_value2_bs, p_value2_cv2)

def fit_3d(x, y, z):
    X = np.column_stack((x, y))
    X = sm.add_constant(X)  
    model = sm.OLS(z, X).fit()  
    intercept, a, b = model.params
    p_values = model.pvalues
    return(intercept, a, b, p_values[0], p_values[1], p_values[2])

def scatter_2d(x, y, color, x_label, y_label):
    fig, ax = plt.subplots(dpi=900)
    sc = plt.scatter(x, y, s=35, c=color, cmap=plt.cm.viridis)
    plt.colorbar(sc)
    plt.xlabel(x_label, fontsize=17)
    plt.ylabel(y_label, fontsize=17)         
    # plt.show()
    plt.clf()
    t_stat, p_value = stats.ttest_ind(x, y)
    return(p_value)

def scatter_interval(data):
    x1 = np.log10(data[1::, 0].astype(float))
    y1 = np.log10(data[1::, 2].astype(float))
    x2 = np.log10(data[1::, 1].astype(float))
    y2 = np.log10(data[1::, 3].astype(float))
    corr_coef1, p_value1 = pearsonr(x1, y1)
    corr_coef2, p_value2 = pearsonr(x2, y2)

    fig, axes = plt.subplots(1, 2, figsize=(4, 1.5), dpi=300) 
    sns.regplot(x=x1, y=y1, ci=95, scatter_kws={'alpha': 0.8, 's': 12}, line_kws={'linewidth': 0.8}, ax=axes[0])
    axes[0].set_xlabel('log10(bf_Fibr.)', fontsize=3.5)
    axes[0].set_ylabel('log10(bf_ES)', fontsize=3.5)
    axes[0].set_xlim(-1.2, 1.7)
    axes[0].set_ylim(-1.2, 1.7)
    axes[0].tick_params(axis='both', which='both', labelsize=3, width=0.3)
    sns.regplot(x=x2, y=y2, ci=95, scatter_kws={'alpha': 0.8, 's': 12}, line_kws={'linewidth': 0.8}, ax=axes[1])
    axes[1].set_xlabel('log10(bs_Fibr.)', fontsize=3.5)
    axes[1].set_ylabel('log10(bs_ES)', fontsize=3.5)
    axes[1].set_xlim(-0.7, 1.5)
    axes[1].set_ylim(-0.7, 1.5)
    axes[1].tick_params(axis='both', which='both', labelsize=3, width=0.3)
    for ax in axes:
        for spine in ax.spines.values(): spine.set_linewidth(0.4)
    plt.tight_layout()   
    plt.show()
    return(p_value1, p_value2)