# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import random, math, warnings
from scipy import stats
from scipy.stats import poisson, norm, gaussian_kde, spearmanr
from scipy.special import j_roots
from scipy.special import beta as beta_fun
from scipy.optimize import minimize
from sklearn.metrics.cluster import normalized_mutual_info_score
import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D

from _utils import _synthetic_data


def global_burst_link(grn_filename, counts_filename, save_filename, verbose1, verbose2, verbose3, test, verbose4):
    """
    Inferface for inferring global dynamic information including gene regulatory interaction and bursting kinetics from scRNA-seq data of two genes.

    Parameters
    ----------
    grn_filename
        The optional input for base-GRN based on scATAC-seq.
    counts_filename
        The input for scRNA-seq data file.
    save_filename
        The saved inferred file name.
    verbose1
        Whether to output images.
        If False, only the inferred results are outputed. If True, the inferred results with visualization images are outputed.
    verbose2
        Whether to save images.
        If False, the images are only shown in the runtime screen. If True, the images are saved as pdf format in a folder under the default path.
    verbose3
        Whether to exhibit all images.
        If 0, all images are exhibited. If 1, only visualization images are exhibited.
    test
        Whether to perform distributional test.
        If False, distribution test is not performed. If True, gibbs samples are generated and distribution test is performed.
    verbose4
        Whether to save the inferred results.
        If False, the inferred results are not saved. If True, the inferred results are saved as csv format in a folder under the default path.
    
    Returns
    -------
    np.array(['gene1', 'gene2', 'feedback1', 'feedback2', 'kon1', 'kon2', 'koff1', 'koff2', 'ksyn1', 'ksyn2', 'w', 'bs1', 'bs2', 'bf1', 'bf2', 'cv1', 'cv2', 'convergence', 'corr', 'NMI', 'cmi_given_x1', 'cmi_given_x2', 'sign_givenx1', 'sign_givenx2'])
   """
    
    grn = pd.read_csv(grn_filename)
    counts = pd.read_csv(counts_filename)
    grn = np.matrix(grn)[:, 1::]
    counts = np.matrix(counts)[:, 2::]
    results =  np.array(['gene1', 'gene2', 'feedback1', 'feedback2', 'kon1', 'kon2', 'koff1', 'koff2', 'ksyn1', 'ksyn2', 'w', 'bs1', 'bs2', 'bf1', 'bf2', 'cv1', 'cv2', 'convergence', 'corr', 'NMI', 'cmi_given_x1', 'cmi_given_x2', 'sign_givenx1', 'sign_givenx2']).reshape([1, 24])
    num = len(grn[:, 0])
    for index in np.arange(num):
        try:
            geneinfo, counts_data = genepair_info(grn, counts, index)
            result_ = genepair_inference(counts_data, geneinfo, index, verbose1, verbose2, verbose3, test)
            result__ = np.hstack([geneinfo, result_.reshape([1, 20])])
            results = np.vstack([results, result__])
        except: continue
    if (verbose4 == True):
        df = pd.DataFrame(results)
        df.to_csv(save_filename) 
        return (results)
    else: return(results)
    
def global_uni_burst_link(grn_filename, counts_filename, genename1, genename2, verbose1, verbose2, verbose3, test):
    """
    Inferface for inferring uni dynamic information including gene regulatory interaction and bursting kinetics from scRNA-seq data of two genes.
    """
    grn = pd.read_csv(grn_filename)
    counts = pd.read_csv(counts_filename)
    grn = np.matrix(grn)[:, 1::]
    counts = np.matrix(counts)[:, 2::]
    try:
        index1 = np.where(np.isin(grn[:, 0], genename1))[0]  
        index2 = np.where(np.isin(grn[:, 1], genename2))[0]   
        index = list(set(index1) & set(index2))[0]
        result = _global_uni_burst_link(grn_filename, counts_filename, index, verbose1, verbose2, verbose3, test)
        return(result)
    except IndexError:
        print('can not find genepair information')
    
def _global_uni_burst_link(grn_filename, counts_filename, index, verbose1, verbose2, verbose3, test):
    grn = pd.read_csv(grn_filename)
    counts = pd.read_csv(counts_filename)
    grn = np.matrix(grn)[:, 1::]
    counts = np.matrix(counts)[:, 2::]
    geneinfo, counts_data =  genepair_info(grn, counts, index)
    result_ = genepair_inference(counts_data, geneinfo, index, verbose1, verbose2, verbose3, test)
    result = np.hstack([geneinfo, result_.reshape([1, 20])])
    return(result)

def genepair_info(grn, counts, index):
    gene1 = grn[index, 0]
    gene2 = grn[index, 1]
    geneinfo = grn[index, :]
    index1 = list(counts[:, 0]).index(gene1)    
    index2 = list(counts[:, 0]).index(gene2)    
    counts_data = np.vstack([counts[index1, 1::].astype(float), counts[index2, 1::].astype(float)])  
    nan_indices1 = np.where(np.isnan(counts_data[0, :]))[1]
    nan_indices2 = np.where(np.isnan(counts_data[1, :]))[1]
    nan_indices = list(set(nan_indices1) | set(nan_indices2))
    counts_data = np.delete(counts_data, nan_indices, axis=1)
    counts_data = np.array(counts_data).astype(int)
    return (geneinfo, counts_data)

def genepair_inference(vals, geneinfo, figflag, verbose1, verbose2, verbose3, test):
    """
    Infer dynamic information including gene regulatory interaction and bursting kinetics from scRNA-seq data of two genes.

    Parameters
    ----------
    vals
        The data matrix of shape 2 × columns.
        Rows correspond to genes or spots and columns to cells.
    geneinfo
        The priori information from the base-GRN built from scATAC-seq data, containing names and regulatory types of two genes.
    figflag
        The index of the outputed images name.
        If choose not to output the images, the default value is 0.
    verbose1
        Whether to output images.
        If False, only the inferred results are outputed. If True, the inferred results with visualization images are outputed.
    verbose2
        Whether to save images.
        If False, the images are only shown in the runtime screen. If True, the images are saved as pdf format in a folder under the default path.
    verbose3
        Whether to exhibit all images.
        If 0, all images are exhibited. If 1, only visualization images are exhibited.
    test
        Whether to perform distributional test.
        If False, distribution test is not performed. If True, gibbs samples are generated and distribution test is performed.
    
    Returns
    -------
    np.array(['gene1', 'gene2', 'feedback1', 'feedback2', 'kon1', 'kon2', 'koff1', 'koff2', 'ksyn1', 'ksyn2', 'w', 'bf1', 'bf2', 'bs1', 'bs2', 'cv1', 'cv2', 'convergence', 'corr', 'NMI', 'cmi_given_x1', 'cmi_given_x2', 'sign_givenx1', 'sign_givenx2'])

   """
   
    truncated_data = pre_truncation(vals)
    bursting_results = genepair_burstinference(truncated_data)
    # bursting_results = [alpha1, alpha2, beta1, beta2, phi1, phi2, w, bf1, bf2, bs1, bs2, cv1, cv2, convergence_success]
    interactions_results = genepair_interactionsinference(truncated_data, geneinfo, 0.4, figflag, verbose1, verbose2, verbose3)
    # interactions_results = [corr, nmi, cmi_given_x1, cmi_given_x2, sign_givenx1, sign_givenx2]
    if (bursting_results[-1] == 1): Binary_PoBe_Interactions(bursting_results[0: 7], truncated_data, figflag, geneinfo[0, 0], geneinfo[0, 1], verbose1, verbose2, verbose3)
    if (test == True):
        if (bursting_results[-1] == 1): ks_distance, h0 = distribution_fit(truncated_data, interactions_results[0: 7])
        elif (bursting_results[-1] == 0): ks_distance, h0 = np.array([0, 0])
        return (np.hstack([bursting_results, interactions_results, np.array(ks_distance, h0)]))
    else: return (np.hstack([bursting_results, interactions_results]))

def pre_truncation(vals):
    boundary_value1 = round(np.mean(vals[0, :]) + 2 * np.std(vals[0, :]))
    boundary_value2 = round(np.mean(vals[1, :]) + 2 * np.std(vals[1, :]))
    indices1 = np.where(vals[0, :] > boundary_value1)[0]
    indices2 = np.where(vals[1, :] > boundary_value2)[0]
    indices = list(set(indices1) & set(indices2))
    truncated_vals = np.delete(vals, indices, axis=1)
    # Select the HDR of the data
    hdr_data = select_hdr(truncated_vals)
    return (hdr_data)

def select_hdr(data):
    kde = gaussian_kde(data)
    density = kde(data)
    # Identify data points with density above the threshold
    density_threshold = np.sort(density)[round(0.01*len(data[0, :]))]
    hdr_indices = np.where(density >= density_threshold)[0]
    hdr_data = data[:, hdr_indices]
    return (hdr_data)

def genepair_burstinference(data):
    """
    Infer bursting kinetics for each gene.
    """
    params0 = initialization(data)
    bnds = ((1e-4, 1e3),)*4 + ((0, 1e4),)*2 + ((-1e3, 1e3),)
    # Constraint conditions
    cons = ({'type': 'ineq', 'fun': lambda params: (params[0] + params[2]) * (params[1] + params[3]) / max(params[0] * params[3], params[1] * params[2]) - params[6]},
            {'type': 'ineq', 'fun': lambda params: (params[0] + params[2]) * (params[1] + params[3]) / max(params[0] * params[1], params[2] * params[3]) + params[6]})
    vals_ = np.copy(data) 
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category = RuntimeWarning)
        ll = minimize(log_likelihood, params0, args = (vals_), method = 'SLSQP', constraints = cons, bounds = bnds)
    estim = ll.x
    success = ll.success.astype(str)
    if (success == 'True'): convergence = 1
    else: convergence = 0
    alpha1, alpha2, beta1, beta2, phi1, phi2, w = estim[0: 7]
    bf1 = 1 / (1 / estim[0] + 1 / estim[2])
    bf2 = 1 / (1 / estim[1] + 1 / estim[3])
    bs1 = estim[4] / estim[2]
    bs2 = estim[5] / estim[3]
    cv1 = np.var(data[0, :]) / ((np.mean(data[0, :]))**2)
    cv2 = np.var(data[1, :]) / ((np.mean(data[1, :]))**2)
    infer_results = np.array([alpha1, alpha2, beta1, beta2, phi1, phi2, w, bf1, bf2, bs1, bs2, cv1, cv2, convergence])
    return (infer_results)

def uni_poisson(vals, m):
    if (max(m) < 1e6): return (poisson.pmf(vals, m))
    else: return (norm.pdf(vals, loc=m, scale=math.sqrt(m)))

def uni_pobe(vals, alpha_, beta_, phi_):
    vals = vals.reshape([len(vals), 1])
    x, w = j_roots(50, beta_ - 1, alpha_ - 1)
    gs_ = np.sum(w * uni_poisson(vals, phi_ * (1 + x) / 2), axis=1)
    prob = 1 / (beta_fun(alpha_, beta_)+1e-12) * 2 ** (-alpha_ - beta_ + 1) * gs_
    return (prob)

def uni_moment_inference(vals):
    m1 = float(np.mean(vals))
    m2 = float(sum(vals * (vals - 1)) / len(vals))
    m3 = float(sum(vals * (vals - 1) * (vals - 2)) / len(vals))
    r1 = m1
    r2 = m2 / m1
    r3 = m3 / m2
    alpha_est = (2*r1*(r3-r2)) / (r1*r2-2*r1*r3+r2*r3)
    beta_est = (2*(r3-r2)*(r1-r3)*(r2-r1)) / ((r1*r2-2*r1*r3+r2*r3) * (r1-2*r2+r3))
    phi_est = (2*r1*r3-r1*r2-r2*r3) / (r1-2*r2+r3)
    return (np.array([alpha_est, beta_est, phi_est]))

def uni_loglikelihood(params, vals):
    alpha_, beta_, phi_ = params[0: 3]
    prob = uni_pobe(vals, alpha_, beta_, phi_)
    uni_loglikelihood_ = -np.sum(np.log(prob + 1e-12))
    return (uni_loglikelihood_)
  
def uni_maximum_likelihood(vals):
    x0 = uni_moment_inference(vals)
    bnds = ((1e-3, 1e3), (1e-3, 1e3), (1, 1e4))
    vals_ = np.copy(vals)  
    try: ll = minimize(uni_loglikelihood, x0, args=(vals_), method='L-BFGS-B', bounds=bnds)
    except: return np.array([np.nan, np.nan, np.nan])
    return (ll.x)

def initialization(vals):
    mle1 = uni_maximum_likelihood(vals[0, :])
    mle2 = uni_maximum_likelihood(vals[1, :])
    alpha1, beta1, phi1 = mle1[0: 3]
    alpha2, beta2, phi2 = mle2[0: 3]
    cov = np.cov(vals[0], vals[1])[0, 1]
    w = float((((alpha1 + beta1)**2) * ((alpha2 + beta2)**2) * (alpha1 + beta1 + 1) * (alpha2 + beta2 + 1) * cov) / (phi1 * phi2 * alpha1 * alpha2 * beta1 * beta2))  
    return (np.array([alpha1, alpha2, beta1, beta2, phi1, phi2, w]))

def log_likelihood(params, data):
    alpha1, alpha2 = params[0: 2]
    beta1, beta2 = params[2: 4]
    phi1, phi2 = params[4: 6]
    w = params[6]
    alpha1_ = alpha1 + 1
    alpha2_ = alpha2 + 1
    mu1 = alpha1 / (alpha1 + beta1)
    mu2 = alpha2 / (alpha2 + beta2)
    prob1 = uni_pobe(data[0, :], alpha1, beta1, phi1)
    prob2 = uni_pobe(data[1, :], alpha2, beta2, phi2)
    prob1_ = uni_pobe(data[0, :], alpha1_, beta1, phi1)
    prob2_ = uni_pobe(data[1, :], alpha2_, beta2, phi2)
    a = w * mu1 * mu2
    prob = prob1 * prob2 + a * (prob1_ - prob1) * (prob2_ - prob2)
    indices = np.where(prob < 1e-12)[0]
    prob_ = prob[indices] = 1.0
    result = (-np.log(prob_)).sum()
    return (result)

def genepair_interactionsinference(data, geneinfo, filter_eps, figflag, verbose1, verbose2, verbose3):
    """
    Infer gene regulatory interactions and regulation visualization.

    Parameters
    ----------
    data
        The data matrix of shape 2 × columns.
        Rows correspond to genes or spots and columns to cells.
    geneinfo
        The priori information from the base-GRN built from scATAC-seq data, containing names and regulatory types of two genes.
    filter_eps
        DREMI filter noise threshold and default value is 0.4
    figflag
        The index of the outputed images name.
        If choose not to output the images, the default value is 0.
    verbose1
        Whether to output images.
        If False, only the inferred results are outputed. If True, the inferred results with visualization images are outputed.
    verbose2
        Whether to save images.
        If False, the images are only shown in the runtime screen. If True, the images are saved as pdf format in a folder under the default path.
    verbose3
        Whether to  exhibit all images.
        If 0, all images are exhibited. If 1, only visualization images are exhibited.
    
    Returns
    -------
    np.array([corr, nmi, cmi_given_x1, cmi_given_x2, sign_givenx1, sign_givenx2])
    
   """
   
    density = density_estimation(data)
    if (verbose2 == True): heatmap(density, figflag, 1, verbose2, geneinfo[0, 0], geneinfo[0, 1])
    # Calculate conditional probability
    px2_given_x1 = density / np.sum(density, axis=0)
    px1_given_x2 = density.T / np.sum(density.T, axis=0)
    rescale_density_given_x1, sign_given_x1 = rescale_density(px2_given_x1, float(geneinfo[0, 2]), figflag, 2, verbose1, verbose2, verbose3, True, geneinfo[0, 0], geneinfo[0, 1])
    rescale_density_given_x2, sign_given_x2 = rescale_density(px1_given_x2, float(geneinfo[0, 3]), figflag, 3, verbose1, verbose2, verbose3, True, geneinfo[0, 1], geneinfo[0, 0])
    # Regulaotry interactions
    sample_given_x1 = resample(rescale_density_given_x1, filter_eps)
    sample_given_x2 = resample(rescale_density_given_x2, filter_eps)
    cmi_given_x1 = conditional_mutualinfo(sample_given_x1, figflag, 4, verbose1, verbose2, verbose3)
    cmi_given_x2 = conditional_mutualinfo(sample_given_x2, figflag, 5, verbose1, verbose2, verbose2)
    corr = np.corrcoef(data[0, :], data[1, :])[0, 1]
    nmi = normalized_mutual_info_score(data[0, :], data[1, :])
    results = np.array([corr, nmi, cmi_given_x1, cmi_given_x2, sign_given_x1, sign_given_x2])
    return (results)

def density_estimation(data):
    kde = gaussian_kde(data, bw_method='scott')
    x_max = np.max(data[0, :])
    y_max = np.max(data[1, :])
    x_grid, y_grid = np.mgrid[0: x_max: 1, 0: y_max: 1]
    grid_points = np.vstack([x_grid.ravel(), y_grid.ravel()])
    density_ = kde(grid_points)
    density = density_ / np.sum(density_)
    density = density.reshape([x_max, y_max]).T
    return (density)

def heatmap(vals, figflag, figindex, verbose2, genename1, genename2):
    fig, ax = plt.subplots(dpi=900)
    ax.set_title('Joint Probability distrbution of data')
    im = ax.pcolormesh(vals, cmap=plt.cm.coolwarm, alpha=0.9)
    cb = fig.colorbar(im)
    plt.xlabel(genename1)
    plt.ylabel(genename2)
    if (verbose2 == True):
        plt.savefig(str(figflag) + '-' + str(figindex) + '.pdf')
        plt.close()
        plt.clf()
    return

def rescale_density(density, info, figflag, figindex, verbose1, verbose2, verbose3, regulation_type, genename1, genename2):
    density_max = np.max(density, axis=0)
    nn = density.shape[1]
    indices = np.where(density_max == 0)[0]
    if (len(indices) > 0):
        for index in indices:
            if (index == nn - 1): density[:, index] = density[:, index - 1]
            else: density[:, index] = (density[:, index - 1] + density[:, index + 1]) / 2
        density_max = np.max(density, axis=0)
        nn = density.shape[1]       
    rescaled_density = np.zeros(density.shape)
    for col in range(nn):
        rescaled_density[:, col] = density[:, col] / density_max[col]  
    if (verbose1 == True) & ((verbose3 == 0) | (verbose3 == 1)):
        heatmap(rescaled_density, figflag, figindex, verbose2, genename1, genename2)
    # Identification for regulation type
    if (regulation_type == True):
        col = rescaled_density.shape[1]
        samples = np.array(['x', 'y']).reshape([2, 1])
        for index in np.arange(col):
            indices = np.where(rescaled_density[:, index] >= 0.5)[0]
            values = np.repeat(index, len(indices))
            samples_ = np.vstack([values, indices])
            samples = np.hstack([samples, samples_])
        samples = samples[:, 1::].astype(float)
        coef, p_value = spearmanr(samples[0, :], samples[1, :])
        if (info == 0) & (p_value >= 0.99): sign_type = 0
        else: sign_type = -1 * (coef < 0) + 1 * (coef > 0)
        return (rescaled_density, sign_type)
    else: return (rescaled_density)

def rejected_samples(prob, num_samples):
    if (np.sum(prob) == 0):
        sample = []
    else:
        probabilities_ = prob / np.sum(prob)
        values = np.where(probabilities_ > 0)[0] 
        probabilities = probabilities_[values]
        samples = []
        while len(samples) < num_samples:
            sample = random.choices(values, probabilities)[0]
            samples.append(sample)
    return (samples)

def resample(density, filter_eps):
    filtered_density = np.where(density > filter_eps, density, 0)
    row, col = np.shape(filtered_density)
    num_samples = 20
    xx = []
    yy = []
    for index in range(col):
        y = rejected_samples(filtered_density[:, index], num_samples)
        if (len(y) == num_samples):
            x = np.repeat(index, num_samples).reshape((num_samples, 1))
            xx.append(x)
            yy.append(y)
    xx = np.array(xx).flatten().T
    yy = np.array(yy).flatten().T
    re_samples = np.vstack((xx, yy))
    return (re_samples)

def density_counts(vals):
    xmax = np.max(vals[0])
    ymax = np.max(vals[1])
    count = np.zeros([xmax, ymax])
    col = vals.shape[1]
    for i in range(xmax):
        for j in range(ymax):
            cnt = 0
            for k in range(col):
                if i <= vals[0, k] and vals[0, k] < i + 1 and j <= vals[1, k] and vals[1, k] < j + 1:
                    cnt = cnt + 1
            count[i, j] = cnt
    count = np.transpose(count)       
    pxy = count / np.sum(np.sum(count))
    return(pxy)

def conditional_mutualinfo(data, figflag, figindex, verbose1, verbose2, verbose3):
    density = density_counts(data)
    px = np.sum(density, axis=0)
    py = np.sum(density.T, axis=0)
    py_givenx = density / px
    rescaled_prob = rescale_density(py_givenx, 0, figflag, figindex, verbose1, verbose2, 3, False, 'X1', 'X2')
    # Calculate CMI
    cmi_ = np.zeros(density.shape)
    for col in np.arange(len(px)):
        for row in np.arange(len(py)):
            if (density[row, col] == 0 or px[col] == 0 or py[row] == 0):
                cmi_[row, col] = 0
            else:
                cmi_[row, col] = density[row, col] / px[col]  * (np.log(density[row, col]) - np.log(px[col]) - np.log(py[row]))
    cmi = np.sum(np.sum(cmi_))
    return(cmi)

def Binary_PoBe_Interactions(params, vals, figflag, genename1, genename2, verbose1, verbose2, verbose3): 
    alpha1, alpha2, beta1, beta2, phi1, phi2, w = params[0: 7]
    mu1 = alpha1 / (alpha1 + beta1)
    mu2 = alpha2 / (alpha2 + beta2)
    alpha1_ = alpha1 + 1
    alpha2_ = alpha2 + 1
    density = density_estimation(vals)
    col = vals[0, :].max()
    row = vals[1, :].max()
    xx, yy = np.mgrid[0:col:1, 0:row:1]
    vals = np.vstack([xx.ravel(), yy.ravel()])
    n = np.shape(vals)[1]
    prob1 = uni_pobe(vals[0, :], alpha1, beta1, phi1)
    prob2 = uni_pobe(vals[1, :], alpha2, beta2, phi2)
    prob1_ = uni_pobe(vals[0, :], alpha1_, beta1, phi1)
    prob2_ = uni_pobe(vals[1, :], alpha2_, beta2, phi2)
    pxy = []
    py_givenx = []
    px_giveny = []
    for index in range(n):
        pxy_ = prob1[index] * prob2[index] + w * mu1 * mu2 * (prob1_[index] - prob1[index]) * (prob2_[index] - prob2[index])
        py_givenx_ = prob2[index] + w * mu1 * mu2 * (prob2_[index] - prob2[index]) * (prob1_[index] / prob1[index] - 1)
        px_giveny_ = prob1[index] + w * mu1 * mu2 * (prob1_[index] - prob1[index]) * (prob2_[index] / prob2[index] - 1)
        pxy.append(pxy_)
        px_giveny.append(px_giveny_)
        py_givenx.append(py_givenx_)
    pxy = np.array(pxy).reshape([col, row]).T
    py_givenx = np.array(py_givenx).reshape([col, row]).T
    px_giveny = np.array(px_giveny).reshape([col, row])
    rescale_py_givenx = rescale_density(py_givenx, 0, figflag, 0, False, False, 3, False, 'genename1', 'genename2')
    rescale_px_giveny = rescale_density(px_giveny, 0, figflag, 0, False, False, 3, False, 'genename2', 'genename1')
    px = np.sum(pxy, axis=0)
    py = np.sum(pxy, axis=1)
    if (verbose1 == True) & ((verbose3 == 0) | (verbose3 == 2)):
        if (verbose2 == False):
            fig, axes = plt.subplots(1, 3, figsize=(9, 3), dpi=300) 
            axes[0].pcolormesh(density, cmap=plt.cm.coolwarm, alpha=0.9)   
            axes[0].set_title("Joint Probability distrbution", fontsize=7)
            axes[0].set_xlabel(genename1)
            axes[0].set_ylabel(genename2)
            axes[1].pcolormesh(rescale_py_givenx, cmap=plt.cm.coolwarm, alpha=0.9)   
            axes[1].set_title("Rescaled Conditional Probability Given $X_1$", fontsize=7)
            axes[1].set_xlabel(genename1)
            axes[1].set_ylabel(genename2)
            axes[2].pcolormesh(rescale_px_giveny, cmap=plt.cm.coolwarm, alpha=0.9)
            axes[2].set_title("Rescaled Conditional Probability Given $X_2$", fontsize=7)
            axes[2].set_xlabel(genename2)
            axes[2].set_ylabel(genename1)
            plt.tight_layout()
            plt.show()
        else:
            fig1, ax1 = plt.subplots(dpi=900)
            ax1.set_title('Theory Inference Rescaled Conditional Probability Given $X_1$')
            im1 = ax1.pcolormesh(rescale_py_givenx, cmap=plt.cm.coolwarm, alpha=0.9)
            fig1.colorbar(im1, ax=ax1)
            plt.xlabel(genename1)
            plt.ylabel(genename2)
            plt.savefig('Rescaled conditional probability1.pdf')  
            plt.close()
            plt.clf()
            fig2, ax2 = plt.subplots(dpi=900)
            ax2.set_title('Theory Inference Rescaled Conditional Probability Given $X_2$')
            im2 = ax2.pcolormesh(rescale_px_giveny, cmap=plt.cm.coolwarm, alpha=0.9)
            fig2.colorbar(im2, ax=ax2)
            plt.xlabel(genename2)
            plt.ylabel(genename1)
            plt.savefig('Rescaled conditional probability2.pdf')  
            plt.close()
            plt.clf() 
    return
    
def distribution_fit(counts_data, stat_params):
    m = 500
    n = 500 + m + 1
    stat_data = _synthetic_data.gibbs_sample(np.array(stat_params), m, n)
    #  KW test
    data1 = np.transpose(counts_data[:, -500:])
    data2 = np.transpose(stat_data)
    S, KS_critical_val, H, p_value = ks_2samp(data1, data2, alpha=0.05)
    return(S, H)

def ks_2samp(x_val, y_val, alpha, asymptotic = False):
    """
    Performs a multivariate two-sample extension of the Kolmogorov-Smirnov test.

    Parameters
    ----------
    x_val
        A numpy array of shape (num_samples_x, dim) representing the first sample.
    y_val
        A numpy array of shape (num_samples_y, dim) representing the second sample.
    alpha
        The significance level.
    asymptotic
        Whether to use the asymptotic approximation or not.
    
    Returns
    -------
    A boolean indicating whether the null hypothesis is rejected.
    
   """
   
    num_samples_x, dim = x_val.shape
    num_samples_y, num_feats_y = y_val.shape
    if dim != num_feats_y:
        raise ValueError("The two samples do not have the same number of features.")
    z = np.zeros((num_samples_x, dim, dim))
    for h in range(dim):
        ind = np.argsort(x_val[:, h])[::-1]
        temp = np.take(x_val, ind, axis=0)
        z[:, :, h] = temp
        for i in range(dim):
            for j in range(num_samples_x - 1, -1, -1):
                if j == num_samples_x - 1:
                    runmax = temp[num_samples_x - 1, i]
                else:
                    runmax = max(runmax, temp[j, i])
                z[j, i, h] = runmax
    diff = np.zeros((num_samples_x, dim))
    for h in range(dim):
        for i in range(num_samples_x):
            val = np.abs(mecdf(x_val, z[i, :, h]) - mecdf(y_val, z[i, :, h])) * (
                round(num_samples_x * mecdf(x_val, z[i, :, h])) == num_samples_x - i
            )
            diff[i, h] = val
            if h == 0:
                diff[i, h] = max(
                    diff[i, h],
                    np.abs(mecdf(x_val, x_val[i, :]) - mecdf(y_val, x_val[i, :])),
                )
    KS = np.max(diff)
    if asymptotic:
        KS_critical_val = np.sqrt(-np.log(alpha / (4 * dim)) * (0.5 / num_samples_x)) + np.sqrt(
            (-1) * np.log(alpha / (4 * dim)) * (0.5 / num_samples_y)
        )
    else:
        KS_critical_val = np.sqrt(
            -np.log(alpha / (2 * (num_samples_x + 1) * dim)) * (0.5 / num_samples_x)
        ) + np.sqrt((-1) * np.log(alpha / (2 * (num_samples_y + 1) * dim)) * (0.5 / num_samples_y))
        
    if (KS > KS_critical_val):
        H = 'True'
        p_value = 1
    else:
        H = 'False'
        phi = -stats.norm.isf(1-alpha)
        sigma = KS_critical_val / phi
        p_value = 1 - stats.norm.sf(x = ((KS- 2*KS_critical_val) / sigma), loc = 0, scale = 1)
    return (KS, KS_critical_val, H, p_value)

def mecdf(x_val, t):
    """
    Computes the multivariate empirical cdf of x_val at t.

    Parameters
    ----------
    x_val
        A numpy array of shape (num_samples_x, dim) representing the sample.
    t
        A numpy array of shape (num_samples_t, dim) representing the point at which to evaluate the cdf.
    
    Returns
    -------
    The multivariate empirical cdf of x_val at t.  
    """

    lower = (x_val <= t) * 1.0
    return (np.mean(np.prod(lower, axis=1)))

def Binary_PoBe(params, vals): 
    alpha1, alpha2, beta1, beta2, phi1, phi2, w = params[0: 7]
    mu1 = alpha1 / (alpha1 + beta1)
    mu2 = alpha2 / (alpha2 + beta2)
    alpha1_ = alpha1 + 1
    alpha2_ = alpha2 + 1
    col = vals[0, :].max()
    row = vals[1, :].max()
    xx, yy = np.mgrid[0:col:1, 0:row:1]
    vals = np.vstack([xx.ravel(), yy.ravel()])
    n = np.shape(vals)[1]
    prob1 = uni_pobe(vals[0, :], alpha1, beta1, phi1)
    prob2 = uni_pobe(vals[1, :], alpha2, beta2, phi2)
    prob1_ = uni_pobe(vals[0, :], alpha1_, beta1, phi1)
    prob2_ = uni_pobe(vals[1, :], alpha2_, beta2, phi2)
    pxy = []
    py_givenx = []
    px_giveny = []
    for index in range(n):
        pxy_ = prob1[index] * prob2[index] + w * mu1 * mu2 * (prob1_[index] - prob1[index]) * (prob2_[index] - prob2[index])
        py_givenx_ = prob2[index] + w * mu1 * mu2 * (prob2_[index] - prob2[index]) * (prob1_[index] / prob1[index] - 1)
        px_giveny_ = prob1[index] + w * mu1 * mu2 * (prob1_[index] - prob1[index]) * (prob2_[index] / prob2[index] - 1)
        pxy.append(pxy_)
        px_giveny.append(px_giveny_)
        py_givenx.append(py_givenx_)
    pxy = np.array(pxy).reshape([col, row]).T
    py_givenx = np.array(py_givenx).reshape([col, row]).T
    px_giveny = np.array(px_giveny).reshape([col, row])
    px = np.sum(pxy, axis=0)
    py = np.sum(pxy, axis=1)
    return (pxy, px, py)