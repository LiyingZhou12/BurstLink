# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import matplotlib.pyplot as plt


def select_genepair_grn(read_filename, save_filename, verbose):
    """ Select gene-pair interactions from base-GRN. """
    parquet_file = pq.ParquetFile(read_filename)
    base_GRN = parquet_file.read().to_pandas()
    tf_genename = base_GRN.columns.tolist()[2::]
    target_genename = base_GRN['gene_short_name'].tolist()
    tf_genename_lower = [s.lower() for s in tf_genename]
    target_genename_lower = [s.lower() for s in target_genename]
    base_grn = np.delete(np.matrix(base_GRN), np.s_[:2], axis=1)
    grn = find_genepair(tf_genename_lower, target_genename_lower, base_grn)
    if (verbose == True):
        df = pd.DataFrame(grn)
        df.to_csv(save_filename) 
        return (grn)
    else: return (grn)

def find_genepair(arr1, arr2, grn_matrix):
    grn = np.array(['gene1', 'gene2', 'feedback1', 'feedback2']).reshape([1, 4])
    genename = list(set(arr1) & set(arr2))
    indices1 = np.where(np.isin(arr1, genename))[0]
    indices2 = np.where(np.isin(arr2, genename))[0]
    for index1 in indices1:
        str1 = arr1[index1]
        for index2 in indices2:
            str2 = arr2[index2]
            index11 = arr2.index(str1)    
            index22 = arr1.index(str2)  
            grn_ = np.array([str1, str2, grn_matrix[index2, index1], grn_matrix[index11, index22]]).reshape([1, 4])
            grn = np.vstack([grn, grn_])        
    return (grn)

def integration_grn(read_filename_atac, read_filename_tf, save_filename, verbose):
    """ Integration for grn from different sources. """
    grn_atac = pd.read_csv(read_filename_atac)
    grn_tf = pd.read_csv(read_filename_tf)
    grn_atac = np.array(grn_atac)[:, 1::]
    grn_tf = np.array(grn_tf)
    filtered_grn_atac = filtered_genepair_grn(grn_atac)
    filtered_grn_tf = filtered_genepair_grn(grn_tf)
    integrated_grn = refiltered_genepair_grn(filtered_grn_atac, filtered_grn_tf)
    if (verbose == True):
        df = pd.DataFrame(integrated_grn)
        df.to_csv(save_filename) 
        return (integrated_grn)
    else:
        return (integrated_grn)
   
def filtered_genepair_grn(grn):
    grn_ = grn.tolist()
    sorted_grn = sorted(grn_, key = lambda x: (x[0], x[1], x[2], x[3]))
    sorted_grn_ = np.matrix(sorted_grn)
    genename1 = sorted_grn_[:, 0]
    genename2 = sorted_grn_[:, 1]
    feedback1 = sorted_grn_[:, 2].astype(float)
    feedback2 = sorted_grn_[:, 3].astype(float)
    # filter for same genens & multi feedback 
    genename1_diff = genename1[0: -1] == genename1[1: :]
    indices1 = np.where(genename1_diff == True)[0]
    genename2_diff = genename2[0: -1] == genename2[1: :]
    indices2 = np.where(genename2_diff == True)[0]
    feedback1_diff = feedback1[0: -1] - feedback1[1: :]
    indices1_ = np.where(feedback1_diff == 0)[0]
    feedback2_diff = feedback2[0: -1] - feedback2[1: :]
    indices2_ = np.where(feedback2_diff == 0)[0]
    indices = list(set(indices1) & set(indices2) & set(indices1_) & set(indices2_))
    filtered_grn = np.delete(sorted_grn_, indices, axis=0)
    # filter for same genes & feedback = 1
    genename1 = filtered_grn[:, 0]
    genename2 = filtered_grn[:, 1]
    genename1_diff = genename1[0: -1] == genename1[1: :]
    indices1 = np.where(genename1_diff == True)[0]
    genename2_diff = genename2[0: -1] == genename2[1: :]
    indices2 = np.where(genename2_diff == True)[0]
    indices = list(set(indices1) & set(indices2))
    filtered_grn = np.delete(filtered_grn, indices, axis=0)
    # filter for A-A same genenames
    genename1 = filtered_grn[:, 0]
    genename2 = filtered_grn[:, 1]
    same_genename = genename1 == genename2
    indices = np.where(same_genename == True)[0]
    filtered_grn = np.delete(filtered_grn, indices, axis=0)
    # filter for [genename1, genename2] & [genename2, genename1]
    genename1 = np.array(filtered_grn[:, 0]).flatten().tolist()
    unique_genename1 = list(set(genename1))
    sorted_unique_genename1 = sorted(unique_genename1)
    for uni_genename1 in sorted_unique_genename1:
        genename1 = filtered_grn[:, 0]
        genename2 = filtered_grn[:, 1]
        indices1 = np.where(np.isin(genename1, uni_genename1))[0]
        indices2 = np.where(np.isin(genename2, uni_genename1))[0]
        indices_ = np.where(np.isin(genename1[indices2], genename2[indices1]))[0]
        indices = indices2[indices_]
        filtered_grn = np.delete(filtered_grn, indices, axis=0)  
    return (filtered_grn)
  
def refiltered_genepair_grn(filtered_grn_atac, filtered_grn_tf):
    genename1_atac = filtered_grn_atac[:, 0]
    genename2_atac = filtered_grn_atac[:, 1]
    genename1_tf = filtered_grn_tf[:, 0]
    genename2_tf = filtered_grn_tf[:, 1]
    intersection1 = list(set(np.array(genename1_atac).flatten().tolist()) & set(np.array(genename1_tf).flatten().tolist()))
    intersection2 = list(set(np.array(genename2_atac).flatten().tolist()) & set(np.array(genename2_tf).flatten().tolist()))
    indices1 = np.where(np.isin(genename1_atac, intersection1))[0]
    indices2 = np.where(np.isin(genename2_atac, intersection2))[0]
    indices = list(set(indices1) & set(indices2))
    refiltered_grn_atac = np.delete(filtered_grn_atac, indices, axis=0)
    integrated_grn = np.vstack([refiltered_grn_atac, filtered_grn_tf]) 
    sorted_integrated_grn = np.matrix(sorted(integrated_grn.tolist(), key = lambda x: (x[0], x[1], x[2], x[3])))
    return (sorted_integrated_grn)
   
def RNAseq_analysis(read_filename1, read_filename2, save_filename1, save_filename2, verbose):
    """ scRNA-seq raw data preprocess. """
    df1 = pd.read_csv(read_filename1)
    df2 = pd.read_csv(read_filename2)
    rnaseq_list1 = df1.values.tolist()
    rnaseq_list2 = df2.values.tolist()
    rnaseq_matrix1 = np.delete(np.matrix(rnaseq_list1), np.s_[:1], axis=1)
    rnaseq_matrix2 = np.delete(np.matrix(rnaseq_list2), np.s_[:1], axis=1)
    genename1 = [row[0] for row in rnaseq_list1]
    genename1 = [s.lower() for s in genename1]
    genename2 = [row[0] for row in rnaseq_list2]
    genename2 = [s.lower() for s in genename2]
    if (genename1 == genename2):
        genename = genename1
        mm = len(genename)
        counts = np.hstack((rnaseq_matrix1, rnaseq_matrix2)).astype(float)
    else:
        print('genenames of two matrices can not match')

    stat = np.array(['mean', 'var', 'cv2']).reshape([1, 3])
    for index in np.arange(mm):
        countdata = counts[index, :]
        countdata_without_nan = countdata[~np.isnan(countdata)]
        mean_ = np.mean(countdata_without_nan)
        var_ = np.var(countdata_without_nan)
        cv2_ = var_ / (mean_**2)
        stat_ = np.array([mean_, var_, cv2_]).reshape([1, 3])
        stat = np.vstack([stat, stat_])
    
    mean = stat[1::, 0].astype(float)
    fig1, ax1 = plt.subplots(dpi=900)
    plt.hist(mean, bins=20, density=True)
    # plt.savefig('rawdata_mean.pdf')  
    plt.show()
    fig2, ax2 = plt.subplots(dpi=900)
    plt.hist(np.log(mean), bins=20, density=True)
    # plt.savefig('rawdata_logmean.pdf')  
    plt.show()
    fig3, ax3 = plt.subplots(dpi=900)
    df_mean = pd.DataFrame(np.log(mean))
    plt.boxplot(df_mean[0])
    # plt.savefig('rawdata_boxplot.pdf')  
    plt.show()
    
    counts_matrix = np.hstack([np.matrix(genename).reshape([mm, 1]), counts])
    stat_matrix = np.hstack([np.matrix(genename).reshape([mm, 1]), np.matrix(stat[1::, :]).astype(float)])
    if (verbose == True):
        df1 = pd.DataFrame(counts_matrix)
        df1.to_csv(save_filename1) 
        df2 = pd.DataFrame(stat_matrix)
        df2.to_csv(save_filename2) 
        return(counts_matrix, stat_matrix)
    else: return(counts_matrix, stat_matrix)
    
def selection_GRNandRNAseq(grn_filename, rnaseq_filename, counts_filename, threshold_value, verbose):
    """ Integration with scRNA-seq data and scATAC-seq data. """
    grn = pd.read_csv(grn_filename)
    grn = np.matrix(grn)[:, 1::]
    # cluster for feedback_GRN: (1, 1)&(1, 0)&(0, 0)
    dv = grn[:, 2] - grn[:, 3]
    indices_uni_feedback = [i for i, x in enumerate(dv) if abs(x) == 1]
    grn_uni_feedback = grn[indices_uni_feedback, :]
    indices1 = [i for i, x in enumerate(dv) if x == 0]
    indices2 = [i for i, x in enumerate(grn[:, 2]) if x == 0]
    indices_no_feedback = list(set(indices1) & set(indices2))
    grn_no_feedback = grn[indices_no_feedback, :]
    indices_ = list(set(indices_no_feedback) | set(indices_uni_feedback))
    grn_bio_feedback = np.delete(grn, indices_, axis=0)

    # selection single-cell data for mean
    stat_matrix = np.matrix(pd.read_csv(rnaseq_filename))
    counts_matrix = np.matrix(pd.read_csv(counts_filename))
    indices_mean = [i for i, x in enumerate(stat_matrix[:, 2]) if x >= threshold_value]
    selected_genename = stat_matrix[indices_mean, 1]
    selected_counts_matrix = counts_matrix[indices_mean, :]
    
    # comparison between grn and single-cell
    selected_genepair_bio_feedback = comparison(grn_bio_feedback, selected_genename)
    selected_genepair_uni_feedback = comparison(grn_uni_feedback, selected_genename)
    selected_genepair_no_feedback = comparison(grn_no_feedback, selected_genename)
    selected_genepair = np.vstack([selected_genepair_bio_feedback, selected_genepair_uni_feedback, selected_genepair_no_feedback]) 

    if (verbose == True):
        df1 = pd.DataFrame(selected_genepair)
        df1.to_csv('selected_genepair.csv') 
        df2 = pd.DataFrame(selected_counts_matrix)
        df2.to_csv('selected_counts_matrix.csv') 
        return(selected_genepair, selected_counts_matrix)
    else: return(selected_genepair, selected_counts_matrix)
    
def comparison(selected_grn, rna_seq_genename):
    """ Screen for gene pair names in both scRNA-seq data and scATAC-seq data. """
    genename1 = selected_grn[:, 0]
    genename2 = selected_grn[:, 1]
    indices1 = np.where(np.isin(genename1, rna_seq_genename))[0]
    indices2 = np.where(np.isin(genename2, rna_seq_genename))[0]
    indices = list(set(indices1) & set(indices2))
    selected_genepair = selected_grn[indices, :]
    return (selected_genepair)