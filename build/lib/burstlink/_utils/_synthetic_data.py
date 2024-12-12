# -*- coding: utf-8 -*-
import random, copy, math
import numpy as np
from types import SimpleNamespace
from scipy.interpolate import interp1d
from scipy.stats import poisson,norm,uniform
from scipy.special import j_roots
from scipy.special import beta as beta_fun
import seaborn as sns
import matplotlib.pyplot as plt


def SSA_coexpression(params, verbose, fig):
    """
    Generate sample data based on gene-gene interactions with bursting mechanistic dynamic model.

    Parameters
    ----------
    params
        The parameters of mechanistic model.
        params = [kon1, kon2, koff1, koff2, ksyn1, ksyn2, kdeg, rK, h1, h2, eps_, T].
    verbose
        Whether to output busrst kinetics.
        If 'burst', return with busrst kinetics.
        The priori information from the base-GRN built from scATAC-seq data, containing names and regulatory types of two genes.
    fig
        Whether to exhibit distribution images.
        If False, no distribution images are shown. If True, the images of joint distribution and marginal distributions are shown.
    
    Returns
    -------
    Sample data based on gene-gene interactions with bursting mechanistic model and its corresponding parameters.

   """
    
    kon1, kon2, koff1, koff2, ksyn1, ksyn2, kdeg = params[0: 7]  
    r_K, h1, h2, eps_ = params[7: 11]  
    T = params[11]
                             #  OFFx  ONx  OFFy  ONy   X    Y
    reactionmatrix = np.matrix([[-1,   1,   0,    0,   0,   0],
                                [ 1,  -1,   0,    0,   0,   0],
                                [ 0,   0,  -1,    1,   0,   0],
                                [ 0,   0,   1,   -1,   0,   0],
                                [ 0,   0,   0,    0,   1,   0],
                                [ 0,   0,   0,    0,   0,   1],
                                [ 0,   0,   0,    0,  -1,   0],
                                [ 0,   0,   0,    0,   0,  -1]])
    m = np.shape(reactionmatrix)[0]
    S = np.matrix([1, 0, 1, 0, 2, 2])  # initial states
    t = [0]
    tt = 0
    interval_time = []   
    bs1 = []
    bs2 = []
    kkonx = []
    kkony = []
    while (tt < T):
        if (h1 < 0):
            if (S[-1, 5] == 0): konx = kon1 * eps_
            else: konx = kon1 * (1 / (1 + 1/((r_K * S[-1, 5]) ** (-h1))) + eps_)
        else:
            konx = (h1 == 0) * kon1 * (1/2 + eps_) + (S[-1, 5] == 0) * (kon1 * (1 + eps_)) * (h1 > 0) + (S[-1, 5] > 0) * kon1 * (1 / (1 + (r_K * S[-1, 5])**h1) + eps_) * (h1 > 0) 
        if (h2 < 0):
            if (S[-1, 4] == 0): kony = kon2 * eps_
            else:  kony = kon2 * (1 / (1 + 1/((r_K * S[-1, 4]) ** (-h2))) + eps_)
        else:
            kony = (h2 == 0) * kon2 * (1/2 + eps_) + (S[-1, 4] == 0) * (kon2 * (1 + eps_)) * (h2 > 0) + (S[-1, 4] > 0) * kon2 * (1 / (1 + (r_K * S[-1, 4])**h2) + eps_) * (h2 > 0)
        PF = np.array([konx * (S[-1, 0]==1), koff1 * (S[-1, 1] == 1), kony * (S[-1, 2] == 1), koff2 * (S[-1,3] == 1), ksyn1 * (S[-1,1] == 1), ksyn2 * (S[-1,3] == 1), kdeg * S[-1, 4], kdeg * S[-1, 5]])
        p0 = np.sum(PF)
        
        r1 = random.random()
        tau = (1/p0) * np.log(1/r1)
        r2 = random.random()
        for index in np.arange(m):
            if (np.sum(PF[0:index+1]) >= r2*p0):
                next_r = index   
                break
        tt = tt + tau
        t.append(tt)
        gain = reactionmatrix[next_r, :]
        last = S[-1, :] + gain
        S = np.concatenate((S, last), axis=0)
        interval_time.append(tau)
        bs2.append(next_r == 5)
        bs1.append(next_r == 4)
        bs2.append(next_r == 5)
        kkonx.append(konx)
        kkony.append(kony)
    t = np.asarray(t)
    tq = np.arange(T-1000, T, 0.1)
    xx = np.array(S[:, 4]).flatten()
    x_interpfunc = interp1d(t, xx, 'previous')
    xx_ = x_interpfunc(tq)  
    yy = np.array(S[:, 5]).flatten()
    y_interpfunc = interp1d(t, yy, 'previous')
    yy_ = y_interpfunc(tq)
    S_stable = np.vstack((xx_, yy_))
    
    if (fig == True):
        fig, ax = plt.subplots(dpi=900)
        g = sns.jointplot(x=S_stable[0, :], y=S_stable[1, :], kind="kde", fill=True, marginal_kws=dict(fill=True), cut=0)
        g.ax_joint.set_xlim([0, np.max(S_stable[0, :])])
        g.ax_joint.set_ylim([0, np.max(S_stable[1, :])])
        g.ax_marg_x.set_xlim([0, np.max(S_stable[0, :])])
        g.ax_marg_y.set_ylim([0, np.max(S_stable[1, :])])
        g.ax_joint.set_xlabel('$X_1$')
        g.ax_joint.set_ylabel('$X_2$')    
        plt.show()
    if (verbose == 'burst'):
        D1 = np.diff(S[:, 1], axis = 0)
        D2 = np.diff(S[:, 3], axis = 0)
        index1 = np.where(D1 == 1)[0] + 1
        index2 = np.where(D2 == 1)[0] + 1
        intervaltime_1 = []
        intervaltime_2 = []
        bs_1 = []
        bs_2 = []
        nn = np.min(len(index1))
        nn0 = int(np.round(nn / 3))
        for idx in np.arange(nn0, nn): 
            intervaltime_1.append(np.sum(interval_time[index1[idx - 1]: index1[idx]]))
            bs_1.append(np.sum(bs1[index1[idx - 1]: index1[idx]]))
        nn = np.min(len(index2))
        nn0 = int(np.round(nn / 3))
        for idx in np.arange(nn0, nn): 
            intervaltime_2.append(np.sum(interval_time[index2[idx - 1]: index2[idx]]))
            bs_2.append(np.sum(bs2[index2[idx - 1]: index2[idx]]))
        mbf1 = 1 / np.mean(intervaltime_1)
        mbf2 = 1 / np.mean(intervaltime_2)
        mbs1 = np.mean(bs_1)
        mbs2 = np.mean(bs_2)
        mkonx = np.mean(kkonx[1000::])
        mkony = np.mean(kkony[1000::])
        mbf1_ = 1 / (1 / mkonx + 1 / koff1)
        mbf2_ = 1 / (1 / mkony + 1 / koff2)
        mbs1_ = ksyn1 / koff1
        mbs2_ = ksyn2 / koff2
        return SimpleNamespace(S_stable=S_stable, mbf1=mbf1, mbf2=mbf2, mbs1=mbs1, mbs2=mbs2)
    else: return (S_stable)
        
def gibbs_sample(params, m, n):
    """
    Generate sample data based on the bivaraite Poisson-Beta distribution via gibbs sampling.

    Parameters
    ----------
    params
        The parameters of statistical model.
    m
        The number of sample.
    n
        The number of iterations.
    
    Returns
    -------
    Sample data based on gene-gene interactions with bursting statistical model.

   """
    
    X = []
    Y = []
    y = [15]
    for i in range(n):
            x = sample_x(params, y[0], 1)
            y = sample_y(params, x[0], 1)
            if (i >= m+1):
                X.append(x)
                Y.append(y)
    X = np.array(X)
    Y = np.array(Y)
    X.shape = (1, len(X))
    Y.shape = (1, len(Y))
    result = np.vstack([X, Y])
    sample = np.array(result)
    return(sample)

def PoissonPMF(A,lam):
    if(max(lam) < 1e6):
        return(poisson.pmf(A,lam))
    else: return(norm.pdf(A,loc=lam,scale=math.sqrt(lam)))

def PoBePDF(A, alpha, beta, phi):    
    row,col = np.shape(A)
    prob = np.zeros((row,col))
    for i in range(row):
        aa = A[i] 
        aa = np.array(aa)
        aa.shape = (len(aa),1)
        x,w = j_roots(col, beta[i]-1, alpha[i]-1)
        GJInt = np.sum(w*PoissonPMF(aa, phi[i]*(1+x)/2), axis=1) 
        prob[i] = 1/beta_fun(alpha[i], beta[i])*2**(-float(alpha[i])-float(beta[i])+1)*GJInt
    return(prob)

def ConditionProb(params, vals): 
        row,col = np.shape(vals)
        a = params[0:2]
        b = params[2:4]
        phi = params[4:6]
        w = params[6]
        mu1 = a[0]/(a[0]+b[0])
        mu2 = a[1]/(a[1]+b[1])
        ac = copy.deepcopy(a)
        for i in range(len(a)):
            ac[i] = a[i] + 1
        uniprob = PoBePDF(vals, a, b, phi)
        uniprob1 = PoBePDF(vals, ac, b, phi)
        p_xgiveny = []
        p_ygivenx = []
        for i in range(col):
            p1 = uniprob[0,i] + w * mu1 * mu2 * (uniprob1[0,i] - uniprob[0,i]) * (uniprob1[1,i] / uniprob[1,i] - 1)
            p2 = uniprob[1,i] + w * mu1 * mu2 * (uniprob1[1,i] - uniprob[1,i]) * (uniprob1[0,i] / uniprob[0,i] - 1)
            p_xgiveny.append(p1)
            p_ygivenx.append(p2)
        px = np.array(p_xgiveny)
        py = np.array(p_ygivenx)
        px.shape = (1,col)
        py.shape = (1,col)
        result = np.vstack([px,py])
        return(result)

def sample_x(params, y, num):
    # Inverse Sampling
    sample = [0 for ss in range(num)]
    Y = [y for ss in range(200)]
    Y = np.array(Y)
    vals = np.vstack((np.arange(200), Y))
    p = ConditionProb(params, vals)[0,]
    for i in range(num):
        u = uniform.rvs(size=1)
        for j in range(200):
            index = np.arange(j+1)
            cdf = np.sum(p[index])
            if (u <= cdf):
                sample[i] = j
                break
    sample = np.array(sample)
    return(sample)

def sample_y(params, x, num): 
    sample = [0 for ss in range(num)]
    X = [x for ss in range(200)]
    X = np.array(X)
    vals = np.vstack([X, np.arange(200)])
    p = ConditionProb(params,vals)[1,]
    for i in range(num):
        u = uniform.rvs(size=1)
        for j in range(200):
            index = np.arange(j+1)
            cdf = np.sum(p[index])
            if (u <= cdf):
                sample[i] = j
                break    
    sample = np.array(sample)
    return(sample)

