import numpy as np
from scipy.special import erf, erfinv
import scipy.stats as stats
import sys
sys.path.append('./')
from correlation import count_events_per_source

def sigma2pval(sigma, oneSided=False):
    r"""Converts gGaussian sigmas into p-values.
    Parameters
    ----------
    sigma : float or array_like
        Gaussian sigmas to be converted in p-values
    oneSided: bool, optional, default=False
        If the sigma should be considered one-sided or two-sided.
    Returns
    -------
    pval : array_like
        p-values.
    """

    pval = 1 - erf(sigma / np.sqrt(2))
    if oneSided:
        pval /= 2.0
    return pval

def pval2sigma(pval, oneSided=False):
    r"""Converts p-values into Gaussian sigmas.
    Parameters
    ----------
    pval : float or array_like
        p-values to be converted in Gaussian sigmas
    oneSided: bool, optional, default=False
        If the p-value should be considered one-sided or two-sided.
    Returns
    -------
    sigma : array_like
        Gaussian sigmas.
    """

    if oneSided:
        pval *= 2.0
    return np.sqrt(2) * erfinv(1 - pval)

def calculate_pvalue_from_trials(trials, threshold):
    """
    Calculate the p-value as the fraction of trials
    above the given threshold.
    Parameters
    --------
    trials (array-like): Trials to be considered.

    threshold (float): Threshold to be compared
    Returns
    -------
    float: The p-value, i.e., the probability of accepting 
        the null hypothesis.
    """
    mask = trials >= threshold
    pvalue = len(trials[mask])/len(trials)

    return pvalue

def calculate_pvalues_from_fit(fit_parameter, threshold):
    """
    Calculate the p-value as the integral of the gaussian 
    parametrised by the given fitted parameters above the given threshold.
    Parameters
    -------
    fit_parameter (array-like): Fitted parameters of the gaussian
    threshold (float): Threshold to be compared
    Returns
    -------
    float: The p-value, i.e., the probability of accepting the null hypothesis.
    """

    mean, std = fit_parameter  # Unpack Gaussian parameters
    pvalue = 1 - stats.norm.cdf(threshold, loc=mean, scale=std)

    return pvalue

def build_llh_distr(trials_AGN, trials_nonAGN):
    llh_dist  = np.zeros(len(trials_AGN))
    for i in np.arange(31250):
        pval_AGN_tmp = np.zeros(30)
        pval_nonAGN_tmp = np.zeros(30)
        for s in np.arange(1,31,1):
            pval_AGN_tmp[s-1] = calculate_pvalue_from_trials(trials_AGN[f'fraction_{s}'],trials_AGN[f'fraction_{s}'][i] )
            pval_nonAGN_tmp[s-1] = calculate_pvalue_from_trials(trials_nonAGN[f'fraction_{s}'],trials_nonAGN[f'fraction_{s}'][i] )
    
        llh_dist[i] = np.amin(pval_AGN_tmp)/np.amin(pval_nonAGN_tmp)

    return llh_dist

def compute_deviation_from_bkg(sources_dict, trials, radius):
    dict_trials = {'src_name': [],
               'counts': [],
               'n_events': [],
               'quantile16':[],
               'quantile84':[],
               'quantile3':[],
               'quantile97':[],
               'quantile0.1':[],
               'quantile99.9':[]}
    for i in np.arange(len(sources_dict['src_name'])):
        print('Running for source', sources_dict['src_name'][i])
        counts_array = []
        dict_trials['src_name'] = np.append(dict_trials['src_name'], sources_dict['src_name'][i])
        for s, seed in enumerate(trials['seed']):
            counts, ra, dec = count_events_per_source(np.radians(sources_dict['src_ra'][i]), np.radians(sources_dict['src_dec'][i]), trials['RA'][s], trials['Dec'][s], np.radians(radius))
            counts_array.append(counts)
            
        dict_trials['counts'] = np.append(dict_trials['counts'], counts_array)
        dict_trials['n_events'] = np.append(dict_trials['n_events'], np.median(counts_array))
        dict_trials['quantile16'] = np.append(dict_trials['quantile16'], np.quantile(counts_array, 0.16))
        dict_trials['quantile84'] = np.append(dict_trials['quantile84'], np.quantile(counts_array, 0.84))
        dict_trials['quantile3'] = np.append(dict_trials['quantile3'], np.quantile(counts_array, 0.03))
        dict_trials['quantile97'] = np.append(dict_trials['quantile97'], np.quantile(counts_array, 0.97))
        dict_trials['quantile0.1'] = np.append(dict_trials['quantile0.1'], np.quantile(counts_array, 0.001))
        dict_trials['quantile99.9'] = np.append(dict_trials['quantile99.9'], np.quantile(counts_array, 0.999))
    
    return dict_trials
