import numpy as np
from scipy.special import erf, erfinv

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
    Parameters:
    trials (array-like): Trials to be considered.
    threshold (float): Threshold to be compared
    """
    pvalue = len(trials[trials >= threshold])/len(trials)
    return pvalue

def calculate_pvalues_from_fit(fit_parameter, fraction):

    return pvalue