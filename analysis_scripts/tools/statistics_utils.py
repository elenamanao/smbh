import numpy as np
from scipy.special import erf, erfinv
import scipy.stats as stats

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