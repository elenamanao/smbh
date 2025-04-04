import numpy as np
import matplotlib.pyplot as plt
import math

from scipy.stats import norm
from scipy.optimize import curve_fit

def gaussian(x, mu, sigma, A):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

def auto_grid_panels(n_panels):
    """Finds the best grid layout for a given number of panels."""
    cols = math.ceil(math.sqrt(n_panels))
    rows = math.ceil(n_panels / cols)
    return rows, cols

def plot_ts_distributions(trials,ax,i,mode, bins, norm = False):
        frac = i + 1
        histo, hedges = np.histogram(trials, bins = bins)
        bin_centers = (hedges[:-1] + hedges[1:]) / 2
        # Plot the fitted Gaussian
        if mode == 'AGN':
            color_histo = 'C0'
        if mode == 'nonAGN':
            color_histo = 'C1'
        if mode == 'all':
            color_histo = 'C2'

        #normalise histo counts
        normalization = 1
        if norm:
            normalization = np.sum(histo)
        histo= histo/normalization
        ax.step(bin_centers, histo, where='mid', alpha=0.3)
        ax.fill_between(bin_centers, histo, step="mid", alpha=0.3, label = f'{mode} trials')
        ax.set_xlabel('f[%]')

        return histo, hedges

def fit_gaussian(trials, bins):
    histo_counts, histo_hedges = np.histogram(trials, bins) 
    bin_centers = (histo_hedges[:-1] + histo_hedges[1:]) / 2
    params, covariance = curve_fit(gaussian, bin_centers, histo_counts, p0=[0, 1, 1])

    return params, covariance
     

def plot_gaussian_fit(params, ax, norm = False):
        x = np.linspace(0,1, 200)
        normalization = 1
        if norm:
            normalization = np.sum(x,gaussian(x, *params))
        ax.plot(x, gaussian(x, *params)/normalization, lw = 1)
