'''Library to treat PAO public data'''

import numpy as np
from scipy.ndimage import gaussian_filter

def GreatCircleDistance(ra_1, dec_1, ra_2, dec_2, unit="rad"):
    """Compute the great circle distance between two events"""
    if unit == "deg":
        ra_1 = np.radians(ra_1)
        ra_2 = np.radians(ra_2)
        dec_1 = np.radians(dec_1)
        dec_2 = np.radians(dec_2)
    delta_dec = np.abs(dec_1 - dec_2)
    delta_ra = np.abs(ra_1 - ra_2)
    x = (np.sin(delta_dec / 2.0)) ** 2.0 + np.cos(dec_1) * np.cos(dec_2) * (
        np.sin(delta_ra / 2.0)
    ) ** 2.0
    return 2.0 * np.arcsin(np.sqrt(x))

def exposure(dec, theta_max=np.radians(60)):
    """
    Compute the Auger exposure as a function of declination, accounting for vertical and inclined events.

    Parameters:
    dec (array-like): Declination values (radians).
    theta_max (float): Maximum zenith angle for detection (radians).
    mode (str): "vertical" (default, up to 60°) or "inclined" (60°-80°).

    Returns:
    array-like: Relative exposure values.
    """
    l = np.radians(-35.23)  # Latitude of Malargüe, Argentina
    arg = (np.cos(theta_max) - np.sin(l) * np.sin(dec)) / (np.cos(l) * np.cos(dec))
    hm = np.arccos(arg.clip(-1, 1))

    exposure = np.cos(l) * np.cos(dec) * np.sin(hm) + hm * np.sin(l) * np.sin(dec)
    
    return exposure

def LoadExposureMap(theta_max_vert, theta_max_incl,  dec):

    # Compute the exposure_map for each pixel regarding its declination.
    exposure_map_total = exposure(dec, np.radians(theta_max_incl)) # Exposure from 0° to 80°
    exposure_map_vert = exposure(dec, np.radians(theta_max_vert))  # Exposure from 0° to 60°

    
    exposure_map_incl = (exposure_map_total - exposure_map_vert) # Exposure from 60° to 80°
    

    exposure_map_incl = exposure_map_incl / np.sum(exposure_map_incl) # Normalized
    
    exposure_map_vert = exposure_map_vert/np.sum(exposure_map_vert) # Normalized
    
    exposure_map = exposure_map_incl + exposure_map_vert # Total normalized exposure
    
    return exposure_map_vert, exposure_map_incl, exposure_map

def dipole_flux(alpha, delta, d, alpha_d, delta_d):
    """
    2D Dipole flux distribution function in RA and Dec.
    
    Parameters:
    alpha (float or np.array): Right ascension (RA) in radians.
    delta (float or np.array): Declination (Dec) in radians.
    d (float): Dipole amplitude.
    alpha_d (float): RA of dipole direction in radians.
    delta_d (float): Dec of dipole direction in radians.
    
    Returns:
    float or np.array: Normalized dipole flux at (alpha, delta).
    """
    # Compute dipole components
    d_x = d * np.cos(alpha_d) * np.cos(delta_d)
    d_y = d * np.sin(alpha_d) * np.cos(delta_d)
    d_z = d * np.sin(delta_d)
    
    dipole_modulation =(1 + d_x * np.cos(alpha) * np.cos(delta) + d_y * np.sin(alpha) * np.cos(delta) + d_z * np.sin(delta))
    return dipole_modulation

def smooth_flux(flux_values, sigma=1):
    """
    Smooth the flux to reduce any sharp concentration at the poles.
    """
    return gaussian_filter(flux_values, sigma=sigma)