'''Library with functions for scrambling sources in the Sky.
All the functions take the coordinates of the sources from a list
and perform a permutation of the r.a. of the array with a given seed.'''

import numpy as np

def GreatCircleDistance(ra_1, dec_1, ra_2, dec_2, unit="rad"):
    """Compute the great circle distance between 
    two points in the Sky"""
    if unit == "deg":
        ra_1 = np.radians(ra_1)
        ra_2 = np.radians(ra_2)
        dec_1 = np.radians(dec_1)
        dec_2 = np.radians(dec_2)
    delta_dec = np.abs(dec_1 - dec_2)
    delta_ra = np.abs(ra_1 - ra_2)
    x = (np.sin(delta_dec / 2.0)) ** 2.0 + np.cos(dec_1) * np.cos(dec_2) * (
        np.sin(delta_ra / 2.0)) ** 2.0
    return 2.0 * np.arcsin(np.sqrt(x))

def get_bin_centers(bins):
    bc = 0.5*(bins[1:]+bins[:-1])
    return bc
            
def do_simulation_bigged_seed(pdf_vertical, pdf_inclined, n_events_vert, n_events_incl, ra_bins, dec_bins, seed, n_vertical_events_final, n_inclined_events_final):
    
    # Compute the cumulative distribution function (CDF) # Normalize to form a probability density
    cdf_2d_vert = np.cumsum(pdf_vertical, axis=1)  # Compute along RA first
    cdf_2d_vert = np.cumsum(cdf_2d_vert, axis=0)  # Then along Dec
    cdf_2d_vert /= cdf_2d_vert[-1, -1]  # Normalize to form a probability density

    # Compute the cumulative distribution function (CDF) # Normalize to form a probability density
    cdf_2d_incl = np.cumsum(pdf_inclined, axis=1)  # Compute along RA first
    cdf_2d_incl = np.cumsum(cdf_2d_incl, axis=0)  # Then along Dec
    cdf_2d_incl /= cdf_2d_incl[-1, -1]  # Normalize to form a probability density

    # Generate uniform random numbers for sampling
    np.random.seed(seed)
    random_ra_vals = np.random.rand(n_events_vert)
    np.random.seed(seed+1)
    random_dec_vals = np.random.rand(n_events_vert)
    # Use interpolation to map random values to RA and Dec
    # Invert the CDF to get RA and Dec
    scrambled_ra_vert= np.interp(random_ra_vals, get_bin_centers(cdf_2d_vert.T[:, -1]), get_bin_centers(ra_bins))  # Solve for RA
    scrambled_dec_vert = np.interp(random_dec_vals, get_bin_centers(cdf_2d_vert.T[-1, :]), get_bin_centers(dec_bins))

    # Generate uniform random numbers for sampling
    np.random.seed(seed)
    random_ra_vals_incl = np.random.rand(n_events_incl)
    np.random.seed(seed+1)
    random_dec_vals_incl = np.random.rand(n_events_incl)
    # Use interpolation to map random values to RA and Dec
    # Invert the CDF to get RA and Dec
    scrambled_ra_incl = np.interp(random_ra_vals_incl, get_bin_centers(cdf_2d_incl.T[:, -1]), get_bin_centers(ra_bins))  # Solve for RA
    scrambled_dec_incl = np.interp(random_dec_vals_incl, get_bin_centers(cdf_2d_incl.T[-1, :]), get_bin_centers(dec_bins))

    # print(f'We simulated {n_events_vert} vertical events and {n_events_incl} inclined events, which matches the expected flux at 8 EeV.')

    #now we randomly sample a poissonian of the number of events above 20 EeV that we observe in our sample, and this will be our final simulation
    np.random.seed(seed)
    n_random_vertical_events = np.random.poisson(n_vertical_events_final)
    scrambled_ra_vert = scrambled_ra_vert[np.random.choice(len(scrambled_ra_vert), size=n_random_vertical_events, replace=False)]
    scrambled_dec_vert = scrambled_dec_vert[np.random.choice(len(scrambled_dec_vert), size=n_random_vertical_events, replace=False)]

    np.random.seed(seed+1)
    n_random_inclined_events = np.random.poisson(n_inclined_events_final)
    scrambled_ra_incl = scrambled_ra_incl[np.random.choice(len(scrambled_ra_incl), size=n_random_inclined_events, replace=False)]
    scrambled_dec_incl = scrambled_dec_incl[np.random.choice(len(scrambled_dec_incl), size=n_random_inclined_events, replace=False)]

    return scrambled_ra_vert, scrambled_dec_vert, scrambled_ra_incl, scrambled_dec_incl

def do_simulation(pdf_vertical, pdf_inclined, n_events_vert, n_events_incl, ra_bins, dec_bins, seed, n_vertical_events_final, n_inclined_events_final):

    # Compute CDFs
    cdf_2d_vert = np.cumsum(np.cumsum(pdf_vertical, axis=1), axis=0)
    cdf_2d_vert /= cdf_2d_vert[-1, -1]

    cdf_2d_incl = np.cumsum(np.cumsum(pdf_inclined, axis=1), axis=0)
    cdf_2d_incl /= cdf_2d_incl[-1, -1]

    # Create dedicated random generators
    rng_ra_vert = np.random.default_rng(seed + 100)
    rng_dec_vert = np.random.default_rng(seed + 101)
    rng_ra_incl = np.random.default_rng(seed + 200)
    rng_dec_incl = np.random.default_rng(seed + 201)
    rng_poisson_vert = np.random.default_rng(seed + 300)
    rng_poisson_incl = np.random.default_rng(seed + 301)
    rng_select_vert = np.random.default_rng(seed + 400)
    rng_select_incl = np.random.default_rng(seed + 401)

    # Sample RA and Dec
    random_ra_vals = rng_ra_vert.random(n_events_vert)
    random_dec_vals = rng_dec_vert.random(n_events_vert)
    scrambled_ra_vert = np.interp(random_ra_vals, get_bin_centers(cdf_2d_vert.T[:, -1]), get_bin_centers(ra_bins))
    scrambled_dec_vert = np.interp(random_dec_vals, get_bin_centers(cdf_2d_vert.T[-1, :]), get_bin_centers(dec_bins))

    random_ra_vals_incl = rng_ra_incl.random(n_events_incl)
    random_dec_vals_incl = rng_dec_incl.random(n_events_incl)
    scrambled_ra_incl = np.interp(random_ra_vals_incl, get_bin_centers(cdf_2d_incl.T[:, -1]), get_bin_centers(ra_bins))
    scrambled_dec_incl = np.interp(random_dec_vals_incl, get_bin_centers(cdf_2d_incl.T[-1, :]), get_bin_centers(dec_bins))

    # Final selection based on Poisson
    n_random_vertical_events = rng_poisson_vert.poisson(n_vertical_events_final)
    indices_vert = rng_select_vert.choice(len(scrambled_ra_vert), size=n_random_vertical_events, replace=False)
    scrambled_ra_vert = scrambled_ra_vert[indices_vert]
    scrambled_dec_vert = scrambled_dec_vert[indices_vert]

    n_random_inclined_events = rng_poisson_incl.poisson(n_inclined_events_final)
    indices_incl = rng_select_incl.choice(len(scrambled_ra_incl), size=n_random_inclined_events, replace=False)
    scrambled_ra_incl = scrambled_ra_incl[indices_incl]
    scrambled_dec_incl = scrambled_dec_incl[indices_incl]

    return scrambled_ra_vert, scrambled_dec_vert, scrambled_ra_incl, scrambled_dec_incl
