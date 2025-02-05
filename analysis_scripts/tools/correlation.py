'''This library contains the function to run the correlation analysis'''

import numpy as np

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

def fraction_of_events_singlesource(ra_src, dec_src, data_ra, data_dec, search_radius):
    count = 0
    for ra_evt, dec_evt in zip(data_ra, data_dec):
        evts_in_source_radius = GreatCircleDistance(ra_src, dec_src, ra_evt, dec_evt, unit = 'deg') < np.radians(search_radius)
        if np.any(evts_in_source_radius):
                count = 1
    
    return count

def fraction_of_events_sourcelist(sources_ra, sources_dec, data_ra, data_dec, search_radius):
    if len(sources_ra) != len(sources_dec):
        print('Check your sources! RA and DEC must have the same length')
    elif len(sources_ra) == len(sources_dec):
        n_sources = len(sources_ra)

    fraction_sourcelist = np.zeros(n_sources)
    for i, (ra, dec) in enumerate(zip(sources_ra, sources_dec)):
        fraction_sourcelist[i] = fraction_of_events_singlesource(ra, dec, data_ra, data_dec, search_radius)

    fraction = np.sum(fraction_sourcelist)/n_sources

    return fraction

def run_correlation(sources_ra, sources_dec, data_ra, data_dec, r_min, r_max, r_step, seed):
    steps = np.linspace(r_min, r_max, r_step)

    results_dtype = [('seed', '<i4'),
                      ('ra_scrambled',list ),
                      ('dec_scrambled',list )]

    for i in steps:
        results_dtype.append((f'fraction_{i}', float))
    
    results = np.zeros(1, dtype = results_dtype)
    results['seed'] = np.ones_like(results['seed'])*seed
    results['ra_scrambled'] = [sources_ra]
    results['dec_scrambled'] = [sources_dec]
    for r in steps:
        fraction = fraction_of_events_sourcelist(sources_ra, sources_dec, data_ra, data_dec, r)
        results[f'fraction_{r}'] = fraction
    
    return results

