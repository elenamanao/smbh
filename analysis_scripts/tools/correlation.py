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

def fraction_of_sources_singleevent(ra_src, dec_src, evt_ra, evt_dec, search_radius):
    count = 0
    for ra_src, dec_src in zip(ra_src, dec_src):
        source_in_event_radius = GreatCircleDistance(ra_src, dec_src, evt_ra, evt_dec) < np.radians(search_radius)
        if np.any(source_in_event_radius):
                count = 1
    
    return count

def fraction_of_sources_allevents(source_ra, source_dec, data_ra, data_dec, search_radius):
    if len(data_ra) != len(data_dec):
        print('Check your data! RA and DEC must have the same length')
    elif len(data_ra) == len(data_ra):
        n_events = len(data_ra)

    fraction_eventlist = np.zeros(n_events)
    for i, (ra, dec) in enumerate(zip(data_ra, data_dec)):
        fraction_eventlist[i] = fraction_of_sources_singleevent(source_ra, source_dec, ra, dec, np.radians(search_radius))

    fraction = np.sum(fraction_eventlist)/n_events

    return fraction

def run_correlation(sources_ra, sources_dec, data_ra, data_dec, r_min, r_max, r_step, seed):
    steps = np.arange(r_min, r_max, r_step)

    results_dtype = [('seed', '<i4'),
                      ('ra',list ),
                      ('dec',list )]

    for i in steps:
        results_dtype.append((f'fraction_{i}', float))
    
    results = np.zeros(1, dtype = results_dtype)
    results['seed'] = np.ones_like(results['seed'])*seed
    results['ra'] = [data_ra]
    results['dec'] = [data_dec]
    for r in steps:
        fraction = fraction_of_sources_allevents(sources_ra, sources_dec, data_ra, data_dec, r)
        results[f'fraction_{r}'] = fraction
    
    return results

