import numpy as np
import pandas as pd
import sys
import os
import argparse

#load library for analysis
sys.path.append('tools')

from tools import simulation, correlation, auger_tools

p = argparse.ArgumentParser(description="Input parameters producing trials for correaltion analysis.")

p.add_argument("--sourcelist", 
               type=str, 
               help="Path to source list (hdf file)")
p.add_argument("--vertical_data_path", 
               type=str, 
               help="Path to vertical data")
p.add_argument("--inclined_data_path", 
               type=str, 
               help="Path to inclined data")
p.add_argument("--outdir", 
               type=str, 
               help="Path to the directory where to save the results")
p.add_argument("--filter_sources",
               type = str,
               default = 'all_sources',
               help = ''''Which sources do you want to analyse:
               - all_sources: analyse all the sources in the list
               - AGN: take only AGN in the list
               - no_AGN: take only non-AGN in the list''')
p.add_argument("--mask_declination", 
               type = bool, 
               default = True, 
               help = 'Remove sources not visible by Pierre Auger Observatory')
p.add_argument("--search_radius",
               type = list,
               nargs='+',
               default = [1,31, 1],
               help = 'Input three values, the starting radius, the final radius and the step')

args = p.parse_args()

# load source list
sourcelist_path = args.sourcelist
sources = pd.read_hdf(sourcelist_path, key = 'values')

if np.isin([args.filter_sources], ['all_sources', 'AGN', 'non_AGN']):
    print('The filter sources entry is not valid! Check the helpstring :)')

filter_sources = args.filter_sources

if filter_sources == 'all_sources':
    print('We run the correlation using all the sources in the list!')
elif filter_sources == 'AGN':
    print('Selecting only AGN in the list')
    mask = sources['Activity'] == True
    sources = sources[mask]
    print('We are left with', len(sources), 'objects in the list')
elif filter_sources == 'non_AGN':
    print('Selecting only non AGN in the list')
    mask = sources['Activity'] == True
    sources = sources[~mask]
    print('We are left with', len(sources), 'objects in the list')

n_sources_initial = len(sources)

inclined_events = np.load(args.inclined_data_path, allow_pickle = True)
vertical_events = np.load(args.vertical_data_path, allow_pickle = True)

all_events = np.concatenate((inclined_events, vertical_events))

mask_energy = all_events['energy'] > 20

all_events_analysis = all_events[mask_energy]

#apply masks
if args.mask_declination:
    print('Removing the sources that are not visible by PAO...')
    mask = sources["DEC_deg"] < 44.8 # maximum declination visible by PAO when including also inclined events
    sources = sources[mask]

# define the range for angular search
r_min, r_max, r_step = args.search_radius

ra_true = sources.RA_deg.values #r.a. of the sources
dec_true = sources.DEC_deg.values #dec of the sources 

results_dtype = [     ('ra',list ),
                      ('dec',list )]


steps = np.arange(r_min, r_max, r_step)
for i in steps:
    results_dtype.append((f'fraction_{i}', float))

results = np.zeros(1, dtype = results_dtype)
sources_ra_rad, sources_dec_rad = np.radians(sources['RA_deg']), np.radians(sources['DEC_deg'])

results['ra'] = all_events_analysis['ra']
results['dec']= all_events_analysis['dec']

for step in steps:
    #first case, just scramble in r.a.
    results[f'fraction_{step}'] = correlation.fraction_of_sources_allevents(sources_ra_rad, 
                                                sources_dec_rad, 
                                                all_events_analysis['ra'], 
                                                all_events_analysis['dec'], 
                                                np.radians(step))

    #save trials
    outfilename = 'correlation_'+filter_sources+f'_results.npy'

    outdir = args.outdir

    #check that directory exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        print("Directory created successfully!")
    else:
        print("Directory already exists!")

    out_path = os.path.join(outdir, outfilename)

    print('Saving file as: ', out_path)
    np.save(out_path ,results )