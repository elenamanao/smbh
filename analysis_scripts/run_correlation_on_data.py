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
p.add_argument("--pao_hotspot_treatment", 
               type = str,  
               help = '''How to mask the PAO hotspot. Here are the options:
               - no_mask: don't mask the hotspot'''
                )
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

pao_hotspot_treatment = args.pao_hotspot_treatment
mode = args.mode 

if np.isin([pao_hotspot_treatment], ['no_mask', 'mask']):
    print("These trials will be produced for this case: ", pao_hotspot_treatment)
else:
    print("The pao_hotspot_treatment you specified is not defined! Try checking for typos and reading the docstrings :))")

#apply masks
if args.mask_declination:
    print('Removing the sources that are not visible by PAO...')
    mask = sources["DEC_deg"] < 44.8 # maximum declination visible by PAO when including also inclined events
    sources = sources[mask]

# define the range for angular search
r_min, r_max, r_step = args.search_radius

#now the scrambling 

pao_hotspot_ra, pao_hotspot_dec, pao_hotspot_r = 201.24634811, -45.37596794, 27

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
    if pao_hotspot_treatment == 'no_mask':
        results[f'fraction_{step}'] = correlation.fraction_of_sources_allevents(sources_ra_rad, 
                                                    sources_dec_rad, 
                                                    all_events_analysis['ra'], 
                                                    all_events_analysis['dec'], 
                                                    np.radians(step))

        if pao_hotspot_treatment == 'mask':
            #we mask both the data and the sources
            gcd_sources = auger_tools.GreatCircleDistance(sources_ra_rad, 
                                        sources_dec_rad,
                                        np.ones_like(sources['RA_deg'])*pao_hotspot_ra,
                                        np.ones_like(sources['DEC_deg'])*pao_hotspot_dec)
            #remove those
            check_distance_sources = gcd_sources > np.deg2rad(pao_hotspot_r)
            sources = sources[check_distance_sources]
            print(f"After cutting the sources in the PAO hotspot we are left with {len(sources)} sources.")

            gcd_events = auger_tools.GreatCircleDistance(results['ra'][i], 
                                        results['dec'][i],
                                        np.ones_like(results['ra'])*pao_hotspot_ra,
                                        np.ones_like(results['dec'])*pao_hotspot_dec)
                 #remove those
            check_distance_events = gcd_events > np.deg2rad(pao_hotspot_r)
            results['ra'] = results['ra'][check_distance_events]
            results['dec']= results['dec'][check_distance_events]
            print(f"After cutting the events in the PAO hotspot we are left with {len(results['dec'])} sources.")

            results[f'fraction_{step}'] = correlation.fraction_of_sources_allevents(sources_ra_rad, 
                                        sources_dec_rad, 
                                        all_events_analysis['ra'], 
                                        all_events_analysis['dec'], 
                                        step)


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