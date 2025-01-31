import numpy as np
import pandas as pd
import sys
import os
import argparse

#load library for analysis
sys.path.append('tools')

from tools import randomize

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


p = argparse.ArgumentParser(description="Input parameters producing trials for correaltion analysis.")

p.add_argument("--sourcelist", 
               type=str, 
               help="Path to source list (hdf file)")
p.add_argument("--outdir", 
               type=str, 
               help="Path to the directory where to save trials")
p.add_argument("--mask_declination", 
               type = bool, 
               default = True, 
               help = 'Remove sources not visible by Pierre Auger Observatory')
p.add_argument("--mask_distance", 
               type = bool, 
               default = False, 
               help = 'Remove sources further than 50 Mpc')
p.add_argument("--pao_hotspot_treatment", 
               type = str,  
               help = '''How to mask the PAO hotspot. Here are the options:
               - no_mask: don't mask the hotspot
               - mask_and_scramble: first mask the hotspot and then scramble, making sure that no sources fall in the hotspot
               - scramble_and_mask_with_candidate_sources: first scramble and then remove sources that fall in the PAO hotspot
               - scramble_and_mask_without_candidate_sources: remove the souurces that are believed to be CR emitters, and then do same as scramble_and_mask_with_candidate_sources  '''
                )
p.add_argument("--n_trials",
               type = int, 
               default = 1000, 
               help = 'Number of trials we want to perform.')

args = p.parse_args()

# load source list

sourcelist_path = args.sourcelist
sources = pd.read_hdf(sourcelist_path, key = 'values')
n_sources_initial = len(sources)
print(f"The provided list contains {n_sources_initial} objects.")
n_trials = args.n_trials

mode = args.pao_hotspot_treatment

if np.isin([mode], ['no_mask', 'mask_and_scramble', 'scramble_and_mask_with_candidate_sources', 'scramble_and_mask_without_candidate_sources']):
    print("These trials will be produced for this case: ", mode)
else:
    print("The pao_hotspot_treatment you specified is not defined! Try checking for typos and reading the docstrings :))")

#apply masks
if args.mask_declination:
    print('Removing the sources that are not visible by PAO...')
    mask = sources["DEC_deg"] < 44.8 # maximum declination visible by PAO when including also inclined events
    sources = sources[mask]

if args.mask_distance:
    D_cut = 50
    print(f"Removing the sources more than {D_cut} Mpc away...")
    mask = sources["D"]< D_cut # Mpc
    sources = sources[mask]

#now the scrambling 
    
#define data type for the arrays where we will store our trials
# depending on the case, we will need to store bot r.a. and dec
# or r.a. only, so we define both and use only one of the two later
scrambled_data_dtype_ra_only = [('seed', '<i4'),
                      ('ra_scrambled',list )]

scrambled_data_dtype_ra_and_dec = [('seed', '<i4'),
                      ('ra_scrambled',list ),
                      ('dec_scrambled',list )]

pao_hotspot_ra, pao_hotspot_dec, pao_hotspot_r = 201.24634811, -45.37596794, 27
ra_true = sources.RA_deg.values #r.a. of the sources
dec_true = sources.DEC_deg.values #dec of the sources 

#first case, just scramble in r.a.
if mode == 'no_mask':
    scrambled_data_trials = np.zeros(n_trials, dtype = scrambled_data_dtype_ra_only)
    for s in np.arange(n_trials):
        ra_scramble = randomize.random_ra(ra_true, s)

        scrambled_data_trials['seed'][s] = s
        scrambled_data_trials['ra_scrambled'][s] = ra_scramble
    print(f"After applying the cuts, the list contains {len(ra_true)} objects.")

elif mode == 'mask_and_scramble':
    scrambled_data_trials = np.zeros(n_trials, dtype = scrambled_data_dtype_ra_only)
    #ra true and dec true will be overwritten by the function 
    for s in np.arange(n_trials):
        ra_scramble, ra_true, dec_true = randomize.mask_hotspot_and_randomize(ra_true, dec_true, s, 
                                                                              pao_hotspot_ra, pao_hotspot_dec, pao_hotspot_r )

        scrambled_data_trials['seed'][s] = s
        scrambled_data_trials['ra_scrambled'][s] = ra_scramble

    print(f"After applying the cuts, the list contains {len(ra_true)} objects.")

elif mode == 'scramble_and_mask_with_candidate_sources':
    scrambled_data_trials = np.zeros(n_trials, dtype = scrambled_data_dtype_ra_and_dec)
    for s in np.arange(n_trials):
        ra_scramble, dec_scramble = randomize.randomize_and_mask_hotspot(ra_true, dec_true, s,
                                                                         pao_hotspot_ra, pao_hotspot_dec, pao_hotspot_r )

        scrambled_data_trials['seed'][s] = s
        scrambled_data_trials['ra_scrambled'][s] = ra_scramble
        scrambled_data_trials['dec_scrambled'][s] = dec_scramble

elif mode == 'scramble_and_mask_without_candidate_sources':
    scrambled_data_trials = np.zeros(n_trials, dtype = scrambled_data_dtype_ra_and_dec)
    #mask the sources that are believed to be CRs emitters
    mask = np.logical_or(sources.Input.values == 'NGC4945', sources.Input.values == 'ESO 97-G13')
    ra_true = ra_true[~mask]
    dec_true = dec_true[~mask]
    print(f"After applying the cuts, the list contains {len(ra_true)} objects.")
    #ra true and dec true will be overwritten by the function 
    for s in np.arange(n_trials):
        ra_scramble, dec_scramble = randomize.randomize_and_mask_hotspot(ra_true, dec_true, s,
                                                                pao_hotspot_ra, pao_hotspot_dec, pao_hotspot_r )

        scrambled_data_trials['seed'][s] = s
        scrambled_data_trials['ra_scrambled'][s] = ra_scramble
        scrambled_data_trials['dec_scrambled'][s] = dec_scramble


#save scrambled trials
outfilename = 'trials_random_'+mode+f'_n_trials_{n_trials}.npy'

outdir = args.outdir

#check that directory exist
if not os.path.exists(outdir):
      os.makedirs(outdir)
      print("Directory created successfully!")
else:
      print("Directory already exists!")

out_path = os.path.join(outdir, outfilename)

print('Saving file as: ', out_path)
np.save(out_path ,scrambled_data_trials )




