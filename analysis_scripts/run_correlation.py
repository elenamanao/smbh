import numpy as np
import pandas as pd
import sys
import os
import argparse

#load library for analysis
sys.path.append('tools')

from tools import scramble, correlation

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
p.add_argument("--cr_dataset", 
               type=str, 
               help="Path to cosmic rays dataset (npy file)")
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
p.add_argument("--mode", 
               type = str,
               default = 'simulation'  ,
               help = ''' Do you want to run the analysis on a simulation or on data?''')
p.add_argument("--search_radius",
               type = list,
               nargs='+',
               default = [1,30, 30],
               help = 'Input three values, the starting radius, the final radius and the step')
p.add_argument('--min_energy',
               type = float,
               default = 20,
               help = 'Minimum energy of the cosmic rays (in EeV)')  
p.add_argument("--seed",
               type = int, 
               default = 1, 
               help = 'Seed for scrambling')

args = p.parse_args()

# load source list
sourcelist_path = args.sourcelist
sources = pd.read_hdf(sourcelist_path, key = 'values')

# load cosmic rays dataset
cr_dataset_path = args.cr_dataset
cr_dataset = np.load(cr_dataset_path, allow_pickle=True)

n_sources_initial = len(sources)
print(f"The provided list contains {n_sources_initial} objects.")

# apply energy cut
mask_energy = cr_dataset["energy"] > args.min_energy
cr_dataset = cr_dataset[mask_energy]

print(f'The number of events in the dataset with energy above {args.min_energy} is: ', len(cr_dataset))

seed = args.seed
pao_hotspot_treatment = args.pao_hotspot_treatment
mode = args.mode 

if np.isin([pao_hotspot_treatment], ['no_mask', 'mask_and_scramble', 'scramble_and_mask_with_candidate_sources', 'scramble_and_mask_without_candidate_sources']):
    print("These trials will be produced for this case: ", pao_hotspot_treatment)
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
    mask = sources["D"] < D_cut # Mpc
    sources = sources[mask]

# define the range for angular search
r_min, r_max, r_step = args.search_radius

#now the scrambling 

pao_hotspot_ra, pao_hotspot_dec, pao_hotspot_r = 201.24634811, -45.37596794, 27
ra_true = sources.RA_deg.values #r.a. of the sources
dec_true = sources.DEC_deg.values #dec of the sources 

#first case, just scramble in r.a.
if pao_hotspot_treatment == 'no_mask':
    ra_scramble = scramble.scramble_ra(ra_true, seed) 
    results = correlation.run_correlation(ra_scramble, dec_true, cr_dataset['ra'], cr_dataset['dec'], r_min, r_max, r_step, seed)

elif pao_hotspot_treatment == 'mask_and_scramble':
    #ra true and dec true will be overwritten by the function 
    ra_scramble, ra_true, dec_true = scramble.mask_hotspot_and_scramble(ra_true, dec_true, seed,
                                                                pao_hotspot_ra, pao_hotspot_dec, pao_hotspot_r )
    
    results = correlation.run_correlation(ra_scramble, dec_true, cr_dataset['ra'], cr_dataset['dec'], r_min, r_max, r_step, seed)


elif pao_hotspot_treatment == 'scramble_and_mask_with_candidate_sources':

    ra_scramble, dec_scramble = scramble.scramble_and_mask_hotspot(ra_true, dec_true, seed,
                                                            pao_hotspot_ra, pao_hotspot_dec, pao_hotspot_r )
    
    results = correlation.run_correlation(ra_scramble, dec_true, cr_dataset['ra'], cr_dataset['dec'], r_min, r_max, r_step, seed)


elif pao_hotspot_treatment == 'scramble_and_mask_without_candidate_sources':
    #mask the sources that are believed to be CRs emitters
    mask = np.logical_or(sources["Input"]== 'NGC4945', sources["Input"] == 'ESO 97-G13')
    ra_true = ra_true[~mask]
    dec_true = dec_true[~mask]
    print(f"After applying the cuts, the list contains {len(ra_true)} objects.")
    #ra true and dec true will be overwritten by the function 
    ra_scramble, dec_scramble = scramble.scramble_and_mask_hotspot(ra_true, dec_true, seed,
                                                                pao_hotspot_ra, pao_hotspot_dec, pao_hotspot_r )

    results = correlation.run_correlation(ra_scramble, dec_true, cr_dataset['ra'], cr_dataset['dec'], r_min, r_max, r_step, seed)


#save scrambled trials
outfilename = 'trials_'+pao_hotspot_treatment+f'_seed_{seed}.npy'

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