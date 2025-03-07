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
p.add_argument("--cr_vertical", 
               type=str, 
               help="Path to cosmic rays vertical dataset (npy file)")
p.add_argument("--cr_inclined", 
               type=str, 
               help="Path to cosmic rays inclined dataset (npy file)")
p.add_argument("--outdir", 
               type=str, 
               help="Path to the directory where to save trials")
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
p.add_argument("--mode", 
               type = str,
               default = 'simulation'  ,
               help = ''' Do you want to run the analysis on a simulation or on data?''')
p.add_argument("--search_radius",
               type = list,
               nargs='+',
               default = [1,30, 1],
               help = 'Input three values, the starting radius, the final radius and the step')
p.add_argument('--min_energy_analysis',
               type = float,
               default = 20,
               help = 'Minimum energy for the analysis')  
p.add_argument('--min_energy_dipole',
               type = float,
               default = 8,
               help = 'Minimum energy for the simulation of the events with the dipole')  
p.add_argument("--seed_initial",
               type = int, 
               default = 1, 
               help = 'Seed for simulating')
p.add_argument("--seed_final",
               type = int, 
               default = 100, 
               help = 'Seed for simulating')

args = p.parse_args()

# load source list
sourcelist_path = args.sourcelist
sources = pd.read_hdf(sourcelist_path, key = 'values')

if np.isin([args.filter_sources], ['all_sources', 'AGN', 'non_AGN']):
    print('The filter sources entry is not valid! Check the helpstring :)')

if args.filter_sources == 'all_sources':
    print('We run the correlation using all the sources in the list!')
elif args.filter_sources == 'AGN':
    print('Selecting only AGN in the list')
    mask = sources['Activity'] == True
    sources = sources[mask]
    print('We are left with', len(sources), 'objects in the list')
elif args.filter_sources == 'non_AGN':
    print('Selecting only non AGN in the list')
    mask = sources['Activity'] == True
    sources = sources[~mask]
    print('We are left with', len(sources), 'objects in the list')

# load cosmic rays dataset
cr_vertical_path = args.cr_vertical
cr_inclined_path = args.cr_inclined
cr_vertical = np.load(cr_vertical_path, allow_pickle=True)
cr_inclined = np.load(cr_inclined_path, allow_pickle=True)

n_sources_initial = len(sources)
# print(f"The provided list contains {n_sources_initial} objects.")

# apply energy cut
min_energy = args.min_energy_dipole
mask_energy_vertical = cr_vertical["energy"] > min_energy
mask_energy_inclined = cr_inclined["energy"] > min_energy
cr_vertical = cr_vertical[mask_energy_vertical]
cr_inclined = cr_inclined[mask_energy_inclined]

min_energy_analysis = args.min_energy_analysis
mask_energy_vertical_analysis = cr_vertical["energy"] > min_energy_analysis
mask_energy_inclined_analysis = cr_inclined["energy"] > min_energy_analysis
n_vertical_events_final = len(cr_vertical[mask_energy_vertical_analysis])
n_inclined_events_final = len(cr_inclined[mask_energy_inclined_analysis])

# print(f'The number of vertical events in the dataset with energy above {args.min_energy_dipole}EeV is: ', len(cr_vertical))
# print(f'The number of inclined events in the dataset with energy above {args.min_energy_dipole}EeV is: ', len(cr_inclined))

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
d,  alpha_d, delta_d = 0.074, np.radians(97), np.radians(-38)

ra_true = sources.RA_deg.values #r.a. of the sources
dec_true = sources.DEC_deg.values #dec of the sources 

ra_bins, dec_bins = np.linspace(0,2*np.pi, 1000), np.linspace(-1, 1, 1000)
ra_grid, dec_grid = np.meshgrid(ra_bins, dec_bins)

exposure_vert, exposure_incl, exposure_tot = auger_tools.LoadExposureMap(60, 80, np.arcsin(dec_grid))
exposure_vert_smooth = auger_tools.smooth_flux(exposure_vert, sigma = 15)
exposure_incl_smooth = auger_tools.smooth_flux(exposure_incl, sigma = 10)

flux_values = auger_tools.dipole_flux(ra_grid, np.arcsin(dec_grid), d,  alpha_d, delta_d)

pdf_vertical = flux_values*exposure_vert
pdf_inclined = flux_values*exposure_incl

seeds = np.arange(args.seed_initial, args.seed_final, 2)
results_dtype = [('seed', '<i4'),
                      ('ra',list ),
                      ('dec',list )]

results = np.zeros(len(seeds), dtype = results_dtype)

for i, seed in enumerate(seeds):

    sim_ra_vert, sim_sindec_vert, sim_ra_incl, sim_sindec_incl = simulation.do_simulation(pdf_vertical, pdf_inclined, len(cr_vertical), len(cr_inclined), ra_bins, dec_bins, seed, n_vertical_events_final, n_inclined_events_final)
    sim_ra = np.concatenate([sim_ra_vert, sim_ra_incl])
    sim_sindec = np.concatenate([sim_sindec_vert, sim_sindec_incl])

    sim_dec = np.arcsin(sim_sindec)

    results['seed'][i] = seed
    results['ra'][i] = sim_ra
    results['dec'][i] = sim_dec

#save trials
outfilename = 'trials_'+pao_hotspot_treatment+f'_minEnergy_{min_energy_analysis}EeV_seed_in_{args.seed_initial}_seed_fin_{args.seed_final}.npy'

outdir = args.outdir

#check that directory exist
if not os.path.exists(outdir):
    os.makedirs(outdir)
    print("Directory created successfully!")
else:
    print("Directory already exists!")

out_path = os.path.join(outdir, outfilename)

print('Saving file as: ', out_path)
np.save(out_path ,results)