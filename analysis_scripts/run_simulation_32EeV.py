import numpy as np
import pandas as pd
import sys
import os
import argparse

#load library for analysis
sys.path.append('tools')

from tools import simulation, correlation, auger_tools

p = argparse.ArgumentParser(description="Input parameters producing trials for correaltion analysis.")

p.add_argument("--cr_file", 
               type=str, 
               help="Path to cosmic rays vertical dataset (dat file)")
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
               help = 'First seed to use')
p.add_argument("--seed_final",
               type = int, 
               default = 100, 
               help = 'Last seed to use')
p.add_argument('--mode',
               type = str,
               default = 'sample',
               help = '''default is sample: simulate only the events in the sample (~600)
                If you want to so a simulation forecasting what would happen with the full sample
                of the events, pass 'all_statistics' as argument 
               ''')

args = p.parse_args()

# load cosmic rays dataset
cr_file = args.cr_file
# read flash.dat to a list of lists
datContent = [i.strip().split() for i in open(cr_file).readlines()]

dataset = pd.DataFrame(datContent[1:], columns=datContent[0])

for col in datContent[0]:
    dataset[col] = dataset[col].astype('float') 

mask = dataset['Th'] >= 60
# print(f"The provided list contains {n_sources_initial} objects.")

cr_vertical = dataset[~mask]
cr_inclined = dataset[mask]

# apply energy cut
min_energy = args.min_energy_dipole
mask_energy_vertical = cr_vertical["energy"] > min_energy
mask_energy_inclined = cr_inclined["energy"] > min_energy
cr_vertical = cr_vertical[mask_energy_vertical]
cr_inclined = cr_inclined[mask_energy_inclined]

min_energy_analysis = args.min_energy_analysis
mask_energy_vertical_analysis = cr_vertical["energy"] > min_energy_analysis
mask_energy_inclined_analysis = cr_inclined["energy"] > min_energy_analysis
# n_vertical_events_final = len(cr_vertical[mask_energy_vertical_analysis])
# n_inclined_events_final = len(cr_inclined[mask_energy_inclined_analysis])

#dipole information
d, alpha_d, delta_d = 0.17, np.radians(144), np.radians(-51)

# define the coordinates grid
ra_bins, dec_bins = np.linspace(0,2*np.pi, 1000), np.linspace(-1, 1, 1000)
ra_grid, dec_grid = np.meshgrid(ra_bins, dec_bins)

#define exposure pdf
exposure_vert, exposure_incl, exposure_tot = auger_tools.LoadExposureMap(60, 80, np.arcsin(dec_grid))
exposure_vert_smooth = auger_tools.smooth_flux(exposure_vert, sigma = 15)
exposure_incl_smooth = auger_tools.smooth_flux(exposure_incl, sigma = 10)

#define dipole pdf
flux_values = auger_tools.dipole_flux(ra_grid, np.arcsin(dec_grid), d,  alpha_d, delta_d)

#define final pdf
pdf_vertical = flux_values*exposure_vert
pdf_inclined = flux_values*exposure_incl

#define seeds (need 2 for each simulation because we don't want r.a. and dec to be correlated)
seeds = np.arange(args.seed_initial, args.seed_final, 2)
results_dtype = [('seed', '<i4'),
                      ('ra',list ),
                      ('dec',list )]

results = np.zeros(len(seeds), dtype = results_dtype)

mode = args.mode

if mode == 'sample':
    print('Simulating the events in the public sample')
    n_events_vertical, n_events_inclined = len(cr_vertical), len(cr_inclined)
    n_vertical_events_final = len(cr_vertical[mask_energy_vertical_analysis])
    n_inclined_events_final = len(cr_inclined[mask_energy_inclined_analysis])

if mode == 'all_statistics':
    print('Producing a simulation of the full sample of events')
    n_events_vertical, n_events_inclined = (len(cr_vertical)-109)*10+109, len(cr_inclined)*10
    n_vertical_events_final = (len(cr_vertical[mask_energy_vertical_analysis])-109)*10 + 109
    n_inclined_events_final = len(cr_inclined[mask_energy_inclined_analysis])*10

for i, seed in enumerate(seeds):
    sim_ra_vert, sim_sindec_vert, sim_ra_incl, sim_sindec_incl = simulation.do_simulation(pdf_vertical, pdf_inclined, n_events_vertical, n_events_inclined, ra_bins, dec_bins, seed, n_vertical_events_final, n_inclined_events_final)
    sim_ra = np.concatenate([sim_ra_vert, sim_ra_incl])
    sim_sindec = np.concatenate([sim_sindec_vert, sim_sindec_incl])

    sim_dec = np.arcsin(sim_sindec)

    results['seed'][i] = seed
    results['ra'][i] = sim_ra
    results['dec'][i] = sim_dec

#save trials
outfilename = f'trials_minEnergy_{min_energy_analysis}EeV_seed_in_{args.seed_initial}_seed_fin_{args.seed_final}.npy'

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