import numpy as np
import pandas as pd
import sys
import os
import argparse

# Load library for analysis
sys.path.append('tools')
from tools import simulation, correlation, auger_tools

p = argparse.ArgumentParser(description="Input parameters producing trials for correlation analysis.")

p.add_argument("--cr_vertical", type=str, help="Path to cosmic rays vertical dataset (npy file)")
p.add_argument("--cr_inclined", type=str, help="Path to cosmic rays inclined dataset (npy file)")
p.add_argument("--outdir", type=str, help="Path to the directory where to save trials")
p.add_argument("--filter_sources", type=str, default='all_sources', help="Which sources to analyse (all_sources, AGN, no_AGN)")
p.add_argument("--mask_declination", type=bool, default=True, help="Remove sources not visible by Pierre Auger Observatory")
p.add_argument('--min_energy_analysis', type=float, default=20, help='Minimum energy for the analysis')  
p.add_argument('--min_energy_dipole', type=float, default=8, help='Minimum energy for the simulation of the events with the dipole')  
p.add_argument("--seed_initial", type=int, default=1, help='First seed to use')
p.add_argument("--seed_final", type=int, default=100, help='Last seed to use')
p.add_argument('--mode', type=str, default='sample', help="sample: simulate public events, all_statistics: forecast full sample")
p.add_argument('--n_signal', type=int, default=0, help="Number of signal events to be distributed among sources")
p.add_argument('--source_list', type=str, help="Path to the .npy or .csv file containing the source list with RA/Dec in radians")
p.add_argument('--mask_source', type=str, help="AGN or non-AGN")


args = p.parse_args()

# Load cosmic rays dataset
cr_vertical = np.load(args.cr_vertical, allow_pickle=True)
cr_inclined = np.load(args.cr_inclined, allow_pickle=True)

# Apply energy cuts
min_energy = args.min_energy_dipole
mask_energy_vertical = cr_vertical["energy"] > min_energy
mask_energy_inclined = cr_inclined["energy"] > min_energy
cr_vertical = cr_vertical[mask_energy_vertical]
cr_inclined = cr_inclined[mask_energy_inclined]

min_energy_analysis = args.min_energy_analysis
mask_energy_vertical_analysis = cr_vertical["energy"] > min_energy_analysis
mask_energy_inclined_analysis = cr_inclined["energy"] > min_energy_analysis

# Load source list if provided
sources = pd.read_hdf(args.source_list, key = 'values')
mask_dec = sources['DEC_deg'] <44.8
sources = sources[mask_dec]

if args.mask_source == 'AGN':
    print(len(sources))
    mask = sources['Activity'] == True
    sources = sources[mask]
    print(f"Number of sources: {len(sources)}")
if args.mask_source == 'non_AGN':
    mask = sources['Activity'] == True
    sources = sources[~mask]
    print(f"Number of sources: {len(sources)}")


# Dipole parameters
d, alpha_d, delta_d = 0.074, np.radians(97), np.radians(-38)

# Define coordinate grid
ra_bins, dec_bins = np.linspace(0, 2 * np.pi, 1000), np.linspace(-1, 1, 1000)
ra_grid, dec_grid = np.meshgrid(ra_bins, dec_bins)

# Define exposure and flux maps
exposure_vert, exposure_incl, _ = auger_tools.LoadExposureMap(60, 80, np.arcsin(dec_grid))
exposure_vert = auger_tools.smooth_flux(exposure_vert, sigma=15)
exposure_incl = auger_tools.smooth_flux(exposure_incl, sigma=10)
flux_values = auger_tools.dipole_flux(ra_grid, np.arcsin(dec_grid), d, alpha_d, delta_d)

pdf_vertical = flux_values * exposure_vert
pdf_inclined = flux_values * exposure_incl

# Seeds and output structure
seeds = np.arange(args.seed_initial, args.seed_final, 2)
results_dtype = [('seed', '<i4'), ('ra', list), ('dec', list)]
results = np.zeros(len(seeds), dtype=results_dtype)

# Event counts
mode = args.mode
if mode == 'sample':
    print('Simulating the events in the public sample')
    n_events_vertical = len(cr_vertical)
    n_events_inclined = len(cr_inclined)
    n_vertical_events_final = len(cr_vertical[mask_energy_vertical_analysis])
    n_inclined_events_final = len(cr_inclined[mask_energy_inclined_analysis])
elif mode == 'all_statistics':
    print('Producing a simulation of the full sample of events')
    n_events_vertical = (len(cr_vertical) - 109) * 10 + 109
    n_events_inclined = len(cr_inclined) * 10
    n_vertical_events_final = (len(cr_vertical[mask_energy_vertical_analysis]) - 109) * 10 + 109
    n_inclined_events_final = len(cr_inclined[mask_energy_inclined_analysis]) * 10
else:
    raise ValueError(f"Unsupported mode: {mode}")

# Simulation loop
for i, seed in enumerate(seeds):
    sim_ra_vert, sim_sindec_vert, sim_ra_incl, sim_sindec_incl = simulation.do_simulation(
        pdf_vertical, pdf_inclined, n_events_vertical, n_events_inclined,
        ra_bins, dec_bins, seed, n_vertical_events_final, n_inclined_events_final
    )
    sim_ra = np.concatenate([sim_ra_vert, sim_ra_incl])
    sim_sindec = np.concatenate([sim_sindec_vert, sim_sindec_incl])
    sim_dec = np.arcsin(sim_sindec)

    # Remove random background events and inject signals if requested
    if args.n_signal > 0:
        rng = np.random.default_rng(seed)
        total_events = len(sim_ra)
        if args.n_signal > total_events:
            raise ValueError("n_signal exceeds the number of total simulated events")

        mask_keep = np.ones(total_events, dtype=bool)
        indices_to_replace = rng.choice(total_events, size=args.n_signal, replace=False)
        mask_keep[indices_to_replace] = False

        sim_ra = sim_ra[mask_keep]
        sim_dec = sim_dec[mask_keep]

        if sources is None:
            raise ValueError("To inject signal events, you must provide a source list via --source_list")

        # Inject clustered signal events
        chosen_sources = rng.choice(sources, size=args.n_signal)
        sigma_rad = np.deg2rad(30)
        signal_ra, signal_dec = [], []

        for s in chosen_sources:
            d_ra = rng.normal(0, sigma_rad)
            d_dec = rng.normal(0, sigma_rad)

            ra_sig = (np.radians(s[5]) + d_ra) % (2 * np.pi)
            dec_sig = np.clip(np.radians(s[6]) + d_dec, -np.pi / 2, np.radians(44.8))

            signal_ra.append(ra_sig)
            signal_dec.append(dec_sig)

        sim_ra = np.concatenate([sim_ra, signal_ra])
        sim_dec = np.concatenate([sim_dec, signal_dec])

    results['seed'][i] = seed
    results['ra'][i] = sim_ra
    results['dec'][i] = sim_dec

# Save output
outfilename = f'trials_{args.mask_source}_minEnergy_{min_energy_analysis}EeV_seed_in_{args.seed_initial}_seed_fin_{args.seed_final}_{args.n_signal}_signalevents.npy'
outdir = args.outdir

if not os.path.exists(outdir):
    os.makedirs(outdir)
    print("Directory created successfully!")
else:
    print("Directory already exists!")

out_path = os.path.join(outdir, outfilename)
print('Saving file as: ', out_path)
np.save(out_path, results)
