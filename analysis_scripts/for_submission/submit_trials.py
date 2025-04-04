import numpy as np
import subprocess
import multiprocessing as mp
import argparse

p = argparse.ArgumentParser(description="Input parameters producing trials for correaltion analysis.")

p.add_argument("--pao_hotspot_treatment", 
               type = str,  
               help = '''How to mask the PAO hotspot. Here are the options:
               - no_mask: don't mask the hotspot
               - mask_and_scramble: first mask the hotspot and then scramble, making sure that no sources fall in the hotspot
               - scramble_and_mask_with_candidate_sources: first scramble and then remove sources that fall in the PAO hotspot
               - scramble_and_mask_without_candidate_sources: remove the souurces that are believed to be CR emitters, and then do same as scramble_and_mask_with_candidate_sources  '''
                )
p.add_argument("--min_energy_analysis", 
               type = float,  
               help = '''Minimum energy of the cosmic rays (in EeV)''')
p.add_argument("--seed_initial",
               type = int, 
               default = 0, 
               help = 'Number of trials we want to perform.')
p.add_argument("--seed_final",
               type = int, 
               default = 1000, 
               help = 'Number of trials we want to perform.')

args = p.parse_args()


def run_script(seed):
    subprocess.run([
        "python", "../smbh/analysis_scripts/run_correlation.py",
        "--sourcelist", "../sources_coords.h5",
        "--cr_vertical", "../data/cr_vertical_events.npy",
        "--cr_inclined", "../data/cr_inclined_events.npy",
        "--outdir", "./trials_simulation",
        "--pao_hotspot_treatment", args.pao_hotspot_treatment,
        "--min_energy_analysis", str(args.min_energy_analysis),
        "--seed_initial", str(seed_initial),
        "--seed_final", str(seed_final)
    ])
if __name__ == "__main__":
    seed_initial = args.seed_initial  # Test with just two seeds
    seed_final = args.seed_final
    num_workers = mp.cpu_count()  # Adjust based on your system

    with mp.Pool(num_workers) as pool:
        pool.map(run_script, np.arange(seed_initial, 2*seed_final, 2))
