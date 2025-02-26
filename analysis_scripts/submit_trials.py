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
p.add_argument("--n_seeds",
               type = int, 
               default = 1000, 
               help = 'Number of trials we want to perform.')

args = p.parse_args()


def run_script(seed):
    subprocess.run([
        "python", "run_correlation.py",
        "--sourcelist", "../sources_coords.h5",
        "--cr_vertical", "../data/cr_vertical_events.npy",
        "--cr_inclined", "../data/cr_inclined_events.npy",
        "--outdir", "./trials_simulation",
        "--pao_hotspot_treatment", args.pao_hotspot_treatment,
        "--min_energy_analysis", str(args.min_energy_analysis),
        "--seed", str(seed)
    ])
if __name__ == "__main__":
    num_seeds = args.n_seeds  # Test with just two seeds
    num_workers = mp.cpu_count()  # Adjust based on your system

    with mp.Pool(num_workers) as pool:
        pool.map(run_script, np.arange(0, 2*num_seeds, 2))
