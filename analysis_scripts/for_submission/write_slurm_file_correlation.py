from matplotlib.pylab import f
import numpy as np
import glob
import argparse

p = argparse.ArgumentParser(description="Input parameters producing trials for correaltion analysis.")

p.add_argument("--file_dir", 
               type = str,  
               help = '''How to mask the PAO hotspot. Here are the options:
               - no_mask: don't mask the hotspot
               - mask_and_scramble: first mask the hotspot and then scramble, making sure that no sources fall in the hotspot
               - scramble_and_mask_with_candidate_sources: first scramble and then remove sources that fall in the PAO hotspot
               - scramble_and_mask_without_candidate_sources: remove the souurces that are believed to be CR emitters, and then do same as scramble_and_mask_with_candidate_sources  '''
                )
p.add_argument("--filter_sources",
               type = str,
               default = 'all_sources',
               help = ''''Which sources do you want to analyse:
               - all_sources: analyse all the sources in the list
               - AGN: take only AGN in the list
               - no_AGN: take only non-AGN in the list''')
p.add_argument("--slurm_basename",
               type = str,
               default = 'correlation_trials',
               help = ''''Name of slurm file''')

args = p.parse_args()


files = glob.glob(args.file_dir + "/*")

print(f"The trials are distributed over {len(files)} files")

if len(files) > 25:
    print('Every slurm job will run the correlation for 25 files, so we need to produce more than one slurm.')

n_slurm = len(files)//25 +1

slurmnames = [f"{args.slurm_basename}_{i}.slurm" for i in range(n_slurm)]

files_per_slurm = np.array_split(files, n_slurm)

filter_sources = args.filter_sources

for slurm, files in zip(slurmnames, files_per_slurm):
    with open(slurm, "w") as fout:
        fout.write(f'''#!/bin/bash
#SBATCH -J my_job_farm
#SBATCH -o /dss/dsshome1/00/ge74xut2/submission/out/output_%A_%a.out
#SBATCH -e /dss/dsshome1/00/ge74xut2/submission/err/error_%A_%a.err
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks=1            
#SBATCH --cpus-per-task=24
#SBATCH --mem=30G       
#SBATCH --time=08:00:00
#SBATCH --array=0-{len(files)-1}%4
module load slurm_setup

# Define file list
FILES=({ " ".join(files) })

FILTER={filter_sources}

# Get file for this task
INPUT_FILE=${{FILES[$SLURM_ARRAY_TASK_ID]}}

# Run the script with the selected file
./run_correlation $INPUT_FILE $FILTER
''')



