from ast import Str
import numpy as np
import sys 
import matplotlib.pyplot as plt
import argparse
sys.path.append('tools')

from analysis_scripts.tools import statistics_utils

p = argparse.ArgumentParser(description="Input parameters producing trials for correaltion analysis.")

p.add_argument("--trials_files_AGN", 
               type = str,
               help = 'Input three values, the starting radius, the final radius and the step')
p.add_argument("--trials_files_nonAGN", 
               type = str,
               help = 'Input three values, the starting radius, the final radius and the step')
p.add_argument("--trials_files_all", 
               type = str,
               help = 'Input three values, the starting radius, the final radius and the step')
p.add_argument("--steps", 
               type = list,
               nargs='+',
               default = [1,31, 1],
               help = 'Input three values, the starting radius, the final radius and the step')
p.add_argument("--pvalue", 
               type=float, 
               help="pvalue to be corrected")
p.add_argument("--ntrials", 
               type=int, 
               help="pvalue to be corrected")


args = p.parse_args()
n = args.ntrials
trial_files = [args.trials_files_AGN, args.trials_files_nonAGN, args.trials_files_all]

pval_distribution_dtype = [('seed', int)]
for i, files in enumerate(trial_files):
    pval_distribution_dtype.append((f'pvalue_{i+1}', float)) # type: ignore

pval_distribution = np.zeros(n, dtype=pval_distribution_dtype)
steps = np.arange(args.steps[0], args.steps[1], args.steps[2])
pretrial_pvalues = args.pvalue

for f, trialfile in enumerate(trial_files):
    trials = np.load(trialfile, allow_pickle=True)
    #for each trial, calculate their p-value for all the cases, and see how often we get a p-value smaller than the one we have
    for i in np.arange(n):
        tmp_pval = []
        for step in steps:
            pval = statistics_utils.calculate_pvalue_from_trials(trials[f'fraction_{step}'],trials[f'fraction_{step}'][i])
            tmp_pval.append(pval)
    
        pval_distribution[f'pvalue_{f+1}'][i] = np.amin(tmp_pval)

#calculate the trial correction
pval_final = np.ones(n)
for i in np.arange(n):
    pval_final[i] = np.amin([pval_distribution[f'pvalue_{f+1}'][i] for f in np.arange(len(trial_files))])



trial_correction = np.sum(pval_final < pretrial_pvalues)/len(pval_final)


print(f'The trial correction for the p-value {pretrial_pvalues} is {trial_correction}, resulting in correction factor of {trial_correction/pretrial_pvalues}')



