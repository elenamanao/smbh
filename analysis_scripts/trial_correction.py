import numpy as np
import sys 
import matplotlib.pyplot as plt
sys.path.append('tools')

from tools import statistics

import argparse

p = argparse.ArgumentParser(description="Input parameters producing trials for correaltion analysis.")

p.add_argument("--trials", 
               type=str, 
               help="Path to trials to calculate trial correction")
p.add_argument("--steps", 
               type = list,
               nargs='+',
               default = [1,31, 1],
               help = 'Input three values, the starting radius, the final radius and the step')
p.add_argument("--pvalue", 
               type=float, 
               help="pvalue to be corrected")

args = p.parse_args()

trial_file = args.trials
trials = np.load(trial_file, allow_pickle=True)
steps = np.arange(args.steps[0], args.steps[1], args.steps[2])
pretrial_pvalues = args.pvalue

#for each trial, calculate their p-value for all the cases, and see how often we get a p-value smaller than the one we have
pval_distribution = np.ones(5000)
for i in np.arange(5000):
    tmp_pval = []
    for step in steps:
        pval = statistics.calculate_pvalue_from_trials(trials[f'fraction_{step}'],trials[f'fraction_{step}'][i])
        tmp_pval.append(pval)
    
    pval_distribution[i] = np.amin(tmp_pval)

#calculate the trial correction

trial_correction = np.sum(pval_distribution < pretrial_pvalues)/len(pval_distribution)


print(f'The trial correction for the p-value {pretrial_pvalues} is {trial_correction}, resulting in correction factor of {trial_correction/pretrial_pvalues}')



