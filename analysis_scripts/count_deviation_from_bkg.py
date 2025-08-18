import numpy as np
import sys
import pandas as pd

sys.path.append('tools')

from tools import correlation

#the faxct is that I need to compare the number of events with the background expectation, because the event rate depends on the coordinate of the source.
#load trials 
trials_file = 'simulated_sky_tobereproduced/trials_minEnergy_20EeV_seed_in_0_seed_fin_100000.npy'
trials = np.load(trials_file, allow_pickle = True)
radius = 23

dec = np.radians(np.linspace(-90, 44.8, 15))
ra = np.radians(np.linspace(0,360, 15))

dict_trials = {'dec': [],
               'ra':[],
            #    'counts': [],
               'n_events': [],
               'quantile16':[],
               'quantile84':[],
               'quantile3':[],
               'quantile97':[],
               'quantile0.1':[],
               'quantile99.9':[]}
for d in dec:
    for r in ra:
        print('Running for coordinates', d, r)
        counts_array = []
        for s, seed in enumerate(trials['seed']):
            counts, ra, dec = correlation.count_events_per_source(r, d, trials['ra'][s], trials['dec'][s], np.radians(radius))
            counts_array.append(counts)
        dict_trials['ra'] = np.append(dict_trials['ra'], r)
        dict_trials['dec'] = np.append(dict_trials['dec'], d)
        # dict_trials['counts'] = np.append(dict_trials['counts'], counts_array)
        dict_trials['n_events'] = np.append(dict_trials['n_events'], np.median(counts_array))
        dict_trials['quantile16'] = np.append(dict_trials['quantile16'], np.quantile(counts_array, 0.16))
        dict_trials['quantile84'] = np.append(dict_trials['quantile84'], np.quantile(counts_array, 0.84))
        dict_trials['quantile3'] = np.append(dict_trials['quantile3'], np.quantile(counts_array, 0.03))
        dict_trials['quantile97'] = np.append(dict_trials['quantile97'], np.quantile(counts_array, 0.97))
        dict_trials['quantile0.1'] = np.append(dict_trials['quantile0.1'], np.quantile(counts_array, 0.001))
        dict_trials['quantile99.9'] = np.append(dict_trials['quantile99.9'], np.quantile(counts_array, 0.999))
    

#save the trials
df_trials = pd.DataFrame(dict_trials)
df_trials.to_csv('trials.csv', index = False)


