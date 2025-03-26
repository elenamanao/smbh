These are the scripts to run the correlation analysis

First, we simulate a number of skymaps using `run_simulation.py`

This will produce a numpy file with the simulated CR sky. Here is an example on how to run it

```python run_simulation.py  --cr_vertical ../data/cr_vertical_events.npy --cr_inclined ../data/cr_inclined_events.npy --outdir ./simulated_sky --seed_initial 0 --seed_final 1000  ```

Afterwards, you will need to evaluate the correlation, this is done with the script `run_correlation.py`. Here's an example on how to run it:

```python run_correlation.py  --sourcelist ../all_sources_coords.h5 --sky_map_file ./simulated_sky/trials_no_mask_minEnergy_20EeV_seed_in_0_seed_fin_1000.npy --outdir corr_trials --pao_hotspot_treatment no_mask --filter_sources AGN```

This will have to be run also for the `non_AGN` case and the `all_sources` case, the only difference will the the `filter_sources` argument.