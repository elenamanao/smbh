'''Library with functions for scrambling sources in the Sky.
All the functions take the coordinates of the sources from a list
and perform a permutation of the r.a. of the array with a given seed.'''

import numpy as np

def GreatCircleDistance(ra_1, dec_1, ra_2, dec_2, unit="rad"):
    """Compute the great circle distance between 
    two points in the Sky"""
    if unit == "deg":
        ra_1 = np.radians(ra_1)
        ra_2 = np.radians(ra_2)
        dec_1 = np.radians(dec_1)
        dec_2 = np.radians(dec_2)
    delta_dec = np.abs(dec_1 - dec_2)
    delta_ra = np.abs(ra_1 - ra_2)
    x = (np.sin(delta_dec / 2.0)) ** 2.0 + np.cos(dec_1) * np.cos(dec_2) * (
        np.sin(delta_ra / 2.0)) ** 2.0
    return 2.0 * np.arcsin(np.sqrt(x))

def scramble_ra(ra_true, seed):
    """Scramble the r.a. with a given seed"""
    np.random.RandomState(seed)
    ra_rad_tmp = np.random.permutation(ra_true)

    return ra_rad_tmp

def mask_hotspot_and_scramble(ra_true, dec_true, seed, hs_ra, hs_dec, hs_radius):
    '''First mask the hotspot, 
    then scramble and make sure that'
    no source falls inside the PAO hotspot 
    after scrambling'''

    #Removing the sources that fall inside PAO hotspot
    gcd = GreatCircleDistance(ra_true, 
                              dec_true, 
                              np.ones_like(ra_true)*hs_ra, 
                              np.ones_like(ra_true)*hs_dec,
                              unit = 'deg'
                              )
    mask = gcd > np.deg2rad(hs_radius)
    ra_true = ra_true[mask]
    dec_true = dec_true[mask]

    #now the scrambling begins
    i = 1
    seed *=10000
    while i > 0:
        #first just scramble the remaining sources
        ra_tmp = scramble_ra(ra_true, seed)
        #check that no sources fall into PAO hotspot after permutation
        gcd = GreatCircleDistance(ra_tmp, 
                                  dec_true,
                                  np.ones_like(ra_tmp)*hs_ra, 
                                  np.ones_like(dec_true)*hs_dec, 
                                  unit = 'deg')
        check_distance = gcd > np.deg2rad(hs_radius)
        #repeat until no sources are inside the hotspot.
        if np.all(check_distance):
            i = 0
            return ra_tmp, ra_true, dec_true
        else:
            seed += i

def scramble_and_mask_hotspot(ra_true, dec_true, seed, hs_ra, hs_dec, hs_radius):
    '''We scramble the sources and subsequently we 
    remove all the sources that fall into the PAO hotspot'''
    #scrambling
    ra_tmp = scramble_ra(ra_true, seed)
    #check sources fall into PAO hotspot after permutation
    gcd = GreatCircleDistance(ra_tmp, 
                              dec_true,
                              np.ones_like(ra_tmp)*hs_ra,
                              np.ones_like(dec_true)*hs_dec, 
                              unit = 'deg')
    #remove those
    check_distance = gcd > np.deg2rad(hs_radius)
    ra_tmp = ra_tmp[check_distance]
    dec_tmp = dec_true[check_distance]

    return ra_tmp, dec_tmp
            
            