#!/usr/bin/env python

import numpy as np
import os
from astropy.io import fits

def choose_dark(exptime, dark_dir='darks', suffix='medcombine_sub',
                explist=[5, 10, 20, 30, 60, 120, 300], verbose=False):
    """Chooses or interpolates a dark file based on exposure time from a set 
    of darks.

    Assumes: integer exposure times of darks, naming schema of
        'dark_<exptime>s_<suffix>.fits'

    Parameters
    ----------
    exptime : numeric
        Exposure time of file to be dark-subtracted.
    dark_dir : string
        Directory with library darks. Default: 'darks'
    suffix : string
        Suffix of dark file names. Default: 'medcombine_sub', indicating
        median-combined and bias-subtracted darks
    explist : list of ints
        List of available dark exposure times.
    verbose : boolean
        Whether to print messages indicating choice of dark. Default: False

    Returns
    -------
    dark : 2d array
        Image data of chosen or interpolated dark files
    """
    # assumes: integer exposure times, naming schema of 'dark_<exptime>s_<suffix>.fits'
    explist = np.array(explist)
    if int(exptime) in explist:
        darkfile = os.path.join(dark_dir, 'dark_{}s_{}.fits'.format(int(exptime), suffix))
        if verbose:
            print('Using', darkfile)
        dark = fits.getdata(darkfile, 0)
    else:
        if exptime < explist.min():
            # scale shortest dark
            # if you're in this regime you probably don't need a dark anyway but you do you
            darkfile = os.path.join(dark_dir, 'dark_{}s_{}.fits'.format(int(explist.min()), suffix))
            dark = fits.getdata(darkfile, 0) * exptime/explist.min()
        elif (exptime > explist.min()) & (exptime < explist.max()):
            # calculate weighted average of two nearest darks 
            lower_bound = explist[explist < exptime].max()
            upper_bound = explist[explist > exptime].min()
            diff = upper_bound - lower_bound
            weight_lower = (exptime - lower_bound)/diff
            weight_upper = (upper_bound - exptime)/diff
            darkfile1 = os.path.join(dark_dir, 'dark_{}s_{}.fits'.format(int(lower_bound), suffix))
            darkfile2 = os.path.join(dark_dir, 'dark_{}s_{}.fits'.format(int(upper_bound), suffix))
            if verbose:
                print('Averaging', darkfile1, '&', darkfile2)
            dark1 = fits.getdata(darkfile1)
            dark2 = fits.getdata(darkfile2)
            dark = dark1 * weight_lower + dark2 * weight_upper
        elif exptime > explist.max():
            # scale longest dark
            darkfile = os.path.join(dark_dir, 'dark_{}s_{}.fits'.format(int(explist.max()), suffix))
            if verbose:
                print('Using', darkfile)
            dark = fits.getdata(darkfile, 0) * exptime/explist.max()
    return dark