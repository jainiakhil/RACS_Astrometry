import os
import sys
import time
import socket
import pickle

import numpy as np
import pandas as pd
import multiprocessing as mp

from astropy.coordinates import SkyCoord 
from astropy import units as u
from astropy.table import Table, vstack

from astroquery.utils.tap.core import Tap, TapPlus
from astroquery.vizier import Vizier


def query_racs(ra, dec, radius):
    success = False
    while not success:
        try:
            tap = TapPlus(url="https://casda.csiro.au/casda_vo_tools/tap")
            job = tap.launch_job_async(f"SELECT * FROM AS110.racs_dr1_gaussians_galacticcut_v2021_08_v02 WHERE 1=CONTAINS(POINT('ICRS', ra, dec),CIRCLE('ICRS', {ra},{dec},{radius}))")
            if len(job.get_results()) == 0:
                job = tap.launch_job_async(f"SELECT * FROM AS110.racs_dr1_gaussians_galacticregion_v2021_08_v02 WHERE 1=CONTAINS(POINT('ICRS', ra, dec),CIRCLE('ICRS', {ra},{dec},{radius}))")
            success = True
        except (TimeoutError, ConnectionError, socket.gaierror, ConnectionRefusedError, ConnectionAbortedError) as e: 
            print(f'Retrying again in 10 seconds...')
            sys.stdout.flush()
            time.sleep(10)

    return job.get_results()


def query_racs_gal(ra, dec, radius):
    success = False
    retries = 0  # Counter for retries
    while not success and retries < 5:  # Retry only for 5 times
        try:
            tap = TapPlus(url="https://casda.csiro.au/casda_vo_tools/tap")
            job1 = tap.launch_job_async(f"SELECT * FROM AS110.racs_dr1_gaussians_galacticcut_v2021_08_v02 WHERE 1=CONTAINS(POINT('ICRS', ra, dec),CIRCLE('ICRS', {ra},{dec},{radius}))")
            job2 = tap.launch_job_async(f"SELECT * FROM AS110.racs_dr1_gaussians_galacticregion_v2021_08_v02 WHERE 1=CONTAINS(POINT('ICRS', ra, dec),CIRCLE('ICRS', {ra},{dec},{radius}))")
            success = True
        except (TimeoutError, ConnectionError, socket.gaierror, ConnectionRefusedError, ConnectionAbortedError) as e: 
            print(f'Retrying again in 10 seconds...')
            sys.stdout.flush()
            time.sleep(10)
            retries += 1  # Increment the counter
    if retries == 5:
        job1 = Table(names=['ra', 'dec', 'peak_flux', 'int_flux', 'err_peak_flux', 'err_int_flux', 'maj_axis', 'min_axis', 'pa', 'rms', 'catalog'])
        job2 = Table(names=['ra', 'dec', 'peak_flux', 'int_flux', 'err_peak_flux', 'err_int_flux', 'maj_axis', 'min_axis', 'pa', 'rms', 'catalog'])

    return vstack([job1.get_results(), job2.get_results()])


def query_racs_mid(ra, dec, radius):
    success = False
    retries = 0  # Counter for retries
    while not success and retries < 5:  # Retry only for 5 times
        try:
            tap = TapPlus(url="https://casda.csiro.au/casda_vo_tools/tap")
            job = tap.launch_job_async(f"SELECT * FROM AS110.racs_mid_components_v01 WHERE 1=CONTAINS(POINT('ICRS', ra, dec),CIRCLE('ICRS', {ra},{dec},{radius}))")
            success = True
        except (TimeoutError, ConnectionError, socket.gaierror, ConnectionRefusedError, ConnectionAbortedError) as e: 
            print(f'Retrying again in 10 seconds...')
            sys.stdout.flush()
            time.sleep(10)
            retries += 1  # Increment the counter
    if retries == 5:
        job = Table(names=['ra', 'dec', 'peak_flux', 'int_flux', 'err_peak_flux', 'err_int_flux', 'maj_axis', 'min_axis', 'pa', 'rms', 'catalog'])
        
    return job.get_results()


def query_wise(ra, dec, radius):
    coord = SkyCoord(ra, dec, unit=u.degree)
    success = False
    retries = 0  # Counter for retries
    while not success and retries < 5:  # Retry only for 5 times
        try:
            v = Vizier(columns=['AllWISE', 'RAJ2000', 'DEJ2000', 'eeMaj', 'eeMin', 'eePA'], row_limit=-1)
            job = v.query_region(coord, radius=radius*u.degree, catalog='II/328')[0]
            success = True
        except (TimeoutError, ConnectionError, socket.gaierror, ConnectionRefusedError, ConnectionAbortedError, pickle.UnpicklingError) as e: 
            print(f'Retrying again in 10 seconds...')
            sys.stdout.flush()
            time.sleep(10)
            retries += 1  # Increment the counter
    if retries == 5:
        job = Table(names=['AllWISE', 'RAJ2000', 'DEJ2000', 'eeMaj', 'eeMin', 'eePA'], dtype=['str', 'f8', 'f8', 'f8', 'f8', 'f4'])
    
    return job


def query_first(ra, dec, radius):
    coord = SkyCoord(ra, dec, unit=u.degree)
    success = False
    retries = 0  # Counter for retries
    while not success and retries < 5:  # Retry only for 5 times
        try:
            v = Vizier(columns=['FIRST', 'RAJ2000', 'DEJ2000', 'fMaj', 'fMin', 'Fpeak', 'Fint', 'Rms'], row_limit=-1)
            job = v.query_region(coord, radius=radius*u.degree, catalog='VIII/92/first14')[0]
            success = True
        except (TimeoutError, ConnectionError, socket.gaierror, ConnectionRefusedError, ConnectionAbortedError, IndexError) as e: 
            print(f'Retrying again in 10 seconds...')
            sys.stdout.flush()
            time.sleep(10)
            retries += 1  # Increment the counter
    if retries == 5:
        job = Table()
    
    return job


def query_vlass(ra, dec, radius):
    coord = SkyCoord(ra, dec, unit=u.degree)
    success = False
    retries = 0  # Counter for retries
    while not success and retries < 5:  # Retry only for 5 times
        try:
            v = Vizier(columns=['CompName', 'RAJ2000', 'DEJ2000', 'e_RAJ2000', 'e_DEJ2000', 'Ftot', 'e_Ftot', 'Fpeak', 'e_Fpeak', 'Maj', 'e_Maj', 'Min', 'e_Min', 'PA', 'e_PA'], row_limit=-1)
            job = v.query_region(coord, radius=radius*u.degree, catalog='VLASS')[0]
            success = True
        except (TimeoutError, ConnectionError, socket.gaierror, ConnectionRefusedError, ConnectionAbortedError, IndexError) as e: 
            print(f'Retrying again in 10 seconds...')
            sys.stdout.flush()
            time.sleep(10)
            retries += 1  # Increment the counter
    if retries == 5:
        job = Table()
    
    return job


def VLASS_lookup_local(ra, dec, cat, radius=1.0, ra_col = 'RA', dec_col = 'DEC'):
    """
    Perform a VLASS lookup based on the given coordinates and catalog.

    Parameters:
    ra (float): Right Ascension in degrees.
    dec (float): Declination in degrees.
    cat (pandas.DataFrame): Catalog containing RA and DEC columns.
    radius (float, optional): Search radius in degrees. Default is 1.0 degrees.
    ra_col (str, optional): RA column name in the catalog. Default is 'RA'.
    dec_col (str, optional): DEC column name in the catalog. Default is 'DEC'.

    Returns:
    pandas.DataFrame: Subset of the catalog containing sources within the specified radius.
    """
    cat_ra, cat_dec = np.array(cat[ra_col]), np.array(cat[dec_col])
    
    phi1 = ra * np.pi / 180
    theta1 = dec * np.pi / 180
    phi2 = cat_ra * np.pi / 180
    theta2 = cat_dec * np.pi / 180
    
    cos_sep_radian = np.sin(theta1) * np.sin(theta2) + np.cos(theta1) * np.cos(theta2) * np.cos(phi1-phi2)
    sep = np.arccos(cos_sep_radian) * 180 / np.pi
    
    select_bool = sep < radius
    
    return cat.iloc[select_bool]


def query_local_cat(ra, dec, cat, radius=1.0, ra_col = 'ra', dec_col = 'dec'):
    """
    Query any local catalog based on the given coordinates and radius.

    Parameters:
    ra (float): Right Ascension in degrees.
    dec (float): Declination in degrees.
    cat (pandas.DataFrame): Catalog containing RA and DEC columns.
    radius (float, optional): Search radius in degrees. Default is 1.0 degree.
    ra_col (str, optional): RA column name in the catalog. Default is 'ra'.
    dec_col (str, optional): DEC column name in the catalog. Default is 'dec'.

    Returns:
    pandas.DataFrame: Subset of the catalog containing sources within the specified radius.
    """
    cat_ra, cat_dec = np.array(cat[ra_col]), np.array(cat[dec_col])
    
    phi1 = ra * np.pi / 180
    theta1 = dec * np.pi / 180
    phi2 = cat_ra * np.pi / 180
    theta2 = cat_dec * np.pi / 180
    
    cos_sep_radian = np.sin(theta1) * np.sin(theta2) + np.cos(theta1) * np.cos(theta2) * np.cos(phi1-phi2)
    sep = np.arccos(cos_sep_radian) * 180 / np.pi
    
    select_bool = sep < radius
    
    return cat.iloc[select_bool]