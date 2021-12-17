# third party
import numpy as np

import scipy
from scipy.ndimage.interpolation import shift

import matplotlib.pyplot as plt

import astropy
import astropy.io.fits as fits
from astropy.cosmology import FlatLambdaCDM, Planck15
import astropy.units as u
from astropy.units import Quantity
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

import astroquery
from astroquery.gaia import Gaia
from astroquery.mast import Observations
from astroquery.sdss import SDSS

#import drizzlepac
#from drizzlepac import tweakreg

import json, requests

from reproject import reproject_interp

import os

import aplpy

# internal
from search import *
from align import *
from plotting import *

# ---------------------------------------------------------------------------------------------------
# FIND THE OBJECT ON MAST
# ---------------------------------------------------------------------------------------------------

def query_object(name, radius, **kwargs):
    '''
    Queries an object within a specified catalog and creates an astroquery Catalogs object. 
    Uses the query_object method from astroquery.
    
    Parameters
    ----------
    name: string
        The name of the object
        
    radius: float
        The radius within to search of that object.
        
    **kwargs
        Other keyword arguments that can narrow down the search. Documentation here: https://mast.stsci.edu/api/v0/_services.html
    
    Results
    -------
    obs: Observations or Catalogs
        Returns an astroquery Observations or Catalogs object of a specific object, depending on whether a catalog is specified.
    '''
    radius_string = str(radius) + " deg" # this is the format required for some reason
    
    obs = Observations.query_object(objectname=name,radius=radius_string,**kwargs) # from astroquery package
    
    obs_df = obs.to_pandas() # easier to work with pandas df then astropy table
        
    print("Number of results:",len(obs))
    return obs_df

# ---------------------------------------------------------------------------------------------------
# INSPECT THE OBSERVATION DATAFRAME
# ---------------------------------------------------------------------------------------------------

def obs_filters(observation):
    '''
    Determine what filters exist for a dataframe of observations for a specific astronomical object.
    
    Parameters
    ----------
    observation: DataFrame
        The DataFrame of observations.
        
    Returns
    -------
    numpy array
        Returns an array of the possible existing filters.
    '''
    return np.unique(np.asarray(observation.filters.values).astype(np.dtype(str)))

def obs_catalogs(observation):
    '''
    Determine what catalogs exist for a dataframe of observations for a specific astronomical object.
    
    Parameters
    ----------
    observation: DataFrame
        The DataFrame of observations.
        
    Returns
    -------
    numpy array
        Returns an array of the possible existing catalogs.
    '''
    return np.unique(np.asarray(observation.obs_collection))

def obs_table_filter_catalog(obs, filtername,catalog):
    '''
    Takes an astropy Table of observations and sorts them by filter name.
    
    Parameters
    ----------
    obs: astropy Table
        The list of observations that are being filtered for the filtername.
        
    filtername: string
        The name of the filter that the observations will be narrowed down by.
        
    catalog: string
        The name of the catalog that the observation will be drawn from.
    
    Returns
    -------
    astropy Table
        A list of observations at only one specific filter for one catalog.
    '''
    # check that both the filter and catalog exist simultaneously for a certain object
    if filtername not in obs_filters(obs):
        print('No observations were made of this object in that filter.')
        return None
    
    if catalog not in obs_catalogs(obs):
        print('No observations were made of this object in that catalog.')
        return None

    if catalog not in obs_catalogs(obs.loc[(obs.filters==filtername)]):
        print('No observations were made in this catalog of that filter.')
        return None
    
    else:
        return obs.loc[(obs.filters==filtername)&(obs.obs_collection==catalog)]

# ---------------------------------------------------------------------------------------------------
# DOWNLOAD FITS FILES
# ---------------------------------------------------------------------------------------------------
    
def download_products(obsID,download_dir=None,download_later = False):
    '''
    Download FITS files of a given observation onto the user's computer.
    
    Parameters
    ----------
    obsID: string
        The specific observation that will be downloaded.
        
    download_dir: string
        The directory to download the products to. Defaults to the current directory.
        
    download_later: Boolean
        Whether or not to download the FITS files now or later. Default is False and the images are downloaded when this function is run.
    
    Returns
    -------
    Downloaded FITS file onto the user's computer.
    '''
    if type(obsID) != str: # in case a user assumed it was a float
        obsID = 'obsID'
        
    return Observations.download_products(obsID,download_dir,curl_flag=download_later) # from astroquery package

def download_object_filter_catalog(name,radius=0.02,filtername='F606W',catalog='HST',obsIndex=0,download_dir = None,download_later=False):
    '''
    All in one! 
    Given that your filter and catalog are spelled correctly/exist, download an observation of an object in some filter from some catalog.
    Default is the first observation but can be changed.
    This function can be used alone but the other io functions can be used to gather more information about different available information.
    
    Parameters
    ----------
    name: string
        The name of the astronomical object being searched for.
    
    radius: float
        How wide a radius to look for observations within. If radius is too large the search will time out.
    
    filtername: string
        The specific filter to search for.
    
    catalog: string
        Which catalog to serach for.
    
    obsIndex: float
        Which observation to download (instead of downloading all of them). Default is the first observation.
        
    download_dir: string
        The directory to download the products to. Defaults to the current directory.
        
    download_later: Boolean
        Whether or not to download the FITS files now or later. Default is False and the images are downloaded when this function is run.
    
    Returns
    -------
    Downloaded FITS files of specified astronomical observations
    '''
    obs = query_object(name,radius)
    obs_filtered = obs_table_filter_catalog(obs, filtername, catalog)
    
    if obs_filtered is None:
        print('No observations were made in this catalog of that filter.')
        return None
    
    else:
        if len(obs_filtered[obs_filtered.t_exptime > 500]) > 0:
            obs_exptime = obs_filtered[obs_filtered.t_exptime > 500] # choosing longer exposures
        else:
            obs_exptime = obs_filtered # only if it exists

        if len(obs_exptime.loc[obs_exptime.target_name == name]) > 0: # choosing the observations that are explicitly of that object if they exist
            obs_final = obs_exptime.loc[obs_exptime.target_name == name]
        else:
            obs_final = obs_exptime

        return download_products(obs_filtered.obsid[obs_filtered.index[obsIndex]],download_dir,download_later)
