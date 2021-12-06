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

def gaia_align(fitspath,objname,radius = 6.0,minobj=2,sigma=5.0,threshold=5., conv_width=6):
    '''
    Aligns a given FITS file to the GAIA catalog. Doesn't return anything; just modifies the file.
    To be used only if a FITS file doesn't appear to have WCS or the WCS seems very off.
    
    Parameters
    ----------
    fitspath: string
        The location of a specific FITS file.
        
    objname: string
        An astronomical object we are aligning an image towards in the GAIA catalog.
        
    radius: string
        The radius in arcseconds that we will search around the given object.
        
    minobj: integer
        The minimum number of objects that need to be found in order for the FITS file to be aligned with the GAIA catalog.
    
    sigma: float
        The allowed error for the alignment process.
        
    threshold: float
        The minimum threshold for finding matchign objects.
        
    conv_width: int
        The minimum convergence width for finding objects.
    '''
    coord = resolve(objname)
    radius = Quantity(radius, u.arcmin)
    coord = SkyCoord(ra=coord[0], dec=coord[1], unit=(u.deg, u.deg))
    gaia_query = Gaia.query_object_async(coordinate=coord, radius=radius)
    
    reduced_query = gaia_query['ra', 'dec', 'phot_g_mean_mag']
    reduced_query.write('gaia.cat', format='ascii.commented_header')
    

    cw = 3.5  # Set to two times the FWHM of the PSF.
    wcsname = 'Gaia'  # Specify the WCS name for this alignment

    # ALIGNING TO GAIA
    if os.path.isdir('gaia.cat'): #rewrite gaia.cat if necessary
        refcat = None
    refcat = 'gaia.cat'
    
    # I got TweakReg from STSci but I wrote the rest of this function
    tweakreg.TweakReg(fitspath,  # Pass input images
                      updatehdr=True,  # update header with new WCS solution
                      imagefindcfg={'threshold':threshold,'conv_width':conv_width},  # Detection parameters
                                                                     # threshold varies for different data
                      refcat=refcat,  # Use user supplied catalog (Gaia)
                      interactive=False,
                      see2dplot=False,
                      minobj = minobj,
                      shiftfile=True,  # Save out shift file (so we can look at shifts later)      
                      wcsname=wcsname,  # Give our WCS a new name
                      reusename=True,
                      sigma=sigma,
                      ylimit=0.2,
                      fitgeometry='general')  # Use the 6 parameter fit  
    pass
    
def gaia_filealign(fitsfilter1, fitsfilter2, minobj=2):
    '''
    Aligns a given FITS file with another fits file. Doesn't return anything; just modifies the file.
    Preferably the chi2_images function should be used over this one.
    
    Parameters
    ----------
    fitsfilter1: string
        The location of the FITS file that will be compared against.
        
    fitsfilter2: string
        The location of the FITS file that will be shifted.
        
    minobj: integer
        The minimum number of objects that need to be found in order for the FITS file to be aligned with the GAIA catalog.
    '''
    tweakreg.TweakReg(fitsfilter2,
                      enforce_user_order=False,
                      imagefindcfg={'threshold': 10, 'conv_width': 3.5, 'dqbits': ~4096},
                      minobj = minobj,
                      refimage=fitsfilter1, 
                      refimagefindcfg={'threshold': 10, 'conv_width': 2.5},
                      shiftfile=True,
                      outshifts='shift657_flc.txt',
                      searchrad=5.0,
                      ylimit=0.6,
                      updatehdr=True,
                      #updatewcs=True,
                      wcsname='UVIS_FLC',
                      reusename=True,
                      interactive=False)
    pass

def query_object(name, radius, **kwargs):
    '''
    Queries an object within a specified catalog and creates an astroquery Catalogs object. Uses the query_object method from astroquery.
    
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
    
    obs = Observations.query_object(objectname=name,radius=radius_string,**kwargs)
    
    obs_df = obs.to_pandas()
        
    # else:
    #     obs = Observations.query_object(objectname=name,radius=radius_string,**kwargs)
    #     obs = obs.loc[(obs['obs_collection'] == catalog)]
        
    print("Number of results:",len(obs))
    return obs_df


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
    #return np.unique(np.asarray(observation['filters']))
    return np.unique(np.asarray(observation.filters.values).astype(np.dtype(str)))


# Determine what catalogs exist for a specific observation
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
    if filtername not in obs_filters(obs):
        return 'No observations were made of this object in that filter.'
        
    if catalog not in obs_catalogs(obs):
        return 'No observations were made of this object in that catalog.'

    if catalog not in obs_catalogs(obs.loc[(obs.filters==filtername)]):
        return 'No observations were made in this catalog of that filter.'
    
    else:
        return obs.loc[(obs.filters==filtername)&(obs.obs_collection==catalog)]

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
        
    return Observations.download_products(obsID,download_dir,curl_flag=download_later)

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
    return download_products(obs_filtered.obsid[obs_filtered.index[obsIndex]],download_dir,download_later)
