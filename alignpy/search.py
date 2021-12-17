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

# these are functions that all three of the modules require
# ---------------------------------------------------------------------------------------------------
# OPEN FITS FILES
# ---------------------------------------------------------------------------------------------------

def open_fits(fitspath,mode='readonly',ext=None):
    '''
    Opens a fits file. Determines which extension to use, assuming it is either 1 or 0.
    
    Parameters
    ----------
    fitspath: string
        The path of the FITS file to be opened.
        
    mode: string
        Whether readonly or update. Default is readonly.
        
    ext: string
        Optionally, open a FITS file with a specific string value instead. Default is None
        
    Returns
    -------
    image: FITS object
    '''
    if ext is None: # check 0 and 1st extensions
        hdu = fits.open(fitspath,mode, ignore_missing_end=True)
        image = hdu[0] # assume science is in the 0th extension 
        if type(image.data) != np.ndarray:
            image = hdu[1] # if not, try 1st; sometimes it's there instead
        elif type(image.data) != np.ndarray:
            print('Cannot open FITS file',fitspath,'with extension 1 or 0')
    if ext is not None: # if extension specified, open the file there
        hdu = fits.open(fitspath,mode,ignore_missing_end=True)
        image = hdu[ext]
        if type(image.data) != np.ndarray:
            print('Cannot open FITS file',fitspath,'with extension',ext)

    # need to remember to close the hdu each time after
    return image, hdu

# ---------------------------------------------------------------------------------------------------
# DETERMINE PROPERTIES OF A GIVEN ASTRONOMICAL OBJECT
# ---------------------------------------------------------------------------------------------------

def mastQuery(request, url='https://mast.stsci.edu/api/v0/invoke'):
    '''
    Perform a MAST query for an astronomical object.

    Parameters
    ----------
    request: dictionary 
        The MAST request json object
    url: string 
        The service URL. Default is MAST

    Returns
    -------
    r.text: string
        Returns the coordinates of requested astronomical object.
    '''
    
    # Encoding the request as a json string
    requestString = json.dumps(request)
    r = requests.post(url, data={'request': requestString})
    r.raise_for_status()
    return r.text

def resolve(name,radius):
    '''
    Get the RA and Dec for an object using the MAST name resolver.
    
    Parameters
    ----------
    name: string
        Name of object
    radius: Boolean
        Whether or not to return the radius of an object

    Returns
    -------
    objRa: float
        RA position
    objDec: float
        Dec position
    objRadius: float
        radius position
    '''

    resolverRequest = {'service':'Mast.Name.Lookup',
                       'params':{'input':name,
                                 'format':'json'
                                },
                      }
    resolvedObjectString = mastQuery(resolverRequest)
    resolvedObject = json.loads(resolvedObjectString)
    # The resolver returns a variety of information about the resolved object, 
    # however for our purposes all we need are the RA, Dec, and radius
    if radius is not True:
        try:
            objRa = resolvedObject['resolvedCoordinate'][0]['ra']
            objDec = resolvedObject['resolvedCoordinate'][0]['decl']
        except IndexError as e:
            raise ValueError("Unknown object '{}'".format(name))
        return objRa, objDec
    if radius == True:
        try:
            objRa = resolvedObject['resolvedCoordinate'][0]['ra']
            objDec = resolvedObject['resolvedCoordinate'][0]['decl']
            objRadius = resolvedObject['resolvedCoordinate'][0]['radius']
        except IndexError as e:
            raise ValueError("Unknown object '{}'".format(name))
        return objRa, objDec, objRadius
