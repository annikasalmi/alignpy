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

from alignpy import *
from io import *
from align import *

# these are functions that aligning and downloading requires

def open_fits(fitspath,mode='readonly'):
    '''
    Opens a fits file. Determines which extension to use, assuming it is either 1 or 0.
    
    Parameters
    ----------
    fitspath: string
        The path of the FITS file to be opened.
        
    mode: string
        Whether readonly or update. Default is readonly.
        
    Returns
    -------
    image: FITS object
    '''
    hdu = fits.open(fitspath,mode, ignore_missing_end=True)
    image = hdu[0]
    if type(image.data) != np.ndarray:
        image = hdu[1]
    elif type(image.data) != np.ndarray:
        print('Cannot open FITS file',fitspath,'with extension 1 or 0')
        
    hdu.close()
    return image


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

def resolve(name):
    '''
    Get the RA and Dec for an object using the MAST name resolver.
    
    Parameters
    ----------
    name: string
        Name of object

    Returns
    -------
    objRa: float
        RA position
    objDec: float
        Dec position
    '''

    resolverRequest = {'service':'Mast.Name.Lookup',
                       'params':{'input':name,
                                 'format':'json'
                                },
                      }
    resolvedObjectString = mastQuery(resolverRequest)
    resolvedObject = json.loads(resolvedObjectString)
    # The resolver returns a variety of information about the resolved object, 
    # however for our purposes all we need are the RA and Dec
    try:
        objRa = resolvedObject['resolvedCoordinate'][0]['ra']
        objDec = resolvedObject['resolvedCoordinate'][0]['decl']
        objRadius = resolvedObject['resolvedCoordinate'][0]['radius']
    except IndexError as e:
        raise ValueError("Unknown object '{}'".format(name))
    return objRa, objDec, objRadius
