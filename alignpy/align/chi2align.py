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

import drizzlepac
from drizzlepac import tweakreg

import json, requests

from reproject import reproject_interp

import os

import aplpy

# from project
from search import *
from plotting import *
from inputoutput import *

def counts_to_flux(fitspath):
    '''
    Turns image counts into flux values. Doesn't return anything; instead just modifies the file
    
    Parameters
    ----------
    fitspath: FITS file
        The fits file that needs to be converted to flux values.
    '''
    # ADD HAS ATTRIBUTE WHEN I TURN INTO A CLASS SO IT ONLY RUNS ONCE
    image = open_fits(fitspath,mode='update')
    try: 
        photflam=image.header['PHOTFLAM']
    except: 
        photflam=2.901E-19
    try: 
        photplam=image.header['PHOTPLAM']
    except: 
        photplam=16030.43
    data = image.data
    zpab = -2.5*np.log10(photflam)-5*np.log10(photplam)-2.408
    magab = -2.5*np.log10(data)+zpab
    data = magab * 606 * 10**9 / (2.9979 * 10**8)
    image.data = data
    
    warnings.filterwarnings("ignore")    
    pass

def reproject(fitsfilter1,fitsfilter2,extfilter1,extfilter2,reproject_name):
    '''
    Reproject two FITS images so that the smaller image is on the same pixel scale as the larger image.
    Can also reproject in the other direction but that is a loss of information.f
    
    Parameters
    ----------
    fitsfilter1: string
        The location of a specific FITS file in one filter. This is the larger image.
        
    fitsfilter2: string
        The location of a specific FITS file in a different filter.
        
    extfilter1: string
        The extension of the first FITS file.
        
    extfilter2: string
        The extensino of the second FITS file.
        
    reproject_name: string
        The name of the reprojected FITS file.
    
    Returns
    -------
    reproject_path: string
        Returns the path location of the reprojected FITS file.
    '''
    filter1hdu = fits.open(fitsfilter1, ignore_missing_end=True)
    filter1image=filter1hdu[extfilter1] # can't use my fxn because i need the 'sci' and 'wht' keywords
    filter2hdu = fits.open(fitsfilter2, ignore_missing_end=True)
    filter2image=filter2hdu[extfilter2]
    
    if '.' in reproject_name:
        reproject_path = reproject_name
        
    else:
        reproject_path = reproject_name + '.fits'
    
    if os.path.isfile(reproject_path) is True:
        print('Reprojected file already exists')
        return reproject_path
    
    if os.path.isfile(reproject_path) is False:
        # reprojects filter2 onto size of filter1
        array, footprint = reproject_interp(filter2image, filter1image.header)

        fits.writeto(reproject_path, array, filter1image.header, overwrite=True)

        filter1hdu.close()
        filter2hdu.close()
    
    return reproject_path

def reproject_science(fitsfilter1,fitsfilter2,reproject_name):
    '''
    Reproject two FITS images so that the smaller image is on the same pixel scale as the larger image.
    Can also reproject in the other direction but that is a loss of information.
    Reprojects only for science data.
    
    Parameters
    ----------
    fitsfilter1: string
        The location of a specific FITS file in one filter. This is the larger image.
        
    fitsfilter2: string
        The location of a specific FITS file in a different filter.
        
    reproject_name: string
        The name of the reprojected FITS file.
        
    Returns
    -------
    reproject_path: string
        Returns the path location of the reprojected FITS file.
    '''
    
    reproject_path = reproject(fitsfilter1,fitsfilter2,'SCI','SCI',reproject_name)
    return reproject_path

def reproject_error(fitsfilter1,fitsfilter2,reproject_name):
    '''
    Reproject two FITS images so that the smaller image is on the same pixel scale as the larger image.
    Can also reproject in the other direction but that is a loss of information.
    Reprojects only for error data.
    
    Parameters
    ----------
    fitsfilter1: string
        The location of a specific FITS file in one filter. This is the larger image.
        
    fitsfilter2: string
        The location of a specific FITS file in a different filter.
        
    reproject_name: string
        The name of the reprojected FITS file.
        
    Returns
    -------
    reproject_path: string
        Returns the path location of the reprojected FITS file.
    '''
    if 'WHT' in np.asarray(fits.open(fitsfilter2).info(), ignore_missing_end=True):
        # adjust weights to errors
        reproject_err_path = reproject(fitsfilter1,fitsfilter2,'WHT','SCI',reproject_name)
        image = fits.open(reproject_name, ignore_missing_end=True)[0]
        data = 1/image.data
        image.data = data
        return reproject_err_path
        
    if 'ERR' in np.asarray(fits.open(fitsfilter1, ignore_missing_end=True).info()):
        reproject_err_path = reproject(fitsfilter1,fitsfilter2,'ERR','SCI',reproject_name)
        return reproject_err_path

    else:
        print('No error file present')
        return None
    
def shift_images(fitsfilter1,fitsfilter2,pixelrange):
    '''
    Aligns a given FITS file with another fits file by subtracting the two images and then minimizing the added flux. 
    Checks along a range of different pixel values for the total chi2 and then 
    
    Parameters
    ----------
    fitsfilter1: string
        The location of the FITS file in a specific filter that will be shifted.
        
    fitsfilter2: string
        The location of the FITS file in a different filter that will be compared against.
        
    pixelrange: integer
        The range of pixels to adjust the search over. Ranges from negative of that value to positive of that value.
        
    Returns
    -------
    image_shift: ndarray
        The first FITS file as an array shifted by the amount necessary to minimize the subtraction with the second FITS file.
        
    x_min: float
        The amount to shift the image by in the x direction.
        
    y_min: float
        The amount to shift the image by in the y direction.
    '''
    # reprojects the second image "fits_comparison"; this should be the smalle
    reproject_path = reproject_science(fitsfilter1,fitsfilter2,fitsfilter2.split('.')[0]+'_reproject.fits')
    reproject_err_path = reproject_error(fitsfilter1,fitsfilter2,fitsfilter2.split('.')[0]+'_reproject_err.fits')
    
    # figure out the extensions
    # can't use my function since can't close files till later
    hdu_comparison = fits.open(fitsfilter1, ignore_missing_end=True)
    image_comparison = hdu_comparison[0]
    if type(image_comparison.data) != np.ndarray:
        image_comparison = hdu_comparison[1]
        
    hdu_shift = fits.open(reproject_path, ignore_missing_end=True)
    image_shift = hdu_shift[0]
    if type(image_shift.data) != np.ndarray:
        image_shift = hdu_shift[1]
    
    # if error path exists, check extensions 
    if reproject_err_path is not None:
        hdu_err = fits.open(reproject_err_path, ignore_missing_end=True)
        err = hdu_err[0]
        if type(err.data) != np.ndarray:
            err = hdu_err[1]      
    else: # if no error path exists just set err to None
        err = None
    
    pixel_grid  = np.arange(-pixelrange,pixelrange,1)

    x_arr,y_arr,chi2 = [],[],[]
    chi2_min = 1e20
    x_min, y_min = 0,0

    image = np.zeros(shape=(len(pixel_grid),len(pixel_grid)))
    chi2 = []
    print('Shifting ',fitsfilter1,'to be aligned with ',fitsfilter2)
    for i,x in enumerate(pixel_grid):
        for j,y in enumerate(pixel_grid):
            
            
            # MOVE THE IMAGE OVER
            image_xshift = np.roll(image_shift.data,pixel_grid[i],axis=1) # shifting image in x range
            image_xshift[:,pixel_grid[i]] = np.mean(image_shift.data) # fills in edges with the mean of the data instead
            image_xyshift = np.roll(image_xshift,pixel_grid[j],axis=0)
            image_xyshift[:,pixel_grid[j]] = np.mean(image_xshift.data)
            
            # CALCULATE CHI2
            subtraction = np.abs(image_comparison.data - image_xyshift) # need absolute value otherwise adding it all up doesn't make sense as negatives cancel each other
            flux_sum = np.sum(subtraction)

            if err is not None:
                c = (flux_sum-err)**2
                
            else:
                c = flux_sum**2

            chi2.append(c)
            x_arr.append(x)
            y_arr.append(y)
            image[i,j] = c
            if c < chi2_min:
                chi2_min = c
                x_min    = x
                y_min    = y

    chi2 = np.array(chi2)
    x_arr = np.array(x_arr)
    y_arr = np.array(y_arr)
    
    final_image_xshift = np.roll(image_shift.data,x_min,axis=1) # fills in edges with zeroes; shifting image in both x and y range
    final_image_xyshift = np.roll(final_image_xshift,y_min,axis=0)
    
    hdu_comparison.close()
    hdu_shift.close()
    if reproject_err_path is not None:
        hdu_err.close()
    
    return final_image_xyshift,x_min,y_min

def shift_save(fitsfilter1,fitsfilter2,pixelrange,outputfilename):
    '''
    Shifts a fits file so that when subtracted from another file, its chi2 is minimized. Saves the shifted array. 
    If update is True, it updates the old file to the new data.
    
    Parameters
    ----------
    fitsfilter1: string
        The location of the FITS file that will be shifted.
        
    fitsfilter2: string
        The location of the FITS file that will be compared against.
        
    pixelrange: integer
        The range of pixels to adjust the search over. Ranges from negative of that value to positive of that value.
        
    outputfilename: string
        What to name the new FITS file.
    '''
    if '.' in outputfilename:
        outputfilename = outputfilename
        
    else:
        outputfilename = outputfilename + '.fits'
        
    if os.path.isfile(outputfilename) is True:
        print('Shifted file already exists')
        return outputfilename
    
    if os.path.isfile(outputfilename) is False:
        output_arr,x_min,y_min = shift_images(fitsfilter1,fitsfilter2,pixelrange)
        hdu = fits.PrimaryHDU(data=output_arr)
        hdu.writeto(outputfilename)

        output = fits.open(outputfilename,mode='update', ignore_missing_end=True)

        image_shift = fits.open(fitsfilter1, ignore_missing_end=True)[0]
        if type(image_shift.data) != np.ndarray:
            image_shift = fits.open(fitsfilter1, ignore_missing_end=True)[1]

        hdu = fits.open(fitsfilter1, ignore_missing_end=True)
        image = hdu[0]
        if type(image.data) != np.ndarray:
            image = hdu[1]
        elif type(image.data) != np.ndarray:
            print('Cannot open FITS file',fitsfilter1,'with extension 1 or 0')
            
        w = WCS(image.header)

        # change the header to the wcs from the other file
        output[0].header.update(w.to_header())

        output.close()
        return outputfilename

