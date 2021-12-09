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

# ---------------------------------------------------------------------------------------------------
# ALIGN TO EXTERNAL CATALOGS 
# ---------------------------------------------------------------------------------------------------

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
    # TweakReg from STSci but I wrote the rest of this function
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

# ---------------------------------------------------------------------------------------------------
# CHANGE FILE VALUES FROM FLUX TO COUNT 
# ---------------------------------------------------------------------------------------------------

def counts_to_flux(fitspath):
    '''
    Turns image counts into flux values. Doesn't return anything; instead just modifies the file
    
    Parameters
    ----------
    fitspath: FITS file
        The fits file that needs to be converted to flux values.
    '''
    image, hdu = search.open_fits(fitspath,mode='update')
    
    try: # make sure I only change the counts to flux once for each file
        # I'm not sure how to do this without try/except, since header['PHOTFLAM'] is not a local var, global var, or an attribute
        if image.header.comments['PHOTFLAM'] == 'Counts changed to flux':
            pass
    except:
        try: 
            photflam=image.header['PHOTFLAM'] # if the photflam value is present, use that; otherwise use a set value
            image.header['PHOTFLAM'] = image.header['PHOTFLAM'],'Counts changed to flux'
        except: 
            photflam=2.901E-19
            image.header.append('PHOTFLAM',2.901E-19,'Counts changed to flux')
        try: 
            photplam=image.header['PHOTPLAM'] # if the photplam value is present, use that; otherwise use a set value
        except: 
            photplam=16030.43
        data = image.data
        zpab = -2.5*np.log10(photflam)-5*np.log10(photplam)-2.408
        magab = -2.5*np.log10(data)+zpab
        data = magab * 606 * 10**9 / (2.9979 * 10**8)
        image.data = data # change data values of the file 
        
        hdu.close()

        pass

# ---------------------------------------------------------------------------------------------------
# REPROJECT FILES TO MATCH EACH OTHER
# ---------------------------------------------------------------------------------------------------
    
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
    filter1image, filter1hdu = search.open_fits(fitsfilter1,ext=extfilter1) 
    filter2image, filter2hdu = search.open_fits(fitsfilter2,ext=extfilter2)
    
    # check that the user properly named the reproject file to include .fits
    if '.fits' in reproject_name:
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
    if 'WHT' in np.asarray(fits.open(fitsfilter2).info()): # can't use open_fits because need the info extension
        # adjust weights to errors
        reproject_err_path = reproject(fitsfilter1,fitsfilter2,'WHT','SCI',reproject_name)
        image,hdu = search.open_fits(reproject_name)
        data = 1/image.data
        image.data = data
        hdu.close
        return reproject_err_path
        
    if 'ERR' in np.asarray(fits.open(fitsfilter1).info()):
        reproject_err_path = reproject(fitsfilter1,fitsfilter2,'ERR','SCI',reproject_name)
        return reproject_err_path

    else:
        print('No error file present')
        return None
    
# ---------------------------------------------------------------------------------------------------
# SHIFT FILES TO MATCH EACH OTHER
# ---------------------------------------------------------------------------------------------------
    
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
    image_comparison, hdu_comparison = search.open_fits(fitsfilter1)
    
    image_shift, hdu_shift = search.open_fits(reproject_path)
    
    # if error path exists, check extensions 
    if reproject_err_path is not None:
        image_err, hdu_err = search.open_fits(reproject_err_path)
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
            subtraction = np.abs(image_comparison.data - image_xyshift) # need absolute value 
                                                                        # otherwise adding it all up doesn't make sense as negatives cancel each other
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
    if '.fits' in outputfilename:
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

        output = fits.open(outputfilename,mode='update')

        image, hdu = search.open_fits(fitsfilter1)
            
        w = WCS(image.header)

        # change the header to the wcs from the other file
        output[0].header.update(w.to_header())

        output.close()
        hdu.close()
        return outputfilename
