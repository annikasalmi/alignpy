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

# internal imports
from search import *
from align import *
from inputoutput import *

# ---------------------------------------------------------------------------------------------------
# COMBINE TWO FITS FILES TOGETHER
# ---------------------------------------------------------------------------------------------------

def combine(fitsfilter1,fitsfilter2,method):
    '''
    Combines two files together other.
    
    Parameters
    ----------
    fitsfilter1: string
        The location of the file to use as one file to multiply.
    
    fitsfilter2: string
        The location of the file to use as a second file to multiply.  
    '''
    chi2align.counts_to_flux(fitsfilter1) # change this here so that the combined flux counts are correct
    chi2align.counts_to_flux(fitsfilter2)
    # only run this once: check that the method is correct
    while method not in ['add', 'subtract', 'multiply','divide','mean']:
        return 'Method specified does not exist.'
    
    # check that the user hasn't already created a combined image
    if os.path.isfile(fitsfilter2.split('.')[0]+method+'_combined_image.fits') is False:
        
        # the second image gets shifted
        pixelrange=30

        # this shifts the second file and saves a new image
        # in my research want to do infrared/optical, shift_save returns optical image
        if os.path.isfile(fitsfilter2.split('.')[0]+'_shifted.fits') is False:
            chi2align.shift_save(fitsfilter1,fitsfilter2,pixelrange,fitsfilter2.split('.')[0]+'_shifted.fits')

        print('Shifted the image')

        image1, hdu1 = search.open_fits(fitsfilter1)

        image2, hdu2 = search.open_fits(fitsfilter2.split('.')[0]+'_shifted.fits')
        
        with np.errstate(divide='ignore', invalid='ignore'): # if two numbers won't go together, ignore it
            if method == 'add':
                image_array = image1.data+image2.data
                
            if method == 'subtract':
                image_array = image1.data-image2.data
                
            if method == 'multiply':
                image_array = image1.data*image2.data
                
            if method == 'divide':
                image_array = image1.data/image2.data
                image_array = np.nan_to_num(image_array)
                
            if method == 'mean':
                image_array = (image1.data+image2.data)/2
                
        # SAVE A NEW FITS FILE
        new_hdu = fits.PrimaryHDU(image_array,header = hdu2[0].header)
        new_hdu.writeto(fitsfilter2.split('.')[0]+method+'_combined_image.fits', overwrite=True)

        combined = fits.open(fitsfilter2.split('.')[0]+method+'_combined_image.fits',mode='update', ignore_missing_end=True)

        w = WCS(image1.header)

        # change the header to the wcs from the other file
        combined[0].header.update(w.to_header())

        combined.close()
        hdu1.close()
        hdu2.close()
        
        # Return the name of the new file
        return fitsfilter2.split('.')[0]+method+'_combined_image.fits'
    
    else:
        return fitsfilter2.split('.')[0]+method+'_combined_image.fits'

# ---------------------------------------------------------------------------------------------------
# PLOTTING OPTIONS
# ---------------------------------------------------------------------------------------------------

def plot(filename, objectname, pmax=97.5,z=10e-8,zoom=1,radiussearch=True): 
    '''
    Plots an astronomical image both zoomed in and zoomed out. 
    The user can set how much to zoom in and also what scale to be zoomed in.
    
    Parameters
    ----------
    filename: string
        The name of the FITS file to plot.
    
    objectname: string
        The name of the astronomical object.
    
    pmax: float
        The amount to scale the image by. Default is 97.5.
        
    z: float
        Redshift of a given astronomical object. If not given, default is 10e-8 (essentially 0 but avoids divide by zero error).
    
    zoom: float
        The amount to zoom in on the image. Default is 1.
    
    radius: Boolean or float
        If True, search for radius. If False, set radius to a specific amount.
    '''
    chi2align.counts_to_flux(filename)  # this only runs if counts haven't been changed to flux already
                                        # change this here as oppposed to elsewhere so the colorbar value is in fluxes
    hdu = fits.open(filename) # can't use search.open_fits because ['FILTER'] header could be not with the image array
    if radiussearch is True:
        ra,dec,radius = resolve(objectname)
    if radiussearch is not True:
        ra,dec = resolve(objectname,radius=False)
        radius = radiussearch
    
    hdr = hdu[0].header
   
    # need to do try/excpet because hdr['FILTER'] is not a variable and therefore can't check if it exists
    # returns a warning otherwise
    try:
        try:
            try: # checking that the header exists
                filtername = str(hdr['FILTER'])

            except:
                filtername = str(hdr['FILTER2'])

        except: 
            hdr = hdu[1].header
            try:
                filtername = str(hdr['FILTER'])

            except:
                filtername = str(hdr['FILTER2'])
    except:
            filtername = 'No filter in FITS file'
            
    hdu.close()
    
    # find arcseconds per kpc
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    asecperkpc=cosmo.arcsec_per_kpc_comoving(z).value

    fig = plt.figure(figsize=(15, 7))

    f1 = aplpy.FITSFigure(filename, hdu=1, north=True, figure=fig)
    f1.recenter(ra,dec,radius=radius/zoom*asecperkpc/3600.) # make sure it's centered at the correct radius and declination
    f1.set_title(objectname,size=24, weight='heavy',color='black', loc='left',family='serif')
    f1.add_label(0.05, 0.9, 'Filter: '+filtername, relative=True,size=20,
                 weight='heavy',color='red', horizontalalignment='left',family='serif')
    f1.show_colorscale(pmin=0,pmax=pmax,cmap='gray')#,stretch='arcsinh')
    
    if z == 10e-8:
        print("No z entered and cannot show kpc.")
    
    # add a scalebar if z exists, since now we can calculate the kpc across of the galaxy
    if z != 10e-8:
        theta = 5/3600. * u.arcsec                   # angle
        r_ang = Planck15.kpc_proper_per_arcmin(z) # phys. dist. per angle
        r = r_ang * theta * 60                    # physical distance
        f1.add_scalebar(5./3600.)
        f1.scalebar.set_label(str('5" or ')+ str(round(r.value,3)*1000) +" pc")
        f1.scalebar.set_font(size=20,weight='heavy',family='serif')
        f1.scalebar.set_color('red')
    
    # colorbar attributes
    f1.add_colorbar()
    f1.colorbar.set_location('right')
    f1.colorbar.set_pad(0.1)
    f1.colorbar.set_axis_label_text('Flux')
    f1.colorbar.set_axis_label_font(family='serif',size=15)
    pass
    
def combined_plot(fitsfilter1,fitsfilter2, objectname, method = None, pmax=97.5,z=10e-8,zoom=1,radiussearch=True):
    '''
    Plots an astronomical image both zoomed in and zoomed out. 
    The user can set how much to zoom in and also what scale to be zoomed in.
    
    Parameters
    ----------
    fitsfilter1: string
        The name of the FITS file to plot. In an attenuation map, this is the IR image.
        
    fitsfilter2: string
        The second name of the FITS file to plot. In an attenuation map, this is the optical image.
    
    objectname: string
        The name of the astronomical object.
        
    method: string
        The way to plot the two objects. Default is None.
    
    pmax: float
        The amount to scale the image by. Default is 97.5.
        
    z: float
        Redshift of a given astronomical object. If not given, default is 0.
    
    zoom: float
        The amount to zoom in on the image. Default is 1.
    '''
    fig = plt.figure(figsize=(15, 7))
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    asecperkpc=cosmo.arcsec_per_kpc_comoving(z).value # find arcseconds per kpc
    
    # determine size of the object as well as its location
    if radiussearch is True:
        ra,dec,radius = resolve(objectname)
    if radiussearch is not True:
        ra,dec = resolve(objectname,radius=radiussearch)
        radius = radiussearch
    
    fileresult = combine(fitsfilter1,fitsfilter2,method)
    f1 = aplpy.FITSFigure(fileresult, hdu=0, north=True, figure=fig)
    f1.recenter(ra,dec,radius=radius/zoom*asecperkpc/3600.)
    f1.set_title(objectname,size=24, weight='heavy',color='black', loc='left',family='serif')
    f1.add_label(0.05, 0.9, method, relative=True,size=20,
                 weight='heavy',color='red', horizontalalignment='left',family='serif')
    f1.show_colorscale(pmin=0,pmax=pmax,cmap='hot')#,stretch='arcsinh')

    if z == 10e-8:
        print("No z entered and cannot show kpc.")

    if z != 10e-8:
        theta = 5/3600. * u.arcsec                   # angle
        r_ang = Planck15.kpc_proper_per_arcmin(z) # phys. dist. per angle
        r     = r_ang * theta * 60                    # physical distance

        f1.add_scalebar(5./3600.)
        f1.scalebar.set_label(str('5" or ')+ str(round(r.value,3)*1000) +" pc")
        f1.scalebar.set_font(size=20,weight='heavy',family='serif')
        f1.scalebar.set_color('red')

    # colorbar settings
    f1.add_colorbar()
    f1.colorbar.set_location('right')
    f1.colorbar.set_pad(0.1)
    f1.colorbar.set_axis_label_text('Flux')
    f1.colorbar.set_axis_label_font(family='serif',size=15)
        
    pass
