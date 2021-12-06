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

def plot(filename, objectname, pmax=97.5,z=10e-8,zoom=1):
    '''
    Plots an astronomical image both zoomed in and zoomed out. The user can set how much to zoom in and also what scale to be zoomed in.
    
    Parameters
    ----------
    filename: string
        The name of the FITS file to plot.
    
    objectname: string
        The name of the astronomical object.
    
    pmax: float
        The amount to scale the image by. Default is 97.5.
        
    z: float
        Redshift of a given astronomical object. If not given, default is 0.
    
    zoom: float
        The amount to zoom in on the image. Default is 1.
    '''
    image = fits.open(filename, ignore_missing_end=True)
    ra,dec,radius = resolve(objectname)
    hdr = image[0].header  
    try:
        try:
            filtername = hdr['FILTER']
        except:
            filtername = hdr['FILTER2']
    except:
        filtername = 'No filter name in FITS file'
        
    if type(image[0].data) != np.ndarray:
        hdr = image[1].header 
        if type(image[1].data) != np.ndarray:
            filtername = hdr['FILTER2']
    
    image.close()
    
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    asecperkpc=cosmo.arcsec_per_kpc_comoving(z).value

    fig = plt.figure(figsize=(15, 7))

    f1 = aplpy.FITSFigure(filename, hdu=1, north=True, figure=fig)
    if z == 10e-8:
        print("Scaled as if redshift were 0.")
    f1.recenter(ra,dec,radius=radius/zoom*asecperkpc/3600.)
    
    f1.set_title(objectname,size=24, weight='heavy',color='black', loc='left',family='serif')

    f1.add_label(0.05, 0.9, filtername, relative=True,size=20,
                 weight='heavy',color='red', horizontalalignment='left',family='serif')

    f1.show_colorscale(pmin=0,pmax=pmax,cmap='gray')#,stretch='arcsinh')
    
    theta = 5/3600. * u.arcsec                   # angle
    r_ang = Planck15.kpc_proper_per_arcmin(z) # phys. dist. per angle
    r     = r_ang * theta * 60                    # physical distance
    
    if z != 0:
        f1.add_scalebar(5./3600.)
        f1.scalebar.set_label(str('5" or ')+ str(round(r.value,3)*1000) +" pc")
        f1.scalebar.set_font(size=20,weight='heavy',family='serif')
        f1.scalebar.set_color('red')
    
    f1.add_colorbar()
    f1.colorbar.set_location('right')
    f1.colorbar.set_pad(0.1)
    f1.colorbar.set_axis_label_font(family='serif',size=15)
    pass

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
    # only run this once
    if method != 'divide' & 'add' &'subtract' &'multiply'&'mean':
        print('Method specified does not exist.')
    
    if os.path.isfile(fitsfilter2.split('.')[0]+method+'_combined_image.fits') is False:
        
        # the second image gets shifted
        pixelrange=20

        # this shifts the second file and saves a new image
        # in my research want to do IR/optical, shift_save returns optical image
        if os.path.isfile(fitsfilter2.split('.')[0]+'_shifted.fits') is False:
            shift_save(fitsfilter1,fitsfilter2,pixelrange,fitsfilter2.split('.')[0]+'_shifted.fits')

        print('Shifted the image')

        hdu1 = fits.open(fitsfilter1, ignore_missing_end=True)
        image1 = hdu1[0]
        if type(image1.data) != np.ndarray:
            image1 = hdu1[1]
        elif type(image1.data) != np.ndarray:
            print('Cannot open FITS file',fitsfilter1,'with extension 1 or 0')

        hdu2 = fits.open(fitsfilter2.split('.')[0]+'_shifted.fits',ignore_missing_end=True)
        image2 = hdu2[0]
        if type(image2.data) != np.ndarray:
            image2 = hdu2[1]
        elif type(image2.data) != np.ndarray:
            print('Cannot open FITS file',fitsfilter2.split('.')[0]+'_shifted.fits','with extension 1 or 0')

        #plt.imshow(numerator_image.data)
        #plt.imshow(image1.data)
        
        if method == 'add':
            with np.errstate(divide='ignore', invalid='ignore'):
                image_array = image1.data+image2.data
                
        if method == 'subtract':
            with np.errstate(divide='ignore', invalid='ignore'):
                image_array = image1.data-image2.data
                
        if method == 'multiply':
            with np.errstate(divide='ignore', invalid='ignore'):
                image_array = image1.data*image2.data
                
        if method == 'divide':
            with np.errstate(divide='ignore', invalid='ignore'):
                image_array = image1.data/image2.data
                
        if method == 'mean':
            with np.errstate(divide='ignore', invalid='ignore'):
                image_array = (image1.data+image2.data)/2

        new_hdu = fits.PrimaryHDU(image_array,header = hdu2[0].header)
        new_hdu.writeto(fitsfilter2.split('.')[0]+method+'_combined_image.fits', overwrite=True)

        combined = fits.open(fitsfilter2.split('.')[0]+method+'_combined_image.fits',mode='update', ignore_missing_end=True)

        image = open_fits(fitsfilter2)

        w = WCS(image.header)

        # change the header to the wcs from the other file
        combined[0].header.update(w.to_header())

        combined.close()

        return fitsfilter2.split('.')[0]+method+'_combined_image.fits'
    
def combined_plot(fitsfilter1,fitsfilter2, objectname, method = None, pmax=97.5,z=10e-8,zoom=1):
    '''
    Plots an astronomical image both zoomed in and zoomed out. The user can set how much to zoom in and also what scale to be zoomed in.
    
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
    asecperkpc=cosmo.arcsec_per_kpc_comoving(z).value
    
    ra,dec,radius = resolve(objectname)
    
    fileresult = combine(fitsfilter1,fitsfilter2,method)
    f1 = aplpy.FITSFigure(fileresult, hdu=0, north=True, figure=fig)
    if z == 10e-8:
        print("Scaled as if redshift were 0.")
        f1.recenter(ra,dec,radius=radius/zoom*asecperkpc/3600.)
        f1.set_title(objectname,size=24, weight='heavy',color='black', loc='left',family='serif')
        f1.add_label(0.05, 0.9, method, relative=True,size=20,
                     weight='heavy',color='red', horizontalalignment='left',family='serif')
        f1.show_colorscale(pmin=0,pmax=pmax,cmap='gray')#,stretch='arcsinh')

        theta = 5/3600. * u.arcsec                   # angle
        r_ang = Planck15.kpc_proper_per_arcmin(z) # phys. dist. per angle
        r     = r_ang * theta * 60                    # physical distance

    if z != 0:
        theta = 5/3600. * u.arcsec                   # angle
        r_ang = Planck15.kpc_proper_per_arcmin(z) # phys. dist. per angle
        r     = r_ang * theta * 60                    # physical distance

        f1.add_scalebar(5./3600.)
        f1.scalebar.set_label(str('5" or ')+ str(round(r.value,3)*1000) +" pc")
        f1.scalebar.set_font(size=20,weight='heavy',family='serif')
        f1.scalebar.set_color('red')

    f1.add_colorbar()
    f1.colorbar.set_location('right')
    f1.colorbar.set_pad(0.1)
    f1.colorbar.set_axis_label_font(family='serif',size=15)
        
    pass
