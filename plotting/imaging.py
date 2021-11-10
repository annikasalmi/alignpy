from astropy.io import fits 
import os
import numpy as np
import matplotlib.pyplot as plt

# for Lab 7
__all__ = ['IdealObservation', 'ArtObservation', 'IdealImager', 'ArtImager']

def implot(image,figsize=(15,13),cmap='gray_r',scale=0.5,colorbar=False,header=None,wcs=None,grid=False,**kwargs):
    ''' 
    A function which requires an image input that is a 2d array or masked array. The function has the optional inputs 
    of the float tuple figsize, the string cmap, and the float scale. A colorbar, header, and wcs can also be included.
    There are two **kwargs, vmin and/or vmax. 
    It plots an image and returns the fig, ax objects. 
    
    Parameters
    ----------
    image: 2-d ndarray, MaskedArray
        A 2d array or masked array of floats that is the image data.
    figsize: tuple
        The figure size of the image.
    cmap: string
        The color map scale. See options in the matplotlib color map documentation. Default is 'gray_r'.
    scale: float
        How the image color is scaled by vmin and vmax. Default is 0.5.
    colorbar: Boolean
        Whether or not a colorbar is included. Default is None.
    header:
        A header from a fits file. The default option is None.
    wcs:
        The world coordinate system option. The default is None.
    grid:
        Adds a grid if set to true.
    **kwargs
        Additional keyword arguments vmin and vmax, both floats, can be mentioned here.
    
    Returns
    -------
    fig, ax
        Returns the fig and ax objects from the subplot call.
    '''
    if wcs is not None:
        fig, ax = plt.subplots(1,figsize=figsize,subplot_kw={'projection':wcs})
        ax.set_xlabel("Right Ascension [hms]",fontsize=15)
        ax.set_ylabel("Declination [degrees]",fontsize=15)
        # if grid is True, add that
        if grid:
            ax.coords.grid(color='gray', alpha=0.5, linestyle='solid')

    if (header is not None) and (wcs is None):
        wcsobject = WCS(strip_SIP(header))
        fig, ax = plt.subplots(1,figsize=figsize,subplot_kw={'projection':wcsobject})
        ax.set_xlabel("Right Ascension [hms]",fontsize=15)
        ax.set_ylabel("Declination [degrees]",fontsize=15)
        # if grid is True, add that
        if grid:
            ax.coords.grid(color='gray', alpha=0.5, linestyle='solid')
        
    if (wcs is None) and (header is None):
        fig, ax = plt.subplots(1,figsize=figsize)
    
    vmin = kwargs.get("vmin",np.mean(image) - scale*np.std(image))
    vmax = kwargs.get("vmax",np.mean(image) + scale*np.std(image))

    ax.tick_params(direction="in", length = 30,labelsize=15)
    
    plot = ax.imshow(image,origin='lower',cmap=cmap,vmin=vmin,vmax=vmax)
    
    if colorbar:
        fig.colorbar(mappable = plot,ax=ax)
    
    return fig, ax 