# Project Title

Welcome to **alignpy**! This project is intended to download FITS files of specific astronomical objects, specified by filter and catalog, onto your computer. Once they are downloaded, the files can be aligned together and plotted.

## Motivation

It is possible to align two different images with star patterns to each other by using astrometry.net and knowing the WCS. The reason alignpy is needed is that sometimes even if the images are nominally aligned, the stars are still off by a few pixels if the images are somehow combined together. This package also can help a user find what observations have been made of a specific object.

## Installation

**alignpy** can be installed by cloning the repo and then running
‘’’
pip install -e .
‘’’

**alignpy** is dependent on **numpy**, **scipy**, **matplotlib**, and **aplpy** to plot images. In order to download images, the packages **astropy**, **astroquery**, **json** and **requests** are necessary so that the MAST catalog can be searched via a script. It also depends on **os** to check the location of specific files. Finally, the astropy dependent packages **reproject** and **drizzlepac** are necessary so that FITS files can be aligned to a greater catalog and one FITS files can be reprojected to another.

Since these are a lot of dependencies, I recommend starting an environment so as to not break the user’s base python environment. To start an environment, run the following command in your terminal:


> conda create -n alignenv python=3.8, astropy = 4.2

**astropy** version 4.2 is required in order to be able to open FITS files that don’t have the header keyword [‘END’] in them. To activate the environment in the terminal, run:

> conda activate alignenv

To install the required packages within the environment, run:

> conda install -n alignenv numpy scipy matplotlib aplpy astropy astroquery json requests reproject drizzlepac os

Hit “y” when prompted.

## Usage
This package can be used to determine what observations have been made of an astronomical object, download specific FITS files, align two FITS files together, and plot FITS files.

```python
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

# Returns all observations made of a specific astronomical object within a specified radius.
obsM81 = query_object(‘M81’,0.02)

# Returns all the catalogs that contain observations of that astronomical object.

obs_catalogs(obsM81)

# Returns all the filters that the astronomical object has been observed in.

obs_filters(obsM81)

# Downloads a FITS file onto the user’s computer. The user searches for a specific astronomical object, and then 
# can optionally add a radius, filter, and catalog to download from. Defaults to optical images from HST with a 
# 0.02 radius around the object. 

alignpy.io.download_object_filter_catalog(‘M81’,0.02,’F606W’,’HST’)

# plots a FITS file. Write in the object name for a label and also to find its RA and Dec online. User can 
# manually input the z; otherwise it will default to close by. The image can be presented as is or zoomed 
# in further. pmax  adjusts the scale of the final plotted image. Returns all the downloaded files, which 
# include txt and other non-FITS files.
alignpy.plot.plot(filename, objectname, pmax=97.5,z=10e-8,zoom=1)

# Aligns two FITS files together and saves a new file. It only shifts the objects together over a specified 
# amount of space.
alignpy.align.shift_save(filename1,filename2,pixelrange,outputfilename)

# plots two FITS files together, either by adding, subtracting, dividing, or multiplying the two files. Aligns  
# the two files automatically.
alignpy.plot.combined_plot(filename1, filename2, objectname, method,pmax=97.5,z=10e-8,zoom=1)
```

## Documentation
Hosted on ReadTheDocs.

## Credits
Thank you to Professor Marla Geha for suggesting the chi^2 function. Thank you to both Professor Marla Geha and Imad Pasha for reading over and suggesting comments to my code. Thank you to Professor Meg Urry and Mike Koss at STSci for inspiring me to create this package.

Part of the functions to pull images directly from the web are from the MAST Queries example, at https://astroquery.readthedocs.io/en/latest/mast/mast.html. The way to align images to the catalog is taken from https://drizzlepac.readthedocs.io/en/latest/tweakreg.html. 

## Contributing
Pull requests are welcome. 

## License
[MIT](https://choosealicense.com/licenses/mit/)
