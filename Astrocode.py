#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 10:23:23 2023

@author: benhale & ellalitjens
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# open astro data file and create an array of pixel values
hdulist = fits.open(r"C:\Users\Administrator\GitHub\Astro\mosaic.fits")
hdulist[0].header
data = hdulist[0].data









