#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 10:23:23 2023

@author: benhale
"""
import numpy as np
import scipy as scipy
import matplotlib.pyplot as plt
from astropy.io import fits

hdulist = fits.open("Fits_data/mosaic.fits")
hdulist[0].header
#fits.getheader("Fits_data/mosaic.fits",0)


data = hdulist[0].data

num_bins = 10
plt.hist(data.flatten(), num_bins)
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.show()


#%%


border_size = 130   #cropping border pixels from edge of image
new_width = data.shape[1] - 2 * border_size
new_height = data.shape[0] - 2 * border_size

trimmed_image = data[border_size:border_size + new_height, border_size:border_size + new_width]


# Assuming you have your data in the variable 'data'
# Create a mask to filter values greater than or equal to 10,000
mask = trimmed_image <= 3700

# Apply the mask to your data to exclude values less than 10,000
filtered_data = trimmed_image[mask]

# Specify the number of bins you want in your histogram
num_bins = 50

# Plot the histogram
plt.hist(filtered_data.flatten(), bins=num_bins, density=True, alpha=0.75, color='blue')

# Add labels and a title
plt.xlabel('Pixel Value')
plt.ylabel('Frequency')
plt.title('Histogram (Values >= 10,000)')

# Show the plot
plt.show()