#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 10:23:23 2023

@author: benhale
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# open astro data file and create an array of pixel values
hdulist = fits.open(r"C:\Users\Administrator\GitHub\Astro\mosaic.fits")
hdulist[0].header
data = hdulist[0].data

# cropping border pixels from edge of image
border_size = 130   
new_width = data.shape[1] - 2 * border_size
new_height = data.shape[0] - 2 * border_size

trimmed_image = data[border_size:border_size + new_height, border_size:border_size + new_width]

# Assuming you have your data in the variable 'data'
# Create a mask to filter values 
low_range = 3350
high_range = 3600
mask = (trimmed_image >= low_range) & (trimmed_image <= high_range) 

# Apply the mask to your data to exclude values 
filtered_data = trimmed_image[mask]
#print(filtered_data)


# Specify the number of bins you want in your histogram
num_bins = 350


# Plot the histogram
plt.hist(filtered_data, bins=num_bins, density=False, alpha=0.75, color='blue')

# Add labels and a title
plt.xlabel('Pixel Value')
plt.ylabel('Frequency')
plt.title('Histogram of Pixel Values')


# Show the plot
plt.show()



