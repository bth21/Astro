#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 10:23:23 2023

@author: benhale
"""
<<<<<<< HEAD
print("getting gitty with it")

print("Github sucks")
=======
#%%
from astropy.io import fits
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import ndimage
from skimage.measure import label, regionprops


hdulist = fits.open("Fits_data/mosaic.fits")
hdulist[0].header
#fits.getheader("Fits_data/mosaic.fits",0)


data = hdulist[0].data

border_size = 130   #cropping border pixels from edge of image
new_width = data.shape[1] - 2 * border_size  #removing pixels from both sides
new_height = data.shape[0] - 2 * border_size #removing pixels from top and bottom

trimmed_image = data[border_size:border_size + new_height, border_size:border_size + new_width]

# Create a mask to filter values
low_range = 3350
high_range = 3500
mask = (trimmed_image >= low_range) & (trimmed_image <= high_range)


# Apply the mask to your data to exclude values less than 3500 and greater than 3350
filtered_data = trimmed_image[mask]

# Specify the number of bins you want in your histogram
num_bins = 50

# Plot the histogram
plt.hist(filtered_data.flatten(), bins = num_bins, alpha=0.75, color='blue')

# Add labels and a title
plt.xlabel('Pixel Value')
plt.ylabel('Frequency')
plt.title('Histogram')


hist, bin_edges = np.histogram(filtered_data, bins = num_bins)  #finding peak value of the histogram with 1000 bins
peak_bin_index = np.argmax(hist)
peak_value = (bin_edges[peak_bin_index] * 2 + 1) / 2 #finding centre of bin
#peak_value is 3418.25

x_centres = (bin_edges[:-1] + bin_edges[1:]) / 2
y_centres = hist

def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-(x - mean)**2 / (2 * stddev**2))

# Initial guess for the parameters (amplitude, mean, and standard deviation)
initial_guess = [800000, 3418, np.std(x_centres)]

# Perform the curve fitting
fit_params, _ = curve_fit(gaussian, x_centres, y_centres, p0=initial_guess)

amplitude_fit, mean_fit, stddev_fit = fit_params

plt.scatter(x_centres, y_centres, label='Data', color='red')
plt.plot(x_centres, gaussian(x_centres, *fit_params), label='Fit', color='red')

plt.xlabel('Pixel Value')
plt.ylabel('Frequency')
#plt.title('Fitting a Gaussian Curve to Data')
plt.legend()
#plt.savefig('Astrohist_fit')

# Show the plot
#plt.show()



# %%
>>>>>>> Ben
