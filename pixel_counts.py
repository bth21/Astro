from astropy.io import fits
import numpy as np
from skimage.measure import label, regionprops
import Astrocode as ac
import blooming as bl
import matplotlib.pyplot as plt
#%%

properties = regionprops(bl.restored_data) # Use regionprops to calculate properties of labeled objects

'''removing background count is actually far more complex than this - code needs drastically
altering'''
background_value = 3418 #define background value from mean

corrected_restored_data = np.zeros_like(bl.restored_data) #creating empty array of same shape
for i in range(bl.restored_data.shape[0]):
    for j in range(bl.restored_data.shape[1]):
        if bl.restored_data[i, j] != 0:  #if pixel values are non-zero then minus background value
            corrected_restored_data[i, j] = bl.restored_data[i, j] - background_value
        else:
            corrected_restored_data[i, j] = 0

# Initialize a dictionary to store pixel counts for each labeled object
pixel_counts = {}

# Iterate through the labeled objects and calculate the pixel counts
for region in properties:
    label_value = region.label
    pixel_count = region.area
    pixel_counts[label_value] = pixel_count

# Use regionprops to calculate properties of labeled objects
properties = regionprops(bl.restored_data_label)


# %%
