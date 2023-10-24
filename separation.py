#%%
import numpy as np
from astropy.io import fits
import numpy as np
from scipy import ndimage
from skimage.measure import label, regionprops
import Astrocode as ac
import matplotlib.pyplot as plt
import binary as bn
import blooming as bl
from skimage import measure
import math
import photometry as pt
# %%
'''attempt to separate joined galaxies'''

# Label the objects in the data
labeled_galx = measure.label(pt.masked_data, connectivity=2)

# Analyze region properties
regions = measure.regionprops(labeled_galx)

# Define thresholds for perimeter-to-area ratio and minimum area
threshold_ratio = 0.1  # Adjust as needed
min_area = 70  # Adjust as needed

# Create a binary mask to retain objects that meet the criteria
binary_mask_galx = np.zeros_like(pt.masked_data, dtype=bool)

for region in regions:
    perimeter = region.perimeter
    area = region.area

    # Calculate the perimeter-to-area ratio
    ratio = perimeter / area

    # Check if the object meets the criteria
    if ratio > threshold_ratio and area > min_area:
        # Object meets criteria, retain it in the mask
        binary_mask_galx[labeled_galx == region.label] = True

# Apply the binary mask to the original data
sep_data = np.where(binary_mask_galx, pt.masked_data, 0)
sep_data_label = label(sep_data, connectivity = 2)

# Save the masked data as a new FIT file
fits.writeto("twogalax_data.fits", sep_data, overwrite=True)

# %%
non_binary = sep_data_label * ac.trimmed_image
# Get region properties for the labeled objects
pixel_regions = measure.regionprops(sep_data_label)

# Initialize a list to store the sum of pixel values for each object
sum_pixel_values = []

# Iterate through each labeled object and calculate the sum of pixel values
for region in pixel_regions:
    object_coords = region.coords  # Coordinates of the object
    pixel_values = non_binary[object_coords[:, 0], object_coords[:, 1] ]
    total_pixel_sum = np.sum(pixel_values)  # Sum of pixel values within the object
    sum_pixel_values.append(total_pixel_sum)
# %%

'''removing local background from pixel sum'''
# Initialize a list to store the corrected sum of pixel values for each object
corrected_sum_pixel_values = []

# Define parameters for the background annulus
background_inner_radius = 10  # Adjust as needed
background_outer_radius = 12  # Adjust as needed

# Iterate through each labeled object and calculate the corrected sum of pixel values
for region in pixel_regions:
    object_coords = region.coords  # Coordinates of the object

    # Calculate the annulus mask for the background
    y, x = np.indices(pt.withbackground_data.shape)
    object_centroid = region.centroid
    annulus_mask = (
        (x - object_centroid[1])**2 + (y - object_centroid[0])**2 >= background_inner_radius**2
    ) & (
        (x - object_centroid[1])**2 + (y - object_centroid[0])**2 <= background_outer_radius**2
    )

    # Extract pixel values within the annulus
    annulus_pixel_values = pt.withbackground_data[annulus_mask]

    # Calculate the mean background pixel value
    mean_background_value = np.mean(annulus_pixel_values)

    # Extract pixel values within the object
    object_pixel_values = pt.withbackground_data[object_coords[:, 0], object_coords[:, 1]]

    # Calculate the sum of pixel values within the object
    total_pixel_sum = np.sum(object_pixel_values)

    # Corrected sum of pixel values by subtracting the mean background
    corrected_sum_value = total_pixel_sum - (mean_background_value * len(object_pixel_values))

    corrected_sum_pixel_values.append(corrected_sum_value)

calibrated_mag_list_sep = []
for i in corrected_sum_pixel_values:
    calibrated_mag_list_sep.append(pt.mag(i, pt.magzpt_value))

