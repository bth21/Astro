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
