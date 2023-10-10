"This module is used to filter the data"

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import Astrocode as ac


# cropping border pixels from edge of image
border_size = 130   
new_width = ac.data.shape[1] - 2 * border_size
new_height = ac.data.shape[0] - 2 * border_size

trimmed_image = ac.data[border_size:border_size + new_height, border_size:border_size + new_width]

# Assuming you have your data in the variable 'data'
# Create a mask to filter values 
low_range = 3350
high_range = 3600
mask1 = (trimmed_image >= low_range) & (trimmed_image <= high_range) 

# Apply the mask to your data to exclude values 
filtered_data1 = trimmed_image[mask1]

# Filter the data to select pixels > 2 std above mean
std_filter = 3418.13 + (5*11.79)
mask2 = trimmed_image >= std_filter
filtered_data2 = trimmed_image[mask2]



