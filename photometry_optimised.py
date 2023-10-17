import numpy as np
from astropy.io import fits
import numpy as np
from scipy import ndimage
from skimage.measure import label, regionprops
import Astrocode as ac
import matplotlib.pyplot as plt
import binary as bn
import blooming as bl
import csv

#%%
"Defining a function to calculate the brightness from the net counts"
hdul = fits.open("Fits_data/mosaic.fits")
header = hdul[0].header
magzpt_value = header['MAGZPT']

def mag(c, z):
    m = z - 2.5 * np.log10(c)
    return m

#%%
'''Converting pixel data into binary data with 2stds above mean
background value to regain only the background'''
# Define the threshold value (3442 is 2std above mean)
std = 2
threshold = 3418.13 + std * 11.79

# Apply the threshold
antibinary_data = (ac.trimmed_image < threshold).astype(np.uint8)

antibinary_image = fits.PrimaryHDU(antibinary_data)
antibinary_image.writeto('antibinary_output.fits', overwrite=True)

antitrimmed_image =  antibinary_data * ac.trimmed_image #getting only background pixel counts

withbackground_data = bl.restored_data + antitrimmed_image

header = fits.Header()
header['COMMENT'] = 'Restored image data'
output_file = 'with_background_image.fits'

# Save the restored data to the FITS file
fits.writeto(output_file, withbackground_data, header=header, overwrite=True)

#%%
"Finding the net counts for each galaxy"

import numpy as np
from skimage import measure

# Object detection using label
labeled_image = measure.label(bl.final_clean_data, connectivity=2)

# Initialize a list to store the net counts for each galaxy
galaxies_info = []

# Define parameters
std = 2
threshold_value = 3418.13 + std * 11.79
local_background_inner_radius = 12
local_background_outer_radius = 15

# Calculate a thresholded image within the aperture
thresholded_data = np.where(withbackground_data >= threshold_value, withbackground_data, 0)

# Iterate through labeled objects (galaxies)
for region in measure.regionprops(labeled_image):
    # Extract the object's centroid and characteristics
    object_centroid = region.centroid
    major_axis_length = region.major_axis_length

    # Calculate the adaptive aperture radius
    adaptive_aperture_radius = 0.5 * major_axis_length

    # Define the aperture mask with the adaptive radius
    y, x = np.indices(withbackground_data.shape)
    aperture_mask = ((x - object_centroid[1])**2 + (y - object_centroid[0])**2 <= (adaptive_aperture_radius / 2)**2)

    # Calculate the background annulus mask
    annulus_mask = ((x - object_centroid[1])**2 + (y - object_centroid[0])**2 >= local_background_inner_radius**2) & (
        (x - object_centroid[1])**2 + (y - object_centroid[0])**2 <= local_background_outer_radius**2)

    # Calculate local background and galaxy counts
    local_background_counts = np.sum(withbackground_data * annulus_mask)
    local_background_area = np.sum(annulus_mask)
    local_background_per_pixel = local_background_counts / local_background_area
    galaxy_counts = np.sum(thresholded_data * aperture_mask)

    # Calculate net counts for the galaxy
    net_counts = galaxy_counts - local_background_per_pixel * np.sum(aperture_mask)

    brightness = mag(net_counts,magzpt_value)

    # Store the galaxy information in the list
    galaxies_info.append({'Galaxy': region.label,'Coordinates': region.centroid, 'NetCounts': net_counts,'Brightness': brightness, 'Area': region.area , 'Eccentricity': region.eccentricity})

    
#%%
"Cataloguing data in an ASCII file"

# Specify the file path where you want to save the ASCII file
file_path = "galaxies_info.csv"

# Define the fieldnames (column headers) for your CSV file
fieldnames = ['Galaxy', 'NetCounts', 'Brightness', 'Area', 'Coordinates', 'Eccentricity']

# Open the file in write mode and write the data
with open(file_path, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    
    # Write the headers
    writer.writeheader()
    
    # Write the data for each galaxy
    for galaxy in galaxies_info:
        writer.writerow(galaxy)





    

# %%
