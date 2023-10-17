import numpy as np
from astropy.io import fits
import numpy as np
from scipy import ndimage
from skimage.measure import label, regionprops
import Astrocode as ac
import matplotlib.pyplot as plt
import binary as bn
import blooming as bl

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

# Calculate background_area once (constant for all galaxies)
background_annulus = (
    (x - object_centroid[1])**2 + (y - object_centroid[0])**2 >= background_inner_radius**2) & (
    (x - object_centroid[1])**2 + (y - object_centroid[0])**2 <= background_outer_radius**2)
background_area = np.sum(background_annulus)

# Iterate through labeled objects (galaxies)
for region in regionprops(labeled_image):
    # Extract the object's centroid
    object_centroid = region.centroid

    # Define aperture mask around the object (a circular aperture)
    aperture_mask = (
        (x - object_centroid[1])**2 + (y - object_centroid[0])**2 <= galaxy_aperture_radius**2)

    # Calculate the background level within the annulus
    background_counts = np.sum(withbackground_data * background_annulus)
    background_per_pixel = background_counts / background_area

    # Calculate the net counts for the galaxy
    galaxy_counts = np.sum(withbackground_data * aperture_mask)
    galaxy_area = np.sum(aperture_mask)
    net_counts = galaxy_counts - (background_per_pixel * galaxy_area)

    # Store the galaxy information in the list
    galaxies_info.append({'Galaxy': region.label, 'NetCounts': net_counts})

#%%
"Creating an ASCII file"

# Define the file path for the existing ASCII file
file_path = "galaxies_data.txt"

# Open the file in append mode and write the new data
with open(file_path, "a") as file:
    # Write data for the new galaxy
    file.write(f"{new_galaxy['name']}\t{new_galaxy['bright']}\t{new_galaxy['magnitude']}\n")



    
