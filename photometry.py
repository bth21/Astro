
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
# Object detection using label
labeled_image = measure.label(bl.final_clean_data, connectivity=2)

# Initialize a list to store the net counts for each galaxy
net_counts_list = []
galaxies_info = []

# Define a threshold value for pixel values within the aperture radius
std = 2
threshold_value = 3418.13 + std * 11.79

# Define the size of the local background annulus
local_background_inner_radius = 12  # Adjust as needed
local_background_outer_radius = 15  # Adjust as needed

# Iterate through labeled objects (galaxies)
for region in measure.regionprops(labeled_image):
    # Extract the object's coordinates or centroid
    object_centroid = region.centroid

    # Calculate the adaptive aperture radius based on the object's characteristics
    major_axis_length = region.major_axis_length
    adaptive_aperture_radius = 0.5 * major_axis_length  # Adjust the factor as needed

    # Define the aperture mask around the object with the adaptive radius
    y, x = np.indices(withbackground_data.shape)
    aperture_radius_squared = (adaptive_aperture_radius / 2)**2
    aperture_mask = ((x - object_centroid[1])**2 + (y - object_centroid[0])**2 <= aperture_radius_squared)

    # Calculate the background annulus mask around the object
    annulus_mask = ((x - object_centroid[1])**2 + (y - object_centroid[0])**2 >= local_background_inner_radius**2) & (
        (x - object_centroid[1])**2 + (y - object_centroid[0])**2 <= local_background_outer_radius**2)

    # Apply the threshold within the adaptive aperture radius
    data_within_aperture = withbackground_data * aperture_mask
    data_within_aperture[data_within_aperture < threshold_value] = 0

    # Calculate the local background within the annulus
    local_background_counts = np.sum(withbackground_data * annulus_mask)  # Total counts in annulus
    local_background_area = np.sum(annulus_mask)  # Number of pixels in annulus
    local_background_per_pixel = local_background_counts / local_background_area

    # Subtract the local background from the galaxy counts
    galaxy_counts = np.sum(data_within_aperture)
    net_counts = galaxy_counts - local_background_per_pixel * np.sum(aperture_mask)

    # Store the net counts for the galaxy
    net_counts_list.append(net_counts)
    # Store the galaxy information in the list
    galaxies_info.append({'Galaxy': region.label, 'NetCounts': net_counts}) #gives pixel value for each galaxy minus the background
 

# Count the detected galaxies
number_of_galaxies = len(measure.regionprops(labeled_image))

# %%

hdul = fits.open("Fits_data/mosaic.fits")
header = hdul[0].header
magzpt_value = header['MAGZPT']

def mag(c, z):
    m = z - 2.5 * np.log10(c)
    return m

net_counts_list = [galaxy['NetCounts'] for galaxy in galaxies_info]

calibrated_mag_list = []
for i in net_counts_list:
    calibrated_mag_list.append(mag(i, magzpt_value))

cleaned_list = []
for value in calibrated_mag_list:
    if not math.isnan(value):
        cleaned_list.append(value)

print(cleaned_list)

# %%
mask = np.zeros_like(labeled_image)

for galaxy_info in galaxies_info:
    if galaxy_info["NetCounts"] >= 0:
        # Set the pixels corresponding to the positive object to 1 in the mask
        mask[labeled_image == galaxy_info["Galaxy"]] = 1

# Apply the mask to your data to mask the positive objects
masked_data = bl.final_clean_data.copy()
masked_data[mask == 1] = 0  # Set positive object pixels to 0 or any desired value

fits.writeto("masked_data.fits", masked_data, overwrite=True)

# Perform further analysis or save the masked data
# %%
