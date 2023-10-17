
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
"Cataloging data to ASCII file"


# Define the file path for the existing ASCII file
file_path = "galaxies_data.txt"

# Open the file in append mode and write the new data
with open(file_path, "a") as file:
    # Write data for the new galaxy
    file.write(f"{new_galaxy['name']}\t{new_galaxy['type']}\t{new_galaxy['magnitude']}\n")






#%%
# Object detection using label (you may use a more sophisticated method)
labeled_image = label(bl.final_clean_data, connectivity=2)

# List to store galaxy information (brightness)
galaxies_info = []

# Define aperture and annulus parameters
galaxy_aperture_radius = 6  # Adjust the radius as needed
background_inner_radius = 10
background_outer_radius = 12

# Iterate through labeled objects (galaxies)
for region in regionprops(labeled_image):
    # Extract the object's coordinates or centroid
    object_coords = region.coords
    object_centroid = region.centroid

   # Define aperture mask around the object (a circular aperture, for example)
    y, x = np.indices(withbackground_data.shape)
    aperture_mask = (
        (x - object_centroid[1])**2 + (y - object_centroid[0])**2 <= galaxy_aperture_radius**2)

    # Define annulus mask around the object (an annular region, for example)
    annulus_mask = (
        (x - object_centroid[1])**2 + (y - object_centroid[0])**2 >= background_inner_radius**2) & (
        (x - object_centroid[1])**2 + (y - object_centroid[0])**2 <= background_outer_radius**2)


    # Calculate the background level within the annulus
    background_counts = np.sum(withbackground_data * annulus_mask)  # Total counts in annulus
    background_area = np.sum(annulus_mask)  # Number of pixels in annulus
    background_per_pixel = background_counts / background_area

    # Calculate the net counts for the galaxy
    galaxy_counts = np.sum(withbackground_data * aperture_mask)  # Total counts in aperture
    galaxy_area = np.sum(aperture_mask)  # Number of pixels in aperture
    net_counts = galaxy_counts - (background_per_pixel * galaxy_area)

    # Store the galaxy information in the list
    galaxies_info.append({'Galaxy': region.label, 'NetCounts': net_counts}) #gives pixel value for each galaxy minus the background
 

# Count the detected galaxies
number_of_galaxies = len(regionprops(labeled_image))

#%%
# Object detection using label (you may use a more sophisticated method)
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

calibrated_mag = []
for value in calibrated_mag_list:
    if value == float:
        calibrated_mag.append(calibrated_mag_list)



