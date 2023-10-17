#%%
from astropy.io import fits
import numpy as np
from skimage.measure import label, regionprops
import Astrocode as ac
import matplotlib.pyplot as plt
import blooming as bl
import photometry as pt

#%%
hdulist = fits.open("test_image_2_circles_bloom.fits")
hdulist[0].header
#fits.getheader("Fits_data/mosaic.fits",0)


data = hdulist[0].data

labeled_data = label(data, connectivity = 2) #assigns unique label to each connected group of pixels
num_labels_binary = np.max(labeled_data) #number of objects in binary data
print(num_labels_binary)
object_properties = regionprops(labeled_data)

blooming_mask = np.zeros_like(data, dtype=bool) #creates mask - array of zeros same shape as binary_data

threshold_area = 90

for prop in object_properties:
    # Adjust the criteria as needed to identify blooming artifacts
    if prop.area > threshold_area:
        blooming_mask[labeled_data == prop.label] = True

cleaned_data = data.copy()
cleaned_data[blooming_mask] = 0  #applies mask
cleaned_data_label = label(cleaned_data, connectivity = 2)
num_labels_clean = np.max(cleaned_data_label) #number of objects in cleaned data
print(num_labels_clean)

output_file = 'test_image.fits'

# Create a FITS header
header = fits.Header()
header['COMMENT'] = 'test image data'

# Save the cleaned data to the FITS file
fits.writeto(output_file, cleaned_data, header=header, overwrite=True)


# %%

border_size = 1000   #cropping border pixels from edge of image
new_width = bl.final_clean_data.shape[1] - 2 * border_size  #removing pixels from both sides
new_height = bl.final_clean_data.shape[0] - 2 * border_size #removing pixels from top and bottom

test_trimmed_image = bl.final_clean_data[border_size:border_size + new_height, border_size:border_size + new_width]
test_background_image = pt.withbackground_data[border_size:border_size + new_height, border_size:border_size + new_width]

output_file = 'test_trimmed.fits'
fits.writeto(output_file, test_trimmed_image, header = header, overwrite=True)

#%%

# Object detection using label (you may use a more sophisticated method)
labeled_image = label(test_trimmed_image, connectivity=2)

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
    y, x = np.indices(test_trimmed_image.shape)
    aperture_mask = ((x - object_centroid[1])**2 + (y - object_centroid[0])**2 <= galaxy_aperture_radius**2)


    # Define annulus mask around the object (an annular region, for example)
    annulus_mask = ((x - object_centroid[1])**2 + (y - object_centroid[0])**2 >= background_inner_radius**2) & (
        (x - object_centroid[1])**2 + (y - object_centroid[0])**2 <= background_outer_radius**2)
    

    # Calculate the background level within the annulus
    background_counts = np.sum(test_background_image * annulus_mask)  # Total counts in annulus
    background_area = np.sum(annulus_mask)  # Number of pixels in annulus
    background_per_pixel = background_counts / background_area

    # Calculate the net counts for the galaxy
    galaxy_counts = np.sum(test_background_image * aperture_mask)  # Total counts in aperture
    galaxy_area = np.sum(aperture_mask)  # Number of pixels in aperture
    net_counts = galaxy_counts - (background_per_pixel * galaxy_area)

    # Store the galaxy information in the list
    galaxies_info.append({'Galaxy': region.label, 'NetCounts': net_counts}) #gives pixel value for each galaxy minus the background
 
    # Store or process the net counts for this galaxy

# Count the detected galaxies
number_of_galaxies = len(regionprops(labeled_image))


# %%

# Load your FITS data and perform necessary pre-processing steps
# ...

# Object detection using label (you may use a more sophisticated method)
labeled_image = label(test_trimmed_image, connectivity=2)

# List to store galaxy information (brightness)
galaxies_info_test = []

# Define aperture and annulus parameters
galaxy_aperture_radius = 6  # Adjust the radius as needed
background_inner_radius = 8
background_outer_radius = 15

# Initialize a list to store the net counts for each galaxy
net_counts_list = []

# Define a threshold value for pixel values in the annulus region
std = 2
threshold_background = 3418.13 + std * 11.79
threshold_aperture = 3418.13 + std * 11.79

# Iterate through labeled objects (galaxies)
for region in regionprops(labeled_image):
    # Extract the object's coordinates or centroid
    object_centroid = region.centroid

    # Create an annulus mask for this region
    annulus_mask = ((x - object_centroid[1])**2 + (y - object_centroid[0])**2 >= background_inner_radius**2) & (
        (x - object_centroid[1])**2 + (y - object_centroid[0])**2 <= background_outer_radius**2)

    # Exclude this region from the annulus mask
    annulus_mask[region.coords[:, 0], region.coords[:, 1]] = False

    # Define aperture mask around the object (a circular aperture, for example)
    aperture_mask = ((x - object_centroid[1])**2 + (y - object_centroid[0])**2 <= galaxy_aperture_radius**2)

    # Thresholding in the annulus region and masking
    annulus_mask_above_threshold = annulus_mask.copy()
    annulus_mask_above_threshold[test_background_image * annulus_mask > threshold_background] = False

     # Calculate the background level within the annulus
    background_counts = np.sum(test_background_image * annulus_mask_above_threshold)  # Total counts in annulus
    background_area = np.sum(annulus_mask)  # Number of pixels in annulus
    background_per_pixel = background_counts / background_area

    aperture_mask_below_threshold = aperture_mask.copy()
    aperture_mask_below_threshold[test_background_image * annulus_mask < threshold_aperture] = False

    # Calculate the net counts for the galaxy
    galaxy_counts = np.sum(test_background_image * aperture_mask)  # Total counts in aperture
    galaxy_area = np.sum(aperture_mask)  # Number of pixels in aperture
    net_counts = galaxy_counts - (background_per_pixel * galaxy_area)

    # Store the net counts for the galaxy
    net_counts_list.append(net_counts)
     # Store the galaxy information in the list
    galaxies_info_test.append({'Galaxy': region.label, 'NetCounts': net_counts}) #gives pixel value for each galaxy minus the background
 

# Count the detected galaxies
number_of_galaxies = len(regionprops(labeled_image))


# %%

import numpy as np
from skimage import measure

# Object detection using label (you may use a more sophisticated method)
labeled_image = measure.label(test_trimmed_image, connectivity=2)

# Initialize a list to store the net counts for each galaxy
net_counts_list = []

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
    y, x = np.indices(test_background_image.shape)
    aperture_radius_squared = (adaptive_aperture_radius / 2)**2
    aperture_mask = ((x - object_centroid[1])**2 + (y - object_centroid[0])**2 <= aperture_radius_squared)

    # Calculate the background annulus mask around the object
    annulus_mask = ((x - object_centroid[1])**2 + (y - object_centroid[0])**2 >= local_background_inner_radius**2) & (
        (x - object_centroid[1])**2 + (y - object_centroid[0])**2 <= local_background_outer_radius**2)

    # Apply the threshold within the adaptive aperture radius
    data_within_aperture = test_background_image * aperture_mask
    data_within_aperture[data_within_aperture < threshold_value] = 0

    # Calculate the local background within the annulus
    local_background_counts = np.sum(test_background_image * annulus_mask)  # Total counts in annulus
    local_background_area = np.sum(annulus_mask)  # Number of pixels in annulus
    local_background_per_pixel = local_background_counts / local_background_area

    # Subtract the local background from the galaxy counts
    galaxy_counts = np.sum(data_within_aperture)
    net_counts = galaxy_counts - local_background_per_pixel * np.sum(aperture_mask)

    # Store the net counts for the galaxy
    net_counts_list.append(net_counts)

# Count the detected galaxies
number_of_galaxies = len(measure.regionprops(labeled_image))

# Print or process the list of net counts
for i, net_counts in enumerate(net_counts_list):
    print(f'Galaxy {i + 1}: Net Counts = {net_counts}')

# %%
