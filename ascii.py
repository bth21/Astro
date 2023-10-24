#%%
import csv
import photometry as pt
import separation as sp
import math
from skimage.measure import label, regionprops
import numpy as np
from astropy.io import fits

hdul = fits.open("Fits_data/mosaic.fits")
header = hdul[0].header
magzpt_error = header['MAGZRR']

hdul = fits.open("Fits_data/mosaic.fits")
header = hdul[0].header
magzpt_value = header['MAGZPT']
# Cataloguing primary data in an ASCII file

cleaned_info_dict = {f'Galaxy{i+1}': value for i, value in enumerate(pt.cleaned_list)}

# Create a list of dictionaries with the information you want
object_data = []
errors_pix = []

for region in regionprops(pt.labeled_image):
    x, y = map(int, region.centroid)  # Get the x, y coordinates and convert to integers
    label = region.label

    # Check if the brightness value is available in the cleaned_info_dict
    brightness = cleaned_info_dict.get(f'Galaxy{label}', None)
    if brightness is None:
        # If not available, use the original galaxies_info dictionary
        brightness = pt.galaxies_info[label - 1]['NetCounts']

    pixel_values = [pt.withbackground_data[c[0], c[1]] for c in region.coords]  # Get the pixel values for the object
    total_pixel_count = len(region.coords)
    total_pixel_values = np.sum(pixel_values)
    pixel_std = np.std(pixel_values) * 100  # Calculate the standard deviation of pixel values
    eccentricity = region.eccentricity
    area = region.area
    error = (pixel_std/total_pixel_values) * brightness

    ########### star and colour classification ############

    solidity = region.solidity
    max_solidity = 0.95
    min_eccentricity = 0.4
    classification = 'empty'
    min_brightness = 20 # change to appropriate values
    min_solidty = 0.8 # change to appropriate values

    if solidity > max_solidity and eccentricity < min_eccentricity:
        classification = 'star'
    if brightness < min_brightness and solidity < min_solidty:
        classification = 'blue galaxy'
    elif brightness > min_brightness and solidity > min_solidty:
        classification = 'red galaxy'
    elif solidity > max_solidity and eccentricity < min_eccentricity:
        classification = 'star'
    else: 
        classification = 'N/A'



    #######################################################

   # print(pixel_std)
    object_data.append({
        'x_coord': y,
        'y_coord': x,
        'Brightness': brightness,
        'Total_Pixel_Count': total_pixel_count,
        'Pixel_Std': pixel_std,  # Include the standard deviation in the data
        'Eccentricity': eccentricity,
        'Area': area,
        'Error': error,
        'Classification': classification
    })

# Define the name of the output file
output_file = 'object_catalog.csv'

# Write the data to the output CSV file
with open(output_file, 'w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=object_data[0].keys())

    # Write the header row
    writer.writeheader()

    # Write data for each object
    writer.writerows(object_data)



# %%
'''catalouging secondary data into ASCII file'''

sec_cleaned_info_dict = {f'Galaxy{i+1}': value for i, value in enumerate(sp.calibrated_mag_list_sep)}

# Create a list of dictionaries with the information you want
# Create a list of dictionaries with the information you want
sec_object_data = []

for region in regionprops(sp.sep_data_label):
    x, y = map(int, region.centroid)  # Get the x, y coordinates and convert to integers
    label = region.label

    # Check if the brightness value is available in the cleaned_info_dict
    brightness = sec_cleaned_info_dict.get(f'Galaxy{label}', None)
    
    total_pixel_count = len(region.coords)
    eccentricity = region.eccentricity
    area = region.area

    sec_object_data.append({
        'x_coord': y,
        'y_coord': x,
        'Brightness': brightness,
        'Total_Pixel_Count': total_pixel_count,
        'Eccentricity': eccentricity,
        'Area': area
        
    })

# Define the name of the output file
sec_output_file = 'secondary_object_catalog.csv'
# Write the data to the output CSV file
with open(sec_output_file, 'w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=object_data[0].keys())

    # Write the header row
    writer.writeheader()

    # Write data for each object
    writer.writerows(sec_object_data)

# %%
