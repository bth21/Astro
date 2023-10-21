#%%
import csv
import photometry as pt
import separation as sp
import math
from skimage.measure import label, regionprops


"Cataloguing primary data in an ASCII file"


cleaned_info_dict = {f'Galaxy{i+1}': value for i, value in enumerate(pt.cleaned_list)}

# Create a list of dictionaries with the information you want
# Create a list of dictionaries with the information you want
object_data = []

for region in regionprops(pt.labeled_image):
    x, y = map(int, region.centroid)  # Get the x, y coordinates and convert to integers
    label = region.label

    # Check if the brightness value is available in the cleaned_info_dict
    brightness = cleaned_info_dict.get(f'Galaxy{label}', None)
    if brightness is None:
        # If not available, use the original galaxies_info dictionary
        brightness = pt.galaxies_info[label - 1]['NetCounts']

    total_pixel_count = len(region.coords)
    eccentricity = region.eccentricity
    area = region.area

    object_data.append({
        'x_coord': y,
        'y_coord': x,
        'Brightness': brightness,
        'Total_Pixel_Count': total_pixel_count,
        'Eccentricity': eccentricity,
        'Area': area
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
