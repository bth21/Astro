#%%
from astropy.io import fits
import numpy as np
from scipy import ndimage
from skimage.measure import label, regionprops
import Astrocode as ac
import matplotlib.pyplot as plt
import binary as bn
from skimage.morphology import binary_opening, square

'''removing blooming pixels'''

labeled_data = label(bn.binary_data, connectivity = 2) #assigns unique label to each connected group of pixels
num_labels_binary = np.max(labeled_data) #number of objects in binary data
print(num_labels_binary)
object_properties = regionprops(labeled_data)

blooming_mask = np.zeros_like(bn.binary_data, dtype=bool) #creates mask - array of zeros same shape as binary_data

threshold_area = 1800

for prop in object_properties:
    # Adjust the criteria as needed to identify blooming artifacts
    if prop.area > threshold_area:
        blooming_mask[labeled_data == prop.label] = True

cleaned_data = bn.binary_data.copy()
cleaned_data[blooming_mask] = 0  #applies mask
cleaned_data_label = label(cleaned_data, connectivity = 2)
num_labels_clean = np.max(cleaned_data_label) #number of objects in cleaned data
print(num_labels_clean)  #removing blooming

output_file = 'cleaned_image.fits'

# Create a FITS header
header = fits.Header()
header['COMMENT'] = 'Cleaned image data'

# Save the cleaned data to the FITS file
fits.writeto(output_file, cleaned_data, header=header, overwrite=True)
#%%
'''removing hot pixels and noise from very bright stars'''

labeled_data1 = label(bn.binary_data, connectivity = 2) #assigns unique label to each connected group of pixels
object_properties1 = regionprops(labeled_data1)
# Define a minimum area for objects to keep
min_area = 10  # Adjust this value based on your requirements

# Create a binary mask based on the object area criteria
mask = np.zeros_like(cleaned_data, dtype=bool)

for region in object_properties1:
    if region.area >= min_area:
        # Set the pixels corresponding to the object in the mask to True
        mask[labeled_data1 == region.label] = True

# Apply the mask to your FITS image
final_clean_data = cleaned_data.copy()
final_clean_data[~mask] = 0  # Set pixels outside the mask to 0

output_file = 'squeakyclean.fits'
fits.writeto(output_file, final_clean_data, header = header, overwrite=True)

labels_clean_noise = label(final_clean_data, connectivity = 2)
num_labels_clean_noise = np.max(labels_clean_noise) #number of objects in cleaned data
print(num_labels_clean_noise) #number of objects after removing lone pixels and blooming



#%%
'''combining clean binary data with original data to get pixel values'''
ac.trimmed_image[ac.trimmed_image == 0] = 0 #replacing the background in the original image with zeros

restored_data = ac.trimmed_image * final_clean_data #combining the original data with the cleaned data
output_file = 'restored_image.fits'
restored_data_label = label(restored_data, connectivity = 2)
num_labels_restored = np.max(restored_data_label)
print(num_labels_restored) #number of objects with blooming and hot pixels removed

header['COMMENT'] = 'Restored image data'

# Save the restored data to the FITS file
fits.writeto(output_file, restored_data, header=header, overwrite=True)


# %%
low_r = 3000  #defining mask parameters to plot histogram of pixel data
high_r = 3800
mask = (restored_data >= low_r) & (restored_data <= high_r)
restored_data_hist = restored_data[mask]

plt.hist(restored_data_hist, bins = 500)  #plotting histogram of pixel data for restored image
plt.xlabel('Pixel Value')
plt.ylabel('Frequency')
#plt.savefig('Restored_hist')
plt.show()

# %%
