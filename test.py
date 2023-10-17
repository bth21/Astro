#%%
from astropy.io import fits
import numpy as np
from scipy import ndimage
from skimage.measure import label, regionprops
import Astrocode as ac
import matplotlib.pyplot as plt

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
