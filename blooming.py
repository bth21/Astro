
from astropy.io import fits
import numpy as np
from scipy import ndimage
from skimage.measure import label, regionprops
import Astrocode as ac
import matplotlib.pyplot as plt


labeled_data = label(ac.binary_data, connectivity = 2) #assigns unique label to each connected group of pixels
num_labels_binary = np.max(labeled_data) #number of objects in binary data
print(num_labels_binary)
object_properties = regionprops(labeled_data)

blooming_mask = np.zeros_like(ac.binary_data, dtype=bool) #creates mask - array of zeros same shape as binary_data

threshold_area = 5000

for prop in object_properties:
    # Adjust the criteria as needed to identify blooming artifacts
    if prop.area > threshold_area:
        blooming_mask[labeled_data == prop.label] = True

cleaned_data = ac.binary_data.copy()
cleaned_data[blooming_mask] = 0  #applies mask
cleaned_data_label = label(cleaned_data, connectivity = 2)
num_labels_clean = np.max(cleaned_data_label) #number of objects in cleaned data
print(num_labels_clean)

output_file = 'cleaned_image.fits'

# Create a FITS header
header = fits.Header()
header['COMMENT'] = 'Cleaned image data'

# Save the cleaned data to the FITS file
fits.writeto(output_file, cleaned_data, header=header, overwrite=True)


ac.trimmed_image[ac.trimmed_image == 0] = 0 #replacing the background in the original image with zeros

restored_data = ac.trimmed_image * cleaned_data #combining the original data with the cleaned data
output_file = 'restored_image.fits'
restored_data_label = label(restored_data, connectivity = 2)
num_labels_restored = np.max(restored_data_label)
print(num_labels_restored)

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
