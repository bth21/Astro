#%%
from astropy.io import fits
import numpy as np
import Astrocode as ac
from skimage.measure import label

'''Converting pixel data into binary data with 2stds above mean
background value'''
# Define the threshold value (3442 is 2std above mean)
std = 5
threshold = 3418.13 + std * 11.79

# Apply the threshold
binary_data = (ac.trimmed_image > threshold).astype(np.uint8)


# Close the FITS file
ac.hdulist.close()

# Save the binary data to a binary file (e.g., a FITS file or a binary image)
# Here, we save it as a new FITS file
binary_image = fits.PrimaryHDU(binary_data)
binary_image.writeto('binary_output.fits', overwrite=True)


# %%
