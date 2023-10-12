import numpy as np
from astropy.io import fits

# Define the dimensions of your image (e.g., 100x100 pixels)
image_shape = (100, 100)

# Create a random 2D array with pixel values between 0 and 255
#image_data = np.random.randint(0, 256, size=image_shape, dtype=np.uint8)
pixel_values = np.zeros(image_shape)

# Create a PrimaryHDU object to store the image data and header
hdu = fits.PrimaryHDU(data=pixel_values)

# Define the filename for your FITS file
fits_filename = 'test_image_zeros.fits'

# Write the FITS file
hdu.writeto(fits_filename, overwrite=True)

