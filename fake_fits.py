import numpy as np
from astropy.io import fits

# Define the dimensions of your image (e.g., 100x100 pixels)
image_shape = (100, 100)

# Create a random 2D array with pixel values between 0 and 255
#image_data = np.random.randint(0, 256, size=image_shape, dtype=np.uint8)
pixel_values = np.zeros(image_shape)

import numpy as np

# Define the center coordinates and radius of 2 circles
center_x_1, center_y_1 = 20, 50  
radius_1 = 5  

center_x_2, center_y_2 = 80, 50  
radius_2 = 5  

# Create a vertical line for blooming
start_x = 20
start_y = 50-15
vert_length = 30
width = 1 


# Create a mesh grid of coordinates
x, y = np.meshgrid(np.arange(100), np.arange(100))

# Calculate the distance from each pixel to the center for each circle
distance_1 = np.sqrt((x - center_x_1)**2 + (y - center_y_1)**2)
distance_2 = np.sqrt((x - center_x_2)**2 + (y - center_y_2)**2)

# Calculate the ending X coordinate
end_x = start_x + width
end_y = start_y + vert_length

# Set the values inside the circle(s) to 251
pixel_values[distance_1 <= radius_1] = 251
pixel_values[distance_2 <= radius_2] = 251
pixel_values[start_y:end_y, start_x:end_x] = 251

# Create a PrimaryHDU object to store the image data and header
hdu = fits.PrimaryHDU(data=pixel_values)

# Define the filename for your FITS file
fits_filename = 'test_image_2_circles_bloom.fits'

# Write the FITS file
hdu.writeto(fits_filename, overwrite=True)

