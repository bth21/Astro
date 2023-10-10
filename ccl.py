"""This module carries out Connected Component Labelling (CCL)"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import Astrocode as ac
import filterdata as fd
import cv2
import numpy as np


# Apply binary thresholding
binary_image = (ac.data > fd.std_filter).astype('uint8') * 255

# Show the binary image
cv2.imshow('Binary Image', binary_image)

# Close window
cv2.waitKey(0)
cv2.destroyAllWindows() 

