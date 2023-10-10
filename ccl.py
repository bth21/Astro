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

# Create a named window (with the WINDOW_NORMAL flag)
cv2.namedWindow('Binary Image', cv2.WINDOW_NORMAL)

# Set the window to be resizable
cv2.setWindowProperty('Binary Image', cv2.WND_PROP_FULLSCREEN, cv2.WINDOW_NORMAL)

# Show the binary image
#cv2.imshow('Binary Image', binary_image)

# Perform connected component labeling
num_labels, labeled_image = cv2.connectedComponents(binary_image)

# Show labelled image
cv2.imshow('Binary Image', binary_image)
print(num_labels)

# Close window
cv2.waitKey(0)
cv2.destroyAllWindows() 




