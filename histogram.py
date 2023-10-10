"""This module plots the histogram of the pixel values"""


from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import Astrocode as ac
import filterdata as fd



# Specify the number of bins you want in your histogram
num_bins = 350


# Plot the histogram
plt.hist(fd.filtered_data, bins=num_bins, density=False, alpha=0.75, color='blue')

# Add labels and a title
plt.xlabel('Pixel Value')
plt.ylabel('Frequency')
plt.title('Histogram of Pixel Values')


# Show the plot
plt.show()