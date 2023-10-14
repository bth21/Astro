"""This module plots the histogram of the pixel values"""



from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import Astrocode as ac
import filterdata as fd



# Specify the number of bins in histogram
num_bins = 3500


# Plot the histogram
plt.hist(fd.filtered_data2, bins=num_bins, density=False, alpha=0.75, color='blue')

# Add labels and a title
plt.xlabel('Pixel Value')
plt.ylabel('Frequency')
plt.title('Histogram of Pixel Values')


# Show the plot
plt.show()

#low_r = 3000  #defining mask parameters to plot histogram of pixel data
#high_r = 3800
#mask = (restored_data >= low_r) & (restored_data <= high_r)
#restored_data_hist = restored_data[mask]

#plt.hist(restored_data_hist, bins = 500)  #plotting histogram of pixel data for restored image
#plt.xlabel('Pixel Value')
#plt.ylabel('Frequency')
#plt.savefig('Restored_hist')
#plt.show()
# %%