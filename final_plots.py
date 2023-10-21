#%%
import photometry as pt
import separation as sp
import numpy as np
import matplotlib.pyplot as plt
#%%
# Convert magnitudes to a NumPy array
magnitudes = np.array(pt.cleaned_list)

# Define bin width and create bins
bin_width = 0.1
bin_edges = np.arange(0, 25 + bin_width, bin_width)

# Clip magnitudes to ensure they fall within the range of bin_edges
magnitudes = np.clip(magnitudes, bin_edges[0], bin_edges[-1])

# Calculate the number of galaxies in each bin
hist, bin_edges = np.histogram(magnitudes, bins=bin_edges)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])  # Calculate the center of each bin

# Create a scatter plot of number of galaxies vs. magnitude
plt.scatter(bin_centers, hist)
plt.xlabel('Magnitude')
plt.ylabel('Number of Galaxies')
plt.title('Number of Galaxies vs. Magnitude (Binned)')
plt.grid(True)

plt.show()

complete_mag_list = pt.cleaned_list.append(sp.calibrated_mag_list_sep)


# %%
