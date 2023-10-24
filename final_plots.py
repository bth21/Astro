#%%
import photometry as pt
import separation as sp
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
#%%
hdul = fits.open("Fits_data/mosaic.fits")
header = hdul[0].header
magzpt_error = header['MAGZRR']

error_data = 2 * np.loadtxt('errors.csv', skiprows = 1, delimiter = ',')

primary_magnitudes = pt.cleaned_list
secondary_magnitudes = sp.calibrated_mag_list_sep

# Define the bin width and range
bin_width = 0.1
bin_range = (5, 25)
sec_bin_range = (10, 15)

# Calculate the histogram using np.histogram
num_counts, bin_edges = np.histogram(primary_magnitudes, bins=np.arange(bin_range[0], bin_range[1] + bin_width, bin_width))
sec_num_counts, sec_bin_edges = np.histogram(secondary_magnitudes, bins=np.arange(sec_bin_range[0], sec_bin_range[1] + bin_width, bin_width))


# Calculate the bin centers
mag_bins = (bin_edges[:-1] + bin_edges[1:]) / 2
sec_mag_bins = (sec_bin_edges[:-1] + sec_bin_edges[1:]) / 2

full_counts = num_counts > 0
filtered_bin_centers = mag_bins[full_counts]
filtered_bin_counts = num_counts[full_counts]

sec_full_counts = sec_num_counts > 0
sec_filtered_bin_centers = sec_mag_bins[sec_full_counts]
sec_filtered_bin_counts = sec_num_counts[sec_full_counts]

bin_error_std = []
for i in range(len(filtered_bin_centers)):
    # Assuming you have a list of errors for each data point, e.g., error_data
    errors_in_bin = error_data[(primary_magnitudes >= bin_edges[i]) & (primary_magnitudes < bin_edges[i + 1])]
    bin_error_std.append(np.std(errors_in_bin))

bin_error_std = np.array(bin_error_std)

# Plot the scatter plot with error bars



# Plot the scatter plot
plt.figure(figsize=(10, 6))
plt.scatter(filtered_bin_centers, filtered_bin_counts, marker='o', color='b', label='Number of Galaxies')
plt.scatter(sec_filtered_bin_centers, sec_filtered_bin_counts, marker='o', color='orange', label='Number of Separated Galaxies')
plt.errorbar(filtered_bin_centers, filtered_bin_counts, xerr=bin_error_std, fmt='o', color='b', label='Number of Galaxies', capsize = 5, linestyle = 'none')

plt.xlabel('Magnitude')
plt.ylabel('Number of Galaxies')
plt.title('Scatter Plot of Number of Galaxies vs. Magnitude')
plt.grid(True)
plt.legend()
plt.show()

# %%

x = np.linspace(9.5, 18.5, 20, endpoint =True)
def why(x):
    y = (0.6 * (x)) - 6
    return y
y = why(x) 


log_num_counts = [np.log10(i) for i in filtered_bin_counts]
plt.figure(figsize=(10, 6))
plt.scatter(filtered_bin_centers, log_num_counts, marker='o', color='b', label='Number of Galaxies')
#plt.plot(x, y)
plt.xlabel('Magnitude')
plt.ylabel('log(N(m))')
plt.title('Scatter Plot of log(Number of Galaxies) vs. Magnitude')
plt.grid(True)
plt.legend()
plt.show()


# %%
'''cumulative graph'''

brightness = np.loadtxt('brightness_data.csv', delimiter=',', skiprows=1)
sorted_brightness = np.sort(brightness)


cumulative_counts = np.arange(1, len(sorted_brightness) + 1)

# Calculate log(N)
log_counts = np.log10(cumulative_counts)  # Use np.log for natural logarithm

# Plot log(N) vs magnitude
plt.figure(figsize=(8, 6))
plt.plot(sorted_brightness, log_counts, marker='o', linestyle='none')
plt.
plt.xlabel('Magnitude')
plt.ylabel('log(N)')
plt.title('Cumulative Number of Sources vs. Magnitude')
plt.grid(True)

plt.show()

# %%
