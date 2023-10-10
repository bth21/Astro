

# Specify the number of bins you want in your histogram
num_bins = 350


# Plot the histogram
plt.hist(filtered_data, bins=num_bins, density=False, alpha=0.75, color='blue')

# Add labels and a title
plt.xlabel('Pixel Value')
plt.ylabel('Frequency')
plt.title('Histogram of Pixel Values')


# Show the plot
plt.show()