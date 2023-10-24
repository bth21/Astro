import numpy as np
from skimage import io, measure
from scipy.ndimage import label
from astropy.convolution import Gaussian2DKernel
import astropy.units as u
import photometry as pt
import blooming as bl


# Load your binary image where objects are labeled
#binary_image = io.imread("your_binary_image.png")

# Label objects in the binary image
#labeled_image, num_features = label(binary_image)

# Define a function to calculate FWHM
def calculate_fwhm(data):
    # Assuming data is your object
    # You may need to adapt this for your specific data
    kernel = Gaussian2DKernel(2)  # Use an appropriate kernel size
    smoothed_data = convolve(data, kernel)
    return 2.355 * (smoothed_data == np.max(smoothed_data))

# Define a function to classify objects
def classify_objects(labeled_image, min_area, max_solidity, min_fwhm):
    props = measure.regionprops(labeled_image)

    stars = []
    galaxies = []

    for prop in props:
        area = prop.area
        solidity = prop.solidity
        fwhm = calculate_fwhm(prop.image)

        if area >= min_area and solidity <= max_solidity and fwhm >= min_fwhm:
            stars.append(prop)
        else:
            galaxies.append(prop)

            

    return stars, galaxies

# Define classification criteria
min_area = 10
max_solidity = 0.95
min_fwhm = 2.0  # You can adjust this value

# Classify objects
stars, galaxies = classify_objects(bl.final_clean_data, min_area, max_solidity, min_fwhm)

# Print the results
print(f"Number of Stars: {len(stars)}")
print(f"Number of Galaxies: {len(galaxies)}")

# classification function that can append data to ascii file
def classify_objects_4ascii(labeled_image, min_area, max_solidity, min_fwhm):
    props = measure.regionprops(labeled_image)

    for prop in props:
        area = prop.area
        solidity = prop.solidity
        fwhm = calculate_fwhm(prop.image)

        if area >= min_area and solidity <= max_solidity and fwhm >= min_fwhm:
            classify = 'star'
        else:
            classify = 'galaxy'

        
    return classify


