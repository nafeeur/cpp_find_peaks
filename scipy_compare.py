import numpy as np
from scipy.signal import find_peaks

# Define the signal
signal = [0, 2, 1, 3, 0, 1, 2, 1, 5, 7, 8, 1, 3, 0]

# Define the parameters
height = (1.5, np.inf)
distance = 1
prominence = (0.5, np.inf)
width = (0.0, np.inf)
rel_height = 0.5

# Find peaks
peaks, properties = find_peaks(
    signal,
    height=height,
    distance=distance,
    prominence=prominence,
    width=width,
    rel_height=rel_height
)

# Print the results
print("Peaks at indices:", peaks)
for i, peak in enumerate(peaks):
    print(f"Peak at index {peak} has height {signal[peak]}, "
          f"prominence {properties['prominences'][i]}, "
          f"width {properties['widths'][i]}")

