import numpy as np
import ROOT
import glob
import re

N = 512
x0 = np.linspace(-25, 25, N)
y0 = np.linspace(-25, 25, N)
z0 = np.linspace(-100, 100, N)

files = glob.glob("./*.npy")
z_s = np.linspace(60., -100., 100, dtype=np.float64)

charges = []
z = []

import matplotlib.pyplot as plt

for file in files:
    data = np.load(file)**2
    match = re.search(r"[-+]?\d*\.\d+|\d+", file)

    z_focus = float(match.group())
    z.append(z_focus)
    print(f"Processing: {file}, z={z_focus}")

    z_min = np.argmin(np.abs(z0 - z_focus))
    z_max = np.argmin(np.abs(z0 - (z_focus + 50.)))

    # plt.vlines(z_min, 0, 512)
    # plt.vlines(z_max, 0, 512)
    # plt.imshow(data[256,:,:])
    # plt.show()    


    charges.append(np.sum(data[256,:,z_min:z_max]))

sorted_pairs = sorted(zip(z, charges))

# Unzip the sorted pairs
sorted_distances, sorted_values = zip(*sorted_pairs)

# Convert back to lists if you want
sorted_distances = np.array(sorted_distances, np.float32)
sorted_values = np.array(sorted_values, np.float32)