import numpy as np
from tqdm import tqdm
from pyevtk.hl import imageToVTK
import skfmm
import csv
from scipy.ndimage import zoom
from scipy.ndimage import gaussian_filter


thickness = 1.5

output_vtk = False

factor = 5
sz2 = 128
sz = int(sz2*factor)
s = sz2/sz
resize_factor = (s,s,s)

sigma = thickness*factor/1.2

a = np.ones((sz,sz,sz))
SD2 = np.ones((sz2,sz2,sz2))

print("Calculating initial Values...")
for i in tqdm(range(sz)):
    for j in range(sz):
        for k in range(sz):
            if (sz/2-i)**2 + (sz/2-j)**2 + (sz/2-k)**2 <= (sz/8)**2:
                a[i][j][k] = 0

print("Calculating porosity value using Gaussian Filter...")
filtered_data = gaussian_filter(a, sigma=[sigma,sigma,sigma], mode='nearest', truncate=int(thickness*5*factor))
porosity = zoom(filtered_data, resize_factor, order=1)
print("Processing complete.")

if output_vtk:
    print("Saving Voxel Data to VTK file...")
    imageToVTK('./output/point_data',pointData={'porosity':porosity})
    print("Processing complete.")


print('Saving porosity Data to CSV file')

csv_file_path = "porosity.csv"

# Save 3D array to CSV file
with open(csv_file_path, mode='w', newline='') as file:
    writer = csv.writer(file)
    # Write header
    writer.writerow([sz2, sz2, sz2, "porosity"])  # Change column names as needed
    # Write data
    for ix in tqdm(range(sz2)):
        for iy in range(sz2):
            for iz in range(sz2):
                writer.writerow([ix, iy, iz, format(porosity[ix, iy, iz],'.5E')])
print(f"porosity data has been saved to {csv_file_path}.")

