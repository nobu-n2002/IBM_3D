import numpy as np
from tqdm import tqdm
from pyevtk.hl import imageToVTK
import skfmm
import csv
from scipy.ndimage import zoom


thickness = 1.5

output_vtk = False

factor = 5
sz2 = 128
sz = int(sz2*factor)
s = sz2/sz
resize_factor = (s,s,s)

a = np.ones((sz,sz,sz))
SD2 = np.ones((sz2,sz2,sz2))

print("Calculating initial Values...")
for i in tqdm(range(sz)):
    for j in range(sz):
        for k in range(sz):
            if (sz/2-i)**2 + (sz/2-j)**2 + (sz/2-k)**2 <= (sz/8)**2:
                a[i][j][k] = 0

print("Calculating SDF using FMM...")
PHI = np.where(a, 0, -1) + 0.5
SD = skfmm.distance(PHI, dx = 1)
print("Processing complete.")

resize_factor = (s,s,s)

print("Resizing...")
SD2 = zoom(SD, resize_factor, order=1)*s
print("Processing complete.")

print("Calculating porosity value...")
porosity = 0.5*np.tanh(SD2/thickness)+0.5
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

# imageToVTK('/home/nobu/workspace/Labo/work/3d/test',pointData={'sphere':a})