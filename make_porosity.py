import vtk
from tqdm import tqdm
import csv
import numpy as np
import math

def add_margin_to_bounds(bounds, margin):
    '''
    Expand the grid region in the zyx direction by margin [m] unit.
    '''
    return [bounds[0] - margin, bounds[1] + margin,
            bounds[2] - margin, bounds[3] + margin,
            bounds[4] - margin, bounds[5] + margin]

def ratio_margin_to_bounds(bounds, factor):
    '''
    Apply a factor [-] to expand the grid region in the zyx direction.
    (ex) factor = 1.5, original_bounds = [-1, 1, -1, 1, -1, 1]
        => bounds = [-1.5, 1.5, -1.5, 1.5, -1.5, 1.5]
    '''
    return [bounds[0] * factor, bounds[1] * factor,
            bounds[2] * factor, bounds[3] * factor,
            bounds[4] * factor, bounds[5] * factor]

# -------------------------------------------
filename = "your_file.stl"
factor = 2.0
# margin = 1.0 # [m]
cell_dims = [100, 100, 100]  # x, y, z
thickness = 2.0 # [px] = Delta/dx
# -------------------------------------------

# Read the stl file
print("Reading STL file...")
reader = vtk.vtkSTLReader()
reader.SetFileName(filename)
reader.Update()

# Extract Poly data
poly_data = reader.GetOutput()

# Create Mesh Grid

# Add margin to the original size of the STL file
original_bounds = poly_data.GetBounds()
expanded_bounds = ratio_margin_to_bounds(original_bounds, factor)

# if you use margin
# expanded_bounds = ratio_margin_to_bounds(original_bounds, margin)

# x_min:0 x_max:1, y_min:2, y_max:3, z_min:4, z_max:5
bounds = expanded_bounds

mesh_pitch = [(bounds[1] - bounds[0]) / (cell_dims[0]-1),
              (bounds[3] - bounds[2]) / (cell_dims[1]-1),
              (bounds[5] - bounds[4]) / (cell_dims[2]-1)]
mins = [bounds[0] - mesh_pitch[0]/2,
        bounds[2] - mesh_pitch[1]/2,
        bounds[4] - mesh_pitch[2]/2]

points = vtk.vtkPoints()
print("Creating Mesh Grid Points...")
for ix in tqdm(range(cell_dims[0] + 1)):
    for iy in range(cell_dims[1] + 1):
        for iz in range(cell_dims[2] + 1):
            x = ix * mesh_pitch[0] + mins[0]
            y = iy * mesh_pitch[1] + mins[1]
            z = iz * mesh_pitch[2] + mins[2]
            points.InsertNextPoint(x, y, z)

structured_base_mesh = vtk.vtkStructuredGrid()
structured_base_mesh.SetExtent(0, cell_dims[0], 0, cell_dims[1], 0, cell_dims[2])
structured_base_mesh.SetPoints(points)

# Convert structured grid data to unstructured grid data
append = vtk.vtkAppendFilter()
append.AddInputData(structured_base_mesh)
append.Update()
base_mesh = append.GetOutput()

# Find the coordinates of the center of each Voxel
cell_centers = vtk.vtkCellCenters()
cell_centers.SetInputData(base_mesh)
cell_centers.Update()

# Create Voxel mesh from STL
center_points = cell_centers.GetOutput().GetPoints()
cell_list = vtk.vtkIdList()
sdf = vtk.vtkImplicitPolyDataDistance()
sdf.SetInput(poly_data)
distance_sdf = vtk.vtkDoubleArray()
distance_sdf.SetName("sdf")
print("Calculating SDF Values...")
for idx in tqdm(range(center_points.GetNumberOfPoints())):
    current_center = center_points.GetPoint(idx)
    distance = sdf.FunctionValue(current_center)
    distance_sdf.InsertNextValue(distance)
    if distance <= 0:
        cell_list.InsertNextId(idx)

base_mesh.GetCellData().SetScalars(distance_sdf)  # Add SDF values

# Save the grid data with SDF values to a VTK file
print("Saving Voxel Data to VTK file...")
writer = vtk.vtkXMLDataSetWriter()
writer.SetFileName("all_voxel.vtu")
writer.SetInputData(base_mesh)
writer.Update()
print("Processing complete.")

# =======================================

# Assuming base_mesh is already defined

# Get the cell data array
cell_data_array = base_mesh.GetCellData().GetArray("sdf")

# Create a list to store cell data
dist_3d_array = np.zeros(cell_dims)

# Fill the 3D array with SDF values
cell_idx = 0
for ix in range(cell_dims[0]):
    for iy in range(cell_dims[1]):
        for iz in range(cell_dims[2]):
            dist_3d_array[ix, iy, iz] = 0.5*math.tanh(cell_data_array.GetValue(cell_idx)/(thickness*mesh_pitch[0]))+0.5
            cell_idx += 1

# Specify the CSV file path
csv_file_path = "porosity.csv"

# Save 3D array to CSV file
with open(csv_file_path, mode='w', newline='') as file:
    writer = csv.writer(file)
    # Write header
    writer.writerow([cell_dims[0],cell_dims[1],cell_dims[2], "porosity"])  # Change column names as needed
    # Write data
    for ix in range(dist_3d_array.shape[0]):
        for iy in range(dist_3d_array.shape[1]):
            for iz in range(dist_3d_array.shape[2]):
                writer.writerow([ix, iy, iz, dist_3d_array[ix, iy, iz]])

print(f"Cell data has been saved to {csv_file_path}.")
