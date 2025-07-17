import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import os
import sys

def plot_3d_mesh(points_file, triangles_file, output_folder=".", file_identifier_str="initial"):
    """
    Reads 3D point and triangle data from text files (VCL3D and Triangles)
    and generates a 3D plot.

    Args:
        points_file (str): Path to the VCL3D file containing 3D point coordinates (time, X, Y, Z).
        triangles_file (str): Path to the Triangles file containing triangle vertex indices (Node1 Node2 Node3, 1-indexed).
        output_folder (str): Folder to save the output image.
        file_identifier_str (str): Identifier for the output filename.
    """
    print(f"Reading points from: {points_file}")
    print(f"Reading triangles from: {triangles_file}")

    try:
        # Load points from VCL3D (time, X, Y, Z)
        # We want to skip the first column (time) and don't expect a '#' comment header.
        # The delimiter is a comma.
        # The data starts from the first line, no header to skip.
        points_data = np.loadtxt(points_file, delimiter=',', usecols=(1, 2, 3))
        
        # Ensure points_data is at least 2D (handle single point case)
        if points_data.ndim == 1:
            points_data = points_data.reshape(1, -1)
        
        # Check if the number of columns is correct (should be 3: X, Y, Z)
        if points_data.shape[1] != 3:
            print(f"Error: VCL3D file '{points_file}' does not have 3 valid coordinate columns (X, Y, Z) after skipping time. Found {points_data.shape[1]}.")
            sys.exit(1) # Exit if the VCL3D format is unexpected
        
        points = points_data # Rename for consistency with original code
        
        print(f"Points shape: {points.shape}")
        print(f"Last 5 points:\n{points[-5:]}")

        # Load triangles from 'Triangles' file
        # No header or comments in this file, directly load integers
        triangles = np.loadtxt(triangles_file, dtype=int)
        
        # Ensure triangles is at least 2D (handle single triangle case)
        if triangles.ndim == 1:
            triangles = triangles.reshape(1, -1)
        
        # Check if the number of columns is correct (should be 3: Node1, Node2, Node3)
        if triangles.shape[1] != 3:
            print(f"Error: Triangles file '{triangles_file}' does not have 3 valid vertex index columns. Found {triangles.shape[1]}.")
            sys.exit(1) # Exit if the Triangles format is unexpected

        print(f"Tri shape: {triangles.shape}")
        print(f"Last 5 tri:\n{triangles[-5:]}")

        # Convert Fortran 1-indexed to Python 0-indexed for array access
        triangles_0_indexed = triangles - 1

        print(f"Loaded {len(points)} points and {len(triangles)} triangles.")

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Create a list of vertices for each triangle for Poly3DCollection
        # Ensure indices are within bounds of the points array
        max_idx_in_triangles = np.max(triangles_0_indexed)
        if max_idx_in_triangles >= len(points):
            print(f"Error: Triangle index {max_idx_in_triangles} is out of bounds for points array of size {len(points)}.")
            print("This indicates a mismatch between the VCL3D and Triangles files, or an off-by-one error in Fortran.")
            sys.exit(1)

        verts = points[triangles_0_indexed]

        # Create a Poly3DCollection (collection of polygons)
        mesh = Poly3DCollection(verts, alpha=0.8, edgecolors='k', linewidths=0.5)
        mesh.set_facecolor('lightgreen')

        ax.add_collection3d(mesh)
        ax.scatter(points[:, 0], points[:, 1], points[:, 2], color='red', s=5)

        ax.set_xlabel('X (mm)')
        ax.set_ylabel('Y (mm)')
        ax.set_zlabel('Z (mm)')
        ax.set_title(f'3D Mesh Plot ({file_identifier_str})')
        ax.view_init(elev=30)

        # Set limits for 3D plot based on your previous settings
        ax.set_xlim([0.4,1.2])
        ax.set_ylim([0.1,0.7])
        ax.set_zlim([1.5,3])
  
        output_filename = os.path.join(output_folder, f"plot_mesh_{file_identifier_str}.png")
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        print(f"Plot saved as {output_filename}")
        plt.show()
        plt.close(fig)

    except FileNotFoundError as e:
        print(f"Error: Required file not found - {e}. Make sure Fortran program ran successfully and generated '{points_file}' and '{triangles_file}'.")
    except Exception as e:
        print(f"An unexpected error occurred during plotting: {e}")
        # print full traceback for debugging
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    # Define the input filenames
    # These are now VCL3D and Triangles
    points_file_to_use = 'VCL3D'
    triangles_file_to_use = 'Triangles'

    # Get arguments from command line
    # sys.argv[0] is the script name itself
    # sys.argv[1] would be the output folder
    # sys.argv[2] would be the step number
    output_folder = "." # Default to current directory
    file_identifier_str = "current_erosion_step" # Default identifier

    if len(sys.argv) < 3:
        print("Python mesh function called incorrectly.")
        print("Usage: python plot_mesh.py <output_folder> <file_identifier_string>")
        sys.exit(1)

    output_folder = sys.argv[1]
    file_identifier_str = sys.argv[2]

    # Run the plotting function with the specified files
    plot_3d_mesh(points_file_to_use, triangles_file_to_use, output_folder, file_identifier_str)

    print("\n Plotting script finished")
    print(f"Look for plots in '{output_folder}'.")