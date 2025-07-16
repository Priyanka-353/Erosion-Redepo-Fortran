import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import os
import sys

def plot_3d_mesh(points_file, triangles_file, output_folder=".", file_identifier_str="initial"):
    """
    Reads 3D point and triangle data from text files and generates a 3D plot.


    Args:
        points_file (str): Path to the file containing 3D point coordinates (X Y Z).
        triangles_file (str): Path to the file containing triangle vertex indices (Node1 Node2 Node3, 1-indexed).
        output_filename (str): Name of the output image file (e.g., "mesh_plot.png").
    """
    print(f"Reading points from: {points_file}")
    print(f"Reading triangles from: {triangles_file}")

    try:
        # Load points (X, Y, Z)
        # np.loadtxt will skip lines starting with '#'
        points = np.loadtxt(points_file, comments='#')
        if points.ndim == 1: # Handle case of single point
            points = points.reshape(1, -1)
        if points.shape[1] != 3:
            print(f"Warning: Points file '{points_file}' does not have 3 columns. Found {points.shape[1]}. Assuming X, Y, Z.")
            # Attempt to use first 3 columns if more exist, or pad if less
            points_temp = np.zeros((points.shape[0], 3))
            points_temp[:, :min(points.shape[1], 3)] = points[:, :min(points.shape[1], 3)]
            points = points_temp

        print(f"Points shape: {points.shape}") 
        print(f"Last 5 points:\n{points[-5:]}")


        # Load triangles (vertex indices) - remember Fortran is 1-indexed
        triangles = np.loadtxt(triangles_file, comments='#', dtype=int)
        if triangles.ndim == 1: # Handle case of single triangle
            triangles = triangles.reshape(1, -1)
        if triangles.shape[1] != 3:
            print(f"Warning: Triangles file '{triangles_file}' does not have 3 columns. Found {triangles.shape[1]}. Assuming Node1, Node2, Node3.")
            triangles_temp = np.zeros((triangles.shape[0], 3), dtype=int)
            triangles_temp[:, :min(triangles.shape[1], 3)] = triangles[:, :min(triangles.shape[1], 3)]
            triangles = triangles_temp

        print(f"Tri shape: {triangles.shape}") 
        print(f"Last 5 tri:\n{triangles[-5:]}")

        # Convert Fortran 1-indexed to Python 0-indexed for array access
        triangles_0_indexed = triangles - 1

        print(f"Loaded {len(points)} points and {len(triangles)} triangles.")

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Create a list of vertices for each triangle for Poly3DCollection
        # Each element in `verts` is a (3, 3) array representing the 3D coordinates
        # of the three vertices of a triangle.
        verts = points[triangles_0_indexed]

        # Create a Poly3DCollection (collection of polygons)
        mesh = Poly3DCollection(verts, alpha=0.8, edgecolors='k', linewidths=0.5)
        # You can set the face color for the mesh
        mesh.set_facecolor('lightgreen')

        ax.add_collection3d(mesh)
        ax.scatter(points[:, 0], points[:, 1], points[:, 2], color='red', s=5)

        # Optional: Plot individual points if desired
        # ax.scatter(points[:, 0], points[:, 1], points[:, 2], color='red', s=5)

        ax.set_xlabel('X (mm)')
        ax.set_ylabel('Y (mm)')
        ax.set_zlabel('Z (mm)')
        ax.set_title('3D Mesh Plot')
        ax.view_init(elev=30)

        ax.set_xlabel('X (mm)', fontsize=12, labelpad=10) # Increase label font size, add padding 
        ax.set_ylabel('Y (mm)', fontsize=12, labelpad=10) 
        ax.set_zlabel('Z (mm)', fontsize=12, labelpad=10)
        # Set equal aspect ratio
        #ax.set_box_aspect([np.ptp(points[:,0]), np.ptp(points[:,1]), np.ptp(points[:,2])]) # Aspect ratio based on data range
        ax.set_xlim([0.4,1.2])
        ax.set_ylim([0.1,0.7])
        ax.set_zlim([1.5,3])
    
        # Construct the output filename 
        # os.path.join handles different OS path separators (e.g., / vs \) 
        output_filename = os.path.join(output_folder, f"plot_mesh_{file_identifier_str}.png") 
        # Save the plot 
        plt.savefig(output_filename, dpi=300, bbox_inches='tight') 
        print(f"Plot saved as {output_filename}") 
        plt.show()
        plt.close(fig) 

    except FileNotFoundError as e:
        print(f"Error: Required file not found - {e}. Make sure Fortran program ran successfully.")
    except Exception as e:
        print(f"An error occurred during plotting: {e}")

if __name__ == "__main__":
    # Define the input and output filenames
    points_file = 'mesh_points_3d.txt'
    triangles_file = 'mesh_triangles_3d.txt'
    # output_image = 'initial_mesh_plot.png' # Name for the output plot image

    # Get arguments from command line
    # sys.argv[0] is the script name itself
    # sys.argv[1] would be the output folder
    # sys.argv[2] would be the step number
    output_folder = "." # Default to current directory
    step_number = "initial" # Default step identifier

    if len(sys.argv) < 3:
        print("Python mesh funciton called incorrectly")
        sys.exit(1)

    output_folder = sys.argv[1]
    file_identifier_str = sys.argv[2]

    # Run the plotting function with arguments
    plot_3d_mesh(points_file, triangles_file, output_folder, file_identifier_str)

    print("\n Plotting script finished")
    print(f"Look for plots in '{output_folder}'.")