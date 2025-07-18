
-------------------------RUN INFO----------------------
To run (assuming Intel Fortran HPC is installed):
Open terminal and run the following.
docker run -it \
  -v "/Users/pbalaji/Desktop/Erosion and Redeposition Fortran Code/ErosionRedepoCode":/workspace \
  -w /workspace \
  intel/oneapi-hpckit:2022.2-devel-ubuntu20.04 \
  /bin/bash

source /opt/intel/oneapi/setvars.sh

apt install -y python3 python3-pip
pip3 install matplotlib numpy

which ifort

make -f makefile
make clean
make --debug=v

./Redepo_v3

-------------------------Generated Files----------------------
cell_bndry_values contains the on_bndry parameter for each cell to know which cells are on the surface and which lie on a boundary 
IonInfo contains informatin about the ion  position, direction, number of reflections from Erosion
triangle_intersection_log contains the output of the triangulation function. Rows show the ion number and locations of 1s indicate which cell the ion hit.  
NodeData_Cumulative stores data about the cells of the triangulation include coordinates and mass loss at each erosion timestep

-------------------------DEBUGGER (GDB) Obsolete----------------------
apt update
apt install gdb -y

gdb ./Redepo_v3
break Redepo_v3.f90:165  # Break at a specific line
run

print vcl3d(:,226) # Inspect variables
next     # Step over
step     # Step into
continue # Continue to next breakpoint

quit 