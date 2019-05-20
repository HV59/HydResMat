#!/bin/bash

# -- Create mesh --
[ -d mesh ] || mkdir mesh
gmsh -3 -format msh2 -o mesh/demo.msh demo.geo

# -- Convert mesh --
dolfin-convert mesh/demo.msh mesh/demo.xml
python3 -m hydresmat.convert2h5 mesh/demo.xml mesh/demo.h5

# -- Run simulations --
#mpirun -np 2 python3 sim.py
python3 sim.py

# -- Calculate hydrodynamic resistance matrix --
python3 calc.py rot
python3 calc.py trans
