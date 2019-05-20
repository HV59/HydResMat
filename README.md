HydResMat - FEM-based Code for Calculating the Hydrodynamic Resistance Matrix
=============================================================================

This code is related to the following article:

*J. Voß and R. Wittkowski, Hydrodynamic resistance matrices of colloidal
particles with various shapes, arXiv:1811.01269 (2018)*

When publishing results obtained by using our code, please cite this article.


About this file
---------------

This Readme file describes the HydResMat code, which allows calculating the
hydrodynamic resistance matrix of an arbitrarily-shaped colloidal particle. The
software is written in python and uses the finite-element software package
FEniCS with version number "2017.2.0" or "2018.1.0". 

We suggest reading the article mentioned at the top of this file for background
information on hydrodynamic resistance matrices. This article contains also
useful equations that allow to calculate, e.g., the following:
1. The hydrodynamic drag force and torque acting on a particle that moves with
prescribed translational and angular velocities as well as vice versa the
translational and angular velocities of a particle for prescribed force and
torque acting on the particle
2. The diffusion tensor of a particle from its hydrodynamic resistance matrix
3. How the hydrodynamic resistance matrix of a particle changes when the
position of the point of reference is changed
4. How the hydrodynamic resistance matrix of a particle changes when its size
is changed

The HydResMat code is released together with a demo case that corresponds to a
cylindrical particle with two spherical caps, where one cap is concave and the
other one is convex.

Calculating the hydrodynamic resistance matrix of a colloidal particle with the
help of this code requires only a few steps that are described in the following.


Mesh generation
---------------

To use the HydResMat code, you first need a finite-element mesh that defines
the particle shape and the simulation domain. For the calculations in the
article mentioned at the top of this file, we used the mesh generator Gmsh in
version "3.0.6", but you should be able to use also other software for
generating the mesh.

You should create a finite-element mesh for a large cuboid-shaped domain
containing a small particle near the center of the cuboid. Please ensure that
only the domain outside of the particle is meshed. This can save a lot of
computation time. The demo case includes a suitable mesh created with Gmsh.
It is important to know the index numbers of the surfaces of the mesh, since
they are needed for the FEniCS simulations later. With these index numbers, you
can assign boundary conditions on the surfaces. They are also required to
calculate surface integrals that determine the forces and torques acting on
the particle, which are needed in the process of obtaining the hydrodynamic
resistance matrix from the results of the FEM simulations (see further below).
You can change the index numbers of the mesh, but then you have to change them
in the code for the FEniCS simulations accordingly.

We used the graphical interface to handle Gmsh. There you have to load the .geo
file that defines the geometry of the simulation domain that shall be meshed.
After that, you have to click "Mesh -> 3D" to mesh the domain in three
spatial dimensions.
When the meshing is done, you can save the mesh by clicking "File -> Export"
and subsequently choosing a path to save the mesh under "name.msh" in the
Version 2 ASCII Format. Please do not select the options "Save all", "Save
parametric coordinates" or "Save on file per partition".

If you want to use Gmsh in the terminal, you have to use the command

```bash
gmsh -3 -format msh2 -o name.msh input.geo
```

where -3 is the number of dimensions you want to mesh, msh2 the format expected
by the FEniCS tool dolfin-convert used in the next step, "name.msh" the output
file and "input.geo" your .geo file.


Conversion to HDF5
-----------------------

Afterwards, you have to convert the mesh file to a format FEniCS can read. This 
is done with the command

```bash
dolfin-convert name.msh name.xml
```

where "name.msh" is the mesh file you created with Gmsh and "name.xml" is the
output file. This command creates three different .xml files.
Be careful to use the same version of FEniCS to convert the mesh file and to
carry out the FEM simulations later, because the dolfin-convert command can
look differently in different versions of FEniCS, which can lead to unphysical
results in the FEM simulations.

If you want to simulate in parallel, which is recommended to speed up the
calculations, you need to further convert from the xml to the h5 format.
This can be done with a program included in our code package. Make sure to
follow the instructions in the INSTALL file first. Then run

```bash
python3 -m hydresmat.convert2h5 name.xml name.h5
```

FEM simulations
---------------

Next, you need to perform FEM simulations with FEniCS based on the Stokes
equation to calculate flow fields in a liquid filling the simulation domain
outside of the particle, where relative motion of the particle and the liquid
is assumed. To obtain the hydrodynamic resistance matrix of a particle, one
needs to simulate six times with different types of relative motion. The
HydResMat code performs three simulations with a purely translational motion
and further three simulations with a purely rotational motion. This enables
to solve a system of linear equations that determine the hydrodynamic
resistance matrix. To consider different types of motion, it is sufficient
to adapt the boundary conditions at the surface of the particle, since the
Stokes equation is not time-dependent. For the purely translational motion and
for the purely rotational motion, our code chooses three linearly independent
(angular) velocity vectors each.

To perform the FEM simulations, you can use the module `hydresmat.sim`.
You can specify both boundary conditions and the type of motion. The simulation
can be run in parallel by starting your script with mpirun, i.e.

```bash
mpirun -np p python3 simulation.py
```

where `p` is the number of parallel processes. An example script where all six
simulations necessary for calculating the hydrodynamic resistance matrix are
performed is given in the `demo/` folder.

Calculation of the hydrodynamic resistance matrix
-------------------------------------------------

A good overview about the calculation of the hydrodynamic resistance matrix
from flow fields given by a velocity field and a pressure field can be found in
the following book chapter:

*J. Happel and H. Brenner, Low Reynolds Number Hydrodynamics: With Special
Applications to Particulate Media, 2nd ed., Mechanics of Fluids and Transport
Processes, Vol. 1 (Kluwer Academic Publishers, Dordrecht, 1991), chapter 5.*

To calculate now the submatrices K, C and Omega of the hydrodynamic resistance
matrix H={{K,C^T},{C,Omega}}, which is a symmetric 2x2-dimensional block matrix,
you can use the module `hydresmat.calc` of our code package. Again, an example
script is given in the `demo/` folder.

Demo
----

To run the demo, switch to the `demo` directory and run `demo.sh`. Please note 
that when performing simulations for scientific purposes, a higher resolution of 
the mesh should be chosen.

For the default resolution of the demo, we needed roughly 27 GB RAM to run the 
FEM simulations serially with 20 minutes run time per FEM simulation using an 
Intel(R) Core(Tm) i7-7700K CPU @ 4.20GHz processor.
For the higher resolution uncommented in the file `demo.geo`, we needed roughly 
52 GB RAM to run the FEM simulations serially with 2 hours and 20 minutes run 
time per FEM simulation using an Intel(R) Xeon(R) Gold 6140 CPU @ 2.30GHz
processor.
