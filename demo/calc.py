""" Copyright (C) 2018-2019 Johannes Voss, Julian Jeggle, Raphael Wittkowski

    This file is part of HydResMat.

    HydResMat is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HydResMat is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with HydResMat. If not, see <http://www.gnu.org/licenses/>."""

import hydresmat
import argparse
import os, os.path
import democonst as demo

# Parse args
parser = argparse.ArgumentParser(
    description="Calculate submatrices from simulation results.")
parser.add_argument("kind", help="Kind of simulation to calculate the"
    "corresponding submatrices from. Must be either 'rot' or 'trans'.")

args = parser.parse_args()

# Load appropriate constants for this kind of simulation
if not args.kind.lower() in ["rot", "trans"]:
    print("Kind must be either 'rot' or 'trans'.")

if args.kind == "rot":
    savepath_simdata = demo.savepaths_simrot
    savepath_solution = demo.savepath_solution_rot
    labels = ["D", "Omega"]
    particle_bc = demo.omegas
else:
    savepath_simdata = demo.savepaths_simtrans
    savepath_solution = demo.savepath_solution_trans
    labels = ["K", "C"]
    particle_bc = demo.U_0s

# Check if solution directory exists and create it if necessary
savedir = os.path.dirname(savepath_solution)
if not os.path.exists(savedir):
    os.mkdir(savedir)

print("Reading mesh data...")
mesh, subdomains, boundaries = hydresmat.loadMeshdata(demo.meshpath)

print("Reading simulation data...")
velocities, pressures = hydresmat.loadSimdata3(savepath_simdata, mesh)

# Calculate submatrices
sm1, sm2 = hydresmat.calc_submatrices(
    particle_bc, mesh, subdomains, boundaries, demo.particle_surface_idx,
    velocities, pressures)

# Write solutions to file
with open(savepath_solution, "w") as f:
    for i in range(3):
        f.write("\t".join(
            ["{}_{}{} = {:.5f}".format(labels[0], i+1, j+1, sm1[i][j])
                for j in range(3)]))
        f.write("\n")
    for i in range(3):
        f.write("\t".join(
            ["{}_{}{} = {:.5f}".format(labels[1], i+1, j+1, sm2[i][j])
                for j in range(3)]))
        f.write("\n")
