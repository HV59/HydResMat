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

import os, os.path
import hydresmat
from hydresmat import print2

import democonst as demo

# Check if simulation directory exists and create it if necessary
savedirs = [os.path.dirname(fn) for fn in
    [demo.savepaths_simrot[0], demo.savepaths_simtrans[0]]]
for savedir in savedirs:
  if not os.path.exists(savedir):
      os.mkdir(savedir)

print2("Reading mesh data...", flush=True)
mesh, subdomains, boundaries = hydresmat.loadMeshdata(demo.meshpath)

print2("Running simulations...")
# "rot" simulations
for i in range(3):
  print2("  Simulation: Omega = ({:.7f}, {:.7f}, {:.7f})".format(
      *demo.omegas[i]))
  u, p = hydresmat.runSimulation(
      mesh, subdomains, boundaries, demo.cube_surface_idxs,
      demo.particle_surface_idx, demo.omegas[i], "rot")
  hydresmat.saveSimdata(demo.savepaths_simrot[i], u, p)

# "trans" simulations
for i in range(3):
  print2("  Simulation: U_0 = ({:.7f}, {:.7f}, {:.7f})".format(
      *demo.omegas[i]))
  u, p = hydresmat.runSimulation(
      mesh, subdomains, boundaries, demo.cube_surface_idxs,
      demo.particle_surface_idx, demo.U_0s[i], "trans")
  hydresmat.saveSimdata(demo.savepaths_simtrans[i], u, p)
