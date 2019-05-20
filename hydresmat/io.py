"""File IO routines for HydResMat."""


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


import dolfin as dol
from hydresmat.common import COMM_WORLD, isOldDolfin

__all__ = ["loadMeshdata", "loadSimdata3", "saveSimdata"]

def loadMeshdata(meshpath):
    """Load mesh, subdomains and boundaries from a given HDF5 file.

    Parameters
    ----------
    meshpath: str
        Path to mesh file (HDF5).

    Returns
    -------
    Mesh data from the "/mesh", "/subdomains" and "/boundaries" section
    respectively.
    """
    mesh = dol.Mesh()
    with dol.HDF5File(COMM_WORLD, meshpath, "r") as hdf:
        hdf.read(mesh, "/mesh", False)
        if isOldDolfin():
            subdomains = dol.MeshFunction("size_t", mesh)
        else:
            subdomains = dol.MeshFunction("size_t", mesh, dim=3)
        hdf.read(subdomains, "/subdomains")
        if isOldDolfin():
            boundaries = dol.MeshFunction("size_t", mesh)
        else:
            boundaries = dol.MeshFunction("size_t", mesh, dim=2)
        hdf.read(boundaries, "/boundaries")
    return mesh, subdomains, boundaries

def loadSimdata3(paths, mesh):
    """Load three simulation results from a given set of three paths. This is
    useful for loading all three simulation results of either translational or
    rotational motion at once. This process requires information on the
    simulation mesh.

    Parameters
    ----------
    paths
        Array of three paths to load.
    mesh
        Mesh data from "/mesh" section.

    Returns
    -------
    Two arrays of velocity and pressured fields respectively.
    """
    V = dol.FunctionSpace(mesh, "CG", 1)
    Q = dol.VectorFunctionSpace(mesh, "CG", 2)
    # Declaring the function space for each velocity and pressure from the
    # simulations
    us = [dol.Function(Q) for i in range(3)]
    ps = [dol.Function(V) for i in range(3)]

    # Reading the solution for all simulations
    for path, u, p in zip(paths, us, ps):
        with dol.HDF5File(COMM_WORLD, path, 'r') as fsim:
            fsim.read(u, "/velocity")
            fsim.read(p, "/pressure")
    return us,ps

#
def saveSimdata(savepath, u, p):
    """Save a single simulation result as a file.

    Parameters
    ----------
    savepath: str
        Path to file to save
    u
        Velocity field
    v
        Pressure field
    """
    with dol.HDF5File(COMM_WORLD, savepath, 'w') as fsim:
        fsim.write(u, "/velocity")
        fsim.write(p, "/pressure")
