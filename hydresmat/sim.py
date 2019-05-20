"""Simulation code for HydResMat."""

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
from dolfin import grad, div, dx
from math import hypot, fabs, log, pi, e, sqrt
import numpy as np

__all__ = ["runSimulation"]

def runSimulation(
    mesh, subdomains, boundaries,
    cube_surface_idxs, particle_surface_idx,
    particle_bc, kind):
    """Perform HydResMat simulations.

    Run a FEM simulation with generated mesh data and particle boundary
    conditions. Both rotational and translational motion can be simulated.

    Parameters
    ----------
    mesh
        Data from the "/mesh" section of the meshfile
    subdomains
        Data from the "/subdomains" section of the meshfile
    boundaries
        Data from the "/boundaries" section of the meshfile
    cube_surface_idxs
        Array of 6 indices for the outer surfaces of the cuboid shaped domain
    particle_surface_idx
        Index of the particle surface
    particle_bc
        Array of three (angular) velocities of the particle
    kind: {"rot", "trans"}
        Type of motion
    """
    dol.parameters['ghost_mode'] = 'shared_facet'
    krylov_method = "minres"
    preconditioner = "petsc_amg"

    # Defining the function space for the calculations
    P2 = dol.VectorElement("Lagrange", mesh.ufl_cell(), 2)
    P1 = dol.FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    TH = P2 * P1
    W = dol.FunctionSpace(mesh, TH)

    # Defining no-slip boundary conditions at the 6 outer surfaces of the
    # cuboid-shaped simulation domain
    noslip = dol.Constant((0.0, 0.0, 0.0))
    bcs = [dol.DirichletBC(W.sub(0), noslip, boundaries, surf_idx)
        for surf_idx in cube_surface_idxs]

    # Defining no-slip boundary conditions at the 6 outer surfaces of the
    # cuboid-shaped simulation domain
    if(kind == "rot"):
        # For a rotational motion of the particle, the boundary condition for
        # the velocity field at the particle surface is given by v = omega x r
        particle_boundary = dol.Expression(
            ("o_y*x[2]-o_z*x[1]", "o_z*x[0]-o_x*x[2]", "o_x*x[1]-o_y*x[0]"),
            degree=2,
            o_x = particle_bc[0], o_y = particle_bc[1], o_z = particle_bc[2])
    elif(kind == "trans"):
        # For a translational motion of the particle, the boundary condition for
        # the velocity field at the particle surface is simply the particle's
        # velocity
        particle_boundary = particle_bc
    else:
        print("Unknown kind of motion."
            "Please use kind = rot or kind = trans. \n")
        print(kind)
        exit()
    # Append last boundary condition
    bcs.append(dol.DirichletBC(
        W.sub(0), particle_boundary, boundaries, particle_surface_idx))

    # Defining the variational problem
    (u, p) = dol.TrialFunctions(W)
    (v, q) = dol.TestFunctions(W)
    f = dol.Constant((0.0, 0.0, 0.0))
    a = dol.inner(grad(u), grad(v))*dx + div(v) * p * dx + q * div(u) * dx
    L = dol.inner(f, v) * dx

    # Definition for use in constructing the preconditioner matrix
    b = dol.inner(grad(u), grad(v))*dx + p*q*dx

    # Assembling the main system
    A, bb = dol.assemble_system(a, L, bcs)

    # Assembling the preconditioner system
    P, btmp = dol.assemble_system(b, L, bcs)

    # Creating the Krylov solver and AMG preconditioner
    solver = dol.KrylovSolver(krylov_method, preconditioner)

    # Associating the operator A and preconditioner matrix P
    solver.set_operators(A, P)

    # Computing the solution
    U = dol.Function(W)
    solver.solve(U.vector(), bb)

    # Getting subfunctions
    (u, p) = U.split()

    return u,p
