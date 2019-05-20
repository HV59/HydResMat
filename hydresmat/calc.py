"""Calculation of quantities from FEM simulation data."""

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
from math import hypot, fabs, log, pi, e, sqrt
import numpy as np
import numpy.linalg as la
from hydresmat.common import COMM_WORLD

__all__ = ["calc_submatrices"]

def det3(A):
    """Returns determinant of 3x3 matrix. given by a row-major multidimensional
    list (necessary for types incompatible with numpy)

    Parameters
    ----------
    A
        3x3 row-major matrix
    """
    return A[0][0]*A[1][1]*A[2][2] + A[0][2]*A[1][0]*A[2][1] + \
           A[2][0]*A[0][1]*A[1][2] - A[0][2]*A[1][1]*A[2][0] - \
           A[0][0]*A[1][2]*A[2][1] - A[2][2]*A[0][1]*A[1][0]

def calc_submatrices(
        particle_bc, mesh, subdomains, boundaries, particle_surface_idx,
        velocities, pressures):
    r"""Calculate submatrices of the hydrodynamic resistance matrix.

    Calcultes submatrices of the hydrodynamic resistance matrix based on FEM
    simulation results by solving a linear system of equations. Returns either
    the transpose of the submatrix C (called D) and the submatrix Omega, when
    given data from a "rot" simulation or the submatrices K and C, when given
    data from a "trans" simulation.

    Parameters
    ----------
    particle_bc
        List of three boundary conditions (each represented as a list of three
        (angular) velocities) used to obtain the simulation results.
    mesh
        Data from the "/mesh" section of the meshfile
    subdomains
        Data from the "/subdomains" section of the meshfile
    boundaries
        Data from the "/boundaries" section of the meshfile
    particle_surface_idx
        Index of the particle surface
    velocities
        Velocity fields obtained by simulations
    pressures
        Pressure fields obtained by simulations
    """
    # Defining the function space for the calculations
    P2 = dol.VectorElement("Lagrange", mesh.ufl_cell(), 2)
    P1 = dol.FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    TH = P2 * P1
    W = dol.FunctionSpace(mesh, TH)

    # Shorthand for velocities and pressures
    us = velocities
    ps = pressures

    # Defining the function space for the derivative
    D = dol.FunctionSpace(mesh, "DG", 1)

    # Getting the spatial coordinates of the mesh
    r = dol.SpatialCoordinate(mesh)

    # Since in the simulations the pressure was defined with the wrong sign to
    # solve the system easier, the sign has to be changed after the simulations
    ps = [-p_ for p_ in ps]

    # Calculating the denominator for the solution of the system of linear
    # equations
    denominator = -la.det(particle_bc)

    # Solve the system of linear equations for the velocity and pressure fields
    # using Cramer's rule
    V = [[None]*3 for i in range(3)]
    p = [None]*3

    for i in range(3):
        for j in range(3):
            A = np.copy(particle_bc.T).tolist()
            A[j] = [u[i] for u in us]
            V[i][j] = -det3(A) / denominator
        A = np.copy(particle_bc.T).tolist()
        A[i] = ps
        p[i] = -det3(A) / denominator

    # Calculating the derivatives of the velocity field with respect to each
    # dimension
    d_V = [[[V[i][j].dx(k)
             for k in range(3)] for j in range(3)] for i in range(3)]

    # Calculating the triadic as defined in the book of Happel and Brenner
    P = [[None for j in range(3)] for i in range(3)]
    # Calculate upper triangular
    P[0][0] = [-p[i] + 2 * d_V[0][i][0] for i in range(3)] # P_11i
    P[1][1] = [-p[i] + 2 * d_V[1][i][1] for i in range(3)] # P_22i
    P[2][2] = [-p[i] + 2 * d_V[2][i][2] for i in range(3)] # P_33i
    P[0][1] = [d_V[1][i][0] + d_V[0][i][1] for i in range(3)] # P_12i
    P[0][2] = [d_V[2][i][0] + d_V[0][i][2] for i in range(3)] # P_13i
    P[1][2] = [d_V[2][i][1] + d_V[1][i][2] for i in range(3)] # P_23i
    # Use symmetry to fill lower triangular
    lower_tri = [(1,0),(2,0),(2,1)]
    for i,j in lower_tri:
        P[i][j] = P[j][i]

    # Getting the normals of the mesh at the particle surface
    # Since the normals are defined to be directed outwards the mesh domain
    # (which means here inwards the particle), their sign has to be changed
    cc = dol.FacetNormal(mesh)
    cc *= -1

    # Defining measures
    ds = dol.Measure("ds")(domain=mesh, subdomain_data=boundaries)

    # Defining the integral for the transpose of the submatrix C (calc_rot)
    # or the integral for the submatrix K (calc_trans)
    b = [[P[0][i][j] * cc[0] +
          P[1][i][j] * cc[1] +
          P[2][i][j] * cc[2] for j in range(3)] for i in range(3)]

    # Index pairs for cross product calculation
    # (i.e. cross = [a[i] * b[j] - a[j] * b[i] for i,j in cross_ind])
    cross_ind = [(1,2),(2,0),(0,1)]

    # Defining the integral for the submatrix Omega (calc_rot) or the integral
    # for the submatrix C (calc_trans)
    c = [[r[i] * b[j][k] - r[j] * b[i][k]
          for k in range(3)] for i,j in cross_ind]

    # Solving the integrals over the particle surface, which is given via the
    # number in ds. This number has to be declared in the creation of the .msh
    # file in Gmsh
    b = [[-dol.assemble(b[i][j] * ds(particle_surface_idx))
        for j in range(3)] for i in range(3)]
    c = [[-dol.assemble(c[i][j] * ds(particle_surface_idx))
        for j in range(3)] for i in range(3)]

    # "rot":   b is the transpose (called D) of the submatrix C and c is the
    #          submatrix Omega
    # "trans": b is the submatrix K and c is the submatrix C
    return b, c
