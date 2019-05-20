# -- Constants for demo --

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


import numpy as np
from numpy import sqrt, e, pi, log
import numpy.linalg as la

# Particle boundary conditions used in the demo

# For "rot" simulations: Angular velocities
omegas = np.array([
    [sqrt(2),        e,              3.0                ],
    [sqrt(2),        e * pi,         log(2) * sqrt(0.2) ],
    [log(5) * 0.154, pi * sqrt(0.1), e * log(0.142)     ]
])
for omega in omegas:
    omega /= la.norm(omega)

# For "trans" simulations: Particle velocities
U_0s = np.array([
    [e,              pi,             log(3)             ],
    [sqrt(2),        e * pi,         log(2) * sqrt(0.2) ],
    [log(5) * 0.154, pi * sqrt(0.1), e * log(0.142)     ]
])
for U_0 in U_0s:
    U_0 /= la.norm(U_0)

# Surface indices
# (these are set in the .geo file)

cube_surface_idxs = [21, 23, 25, 27, 29, 31]
particle_surface_idx = 32

# Paths used in the demo

meshpath = "mesh/demo.h5"

savepaths_simrot = [
    "simulation/particleSimRot1.h5",
    "simulation/particleSimRot2.h5",
    "simulation/particleSimRot3.h5"
]
savepath_solution_rot = "solution/rot.txt"


savepaths_simtrans = [
    "simulation/particleSimTrans1.h5",
    "simulation/particleSimTrans2.h5",
    "simulation/particleSimTrans3.h5"
]
savepath_solution_trans = "solution/trans.txt"
