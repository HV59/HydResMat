"""Common useful helpers, e.g. for dealing with DOLFIN versions or MPI
logging"""

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

def isOldDolfin():
  """Returns whether DOLFIN is older than version 2018.1 or not."""
  # Somewhat ironically, the method of retrieving the version string is
  # dependant on the DOLFIN version
  try:
    versionstr = dol.dolfin_version()
  except AttributeError:
    versionstr = dol.__version__
  return int(versionstr.split(".")[0]) < 2018

COMM_WORLD = None
"""DOLFIN version independent handle to MPI_COMM_WORLD communicator."""
if isOldDolfin():
  COMM_WORLD = dol.mpi_comm_world()
else:
  COMM_WORLD = dol.MPI.comm_world

def print2(*args, **kwargs):
  print(*args, **kwargs)
try:
  from mpi4py import MPI
  if MPI.COMM_WORLD.size > 1:
    def printMPI(*arg, **kwarg):
      if MPI.COMM_WORLD.rank != 0:
        return
      if not "flush" in kwarg:
        kwarg["flush"] = True
      print(*arg, **kwarg)
    print2 = printMPI
except ImportError:
  pass

print2.__doc__ = """Print wrapper that will try to only print on MPI rank 0 if
MPI is enabled."""
