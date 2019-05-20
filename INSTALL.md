Installation instructions
=========================

Dependencies
------------
The following software packages are required to run this code:
* Python 3
* FEniCS (version 2017.2 or above)
* DOLFIN (this may be included in your FEniCS distribution)

Make sure FEniCS is built with support for the HDF5 file format. For building
the documentation, pdoc3 is required (can be installed via pip).

Installation procedure
----------------------
Make the folder "hydresmat" available to Python by copying it to a location in
the Python searchpath. You can print all locations in your current search path
by running the following code in your shell

```bash
python3 -c "import sys; print('\n'.join(sys.path))"
```

You can also add a folder to the search path manually by appending it to the
PYTHONPATH environment variable

```bash
export PYTHONPATH=/path/to/folder/containing/hydresmat:$PYTHONPATH
```

The installation can be tested by running

```bash
python3 -m hydresmat.convert2h5 -h
```

which should print a command help page.
