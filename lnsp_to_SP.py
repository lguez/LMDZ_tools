#!/usr/bin/env python3

"""This script replaces the NetCDF variable lnsp by its exponential, SP.

Author: Lionel GUEZ"""

import netCDF4
import numpy as np
import sys

if len(sys.argv) != 2: sys.exit("Required argument: FILE")

with netCDF4.Dataset(sys.argv[1], "r+") as f:
    f.renameVariable("lnsp", "SP")
    SP = f.variables["SP"]
    SP[:] = np.exp(SP[:])
    SP.long_name = "surface pressure"
    SP.units = "Pa"
    SP.standard_name = "surface_air_pressure"
    for my_attr in ["num_GRIB", "level_desc", "dataset", "actual_range", 
                    "cell_methods"]:
        if my_attr in SP.ncattrs():
            SP.__delattr__(my_attr)
