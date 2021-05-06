#!/usr/bin/env python3

"""This script plots data from a file "limit.nc" or "(re)startphy.nc".
It reshapes a NetCDF primary variable, as well as longitude and
latitude, from the "physics" grid to the "dynamics" grid.

Author: Lionel Guez

"""

import xarray as xr
import numpy as np

def get_lon_lat(f):
    klon = f.dims["points_physiques"]

    # Longitude:
    iim = np.argwhere((f.longitude[2:] == f.longitude[1]).values).squeeze()[0] \
        + 1
    print("iim = ", iim)
    longitude = f.longitude.values[1:iim + 2].copy()
    longitude[- 1] = longitude[0] + 360

    # Latitude:
    jjm = (klon - 2) // iim + 1
    print("jjm = ", jjm)
    latitude = np.empty(jjm + 1)
    latitude[0] = f.variables["latitude"][0]
    latitude[1:] = f.variables["latitude"][1::iim]
    
    return longitude, latitude

def gr_fi_dyn(pfi, longitude, latitude):
        i = pfi.dims.index('points_physiques')
        klon = pfi.points_physiques.size
        iim = longitude.size - 1
        jjm = latitude.size - 1
        
        v = np.empty(pfi.shape[:i] + (jjm + 1, iim + 1) + pfi.shape[i + 1:])

        # North pole:
        v[(slice(None),) * i  + (0, slice(None, - 1))] \
            = pfi.values[(slice(None),) * i  + (0, np.newaxis)]

        # Points other than poles:
        v[(slice(None),) * i  + (slice(1, jjm), slice(None, - 1))] \
            = pfi.values[(slice(None),) * i  + (slice(1, klon - 1),)]\
                    .reshape(pfi.shape[:i] + (jjm - 1, iim))

        # South pole:
        v[(slice(None),) * i + (jjm, slice(None, - 1))] \
            = pfi.values[(slice(None),) * i  + (klon - 1, np.newaxis)]

        # Duplicated longitude:
        v[(slice(None),) * (i + 1) + (- 1,)] \
            = v[(slice(None),) * (i + 1) + (0,)]

        coords = dict({d: pfi[d] for d in pfi.dims[:i] + pfi.dims[i + 1:]},
                      latitude = latitude, longitude = longitude)
        return xr.DataArray(v,  dims = pfi.dims[:i] + ("latitude", "longitude")
                            + pfi.dims[i + 1:], coords = coords)

if __name__ == "__main__":
    import cartopy.crs as ccrs
    import sys
    import jumble
    import matplotlib
    from matplotlib import pyplot as plt, colors, cm, ticker

    if len(sys.argv) != 2: sys.exit("Required argument: path of input file (a "
                                    "limit or startphy NetCDF file)")
    matplotlib.interactive(True)
    cmap = cm.autumn
    src_crs = ccrs.PlateCarree()

    # Whole world:
    extents = None
    fig_list = [(ccrs.Robinson(), extents)]

    # Northern hemisphere:
    extents = (-180, 180, 10, 90)
    fig_list.append((ccrs.NorthPolarStereo(), extents))

    # Southern hemisphere:
    extents = (-180, 180, -90, -10)
    fig_list.append((ccrs.SouthPolarStereo(), extents))

    with xr.open_dataset(sys.argv[1]) as f:
        longitude, latitude = get_lon_lat(f)

        # Longitude bounds:
        long_edge = jumble.edge(longitude)
        long_edge[0] = -180
        long_edge[- 1] = 180

        # Latitude bounds:
        lat_edge = jumble.edge(latitude)
        lat_edge[0] = 90
        lat_edge[- 1] = - 90

        long_edge_mesh, lat_edge_mesh = np.meshgrid(long_edge, lat_edge)

        while True:
            var_name = input("Name of NetCDF primary variable? ")
            if len(var_name) == 0 or var_name.isspace(): break

            try:
                my_var_nc = f[var_name]
            except KeyError:
                print("Not found")
                print("Variables are:")
                print(list(f.variables))
                continue

            for dim in f[var_name].dims:
                if dim != 'points_physiques':
                    l = input(f"Subscript of {dim} (0-based)? ")
                    l = int(l)
                    my_var_nc = my_var_nc[{dim: l}]

            my_var = gr_fi_dyn(my_var_nc, longitude, latitude)

            # Colorbar levels:
            level_min = my_var.values.min()
            level_max = my_var.values.max()
            levels = ticker.MaxNLocator(nbins = 5).tick_values(level_min,
                                                               level_max)
            norm = colors.BoundaryNorm(levels, cmap.N)

            for projection, extents in fig_list:
                plt.figure()
                plt.title(var_name)
                ax = plt.axes(projection = projection)
                pcolormesh_return = ax.pcolormesh(long_edge_mesh, lat_edge_mesh,
                                                   my_var, transform = src_crs,
                                                   cmap = cmap, norm = norm)
                if extents: ax.set_extent(extents, crs = src_crs)
                plt.colorbar(pcolormesh_return, orientation = 'horizontal',
                             shrink = 0.5)
                ax.coastlines()
                ax.gridlines(draw_labels=True)
                plt.suptitle(sys.argv[1])

    # No plt.show() since matplotlib.interactive(True)
