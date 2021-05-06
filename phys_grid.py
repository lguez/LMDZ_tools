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
    import sys
    from mpl_toolkits import basemap
    import jumble
    import matplotlib
    from matplotlib import pyplot as plt, colors, cm, ticker

    if len(sys.argv) != 2: sys.exit("Required argument: path of input file (a "
                                    "limit or startphy NetCDF file)")
    matplotlib.interactive(True)
    cmap = cm.autumn

    # Whole world:
    b = basemap.Basemap(projection = "robin", lon_0 = 0)
    method = b.pcolormesh
    fig_list = [(b, method)]

    # Northern hemisphere:
    b = basemap.Basemap(projection="nplaea", boundinglat = 10, lon_0 = 0, 
                        round=True)
    method = b.pcolor
    fig_list.append((b, method))

    # Southern hemisphere:
    b = basemap.Basemap(projection="splaea", boundinglat = - 10, lon_0 = 0, 
                        round=True)
    method = b.pcolor
    fig_list.append((b, method))

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
            level_min = my_var.min()
            level_max = my_var.max()
            levels = ticker.MaxNLocator(nbins = 5).tick_values(level_min,
                                                               level_max)
            norm = colors.BoundaryNorm(levels, cmap.N)

            for b, method in fig_list:
                plt.figure()
                method(long_edge_mesh, lat_edge_mesh, my_var, latlon = True, 
                       cmap = cmap, norm = norm)
                b.colorbar()
                b.drawparallels(range(-90,91,45), labels=[1,0,0,0])
                b.drawmeridians(range(-180,181,60), labels=(0,0,0,1))
                b.drawcoastlines()
                plt.title(var_name)
                plt.suptitle(sys.argv[1])

    # No plt.show() since matplotlib.interactive(True)
