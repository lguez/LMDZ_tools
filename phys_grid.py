#!/usr/bin/env python3

"""This script plots data from a file "limit.nc" or "(re)startphy.nc".
It reshapes a NetCDF primary variable, as well as longitude and
latitude, from the "physics" grid to the "dynamics" grid.

Author: Lionel Guez

"""

if __name__ == "__main__":
    import xarray as xr
    import numpy as np
    import sys
    from mpl_toolkits import basemap
    import jumble
    import matplotlib
    from matplotlib import pyplot as plt, colors, cm, ticker

    if len(sys.argv) != 2: sys.exit("Required argument: path of input file (a "
                                    "limit or startphy NetCDF file)")
    matplotlib.interactive(True)

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
    cmap = cm.autumn

    with xr.open_dataset(sys.argv[1]) as f:
        klon = f.dims["points_physiques"]

        longitude = f.variables["longitude"][:]
        iim = np.argwhere((f.longitude[2:] == f.longitude[1]).values).squeeze()[0] \
            + 1
        print("iim = ", iim)
        longitude = f.longitude.values[1:iim + 2]
        longitude[- 1] = longitude[0] + 360

        jjm = (klon - 2) // iim + 1
        print("jjm = ", jjm)

        latitude = np.empty(jjm + 1)
        latitude[0] = f.variables["latitude"][0]
        latitude[1:] = f.variables["latitude"][1::iim]

        long_edge = jumble.edge(longitude)
        long_edge[0] = -180
        long_edge[- 1] = 180

        lat_edge = jumble.edge(latitude)
        lat_edge[0] = 90
        lat_edge[- 1] = - 90

        long_edge_mesh, lat_edge_mesh = np.meshgrid(long_edge, lat_edge)

        my_var = np.empty((jjm + 1, iim + 1))

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

            my_var_nc = my_var_nc.values
            my_var[0, :- 1]= my_var_nc[0]
            my_var[1:jjm, :- 1] = my_var_nc[1:klon - 1].reshape(jjm - 1, iim)
            my_var[jjm, :- 1] = my_var_nc[klon - 1]

            # Duplicated longitude:
            my_var[:, - 1] = my_var[:, 0]

            level_min = my_var.min()
            level_max = my_var.max()
            levels = ticker.MaxNLocator(nbins = 5).tick_values(level_min, level_max)
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
