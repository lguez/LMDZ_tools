#!/usr/bin/env python3

"""This script plots data from a file "limit.nc" or "(re)startphy.nc".
It reshapes a NetCDF primary variable, as well as longitude and
latitude, from the "physics" grid to the "dynamics" grid.

Author: Lionel Guez

"""

import xarray as xr
import numpy as np


def get_lon_lat(my_dataset):
    klon = my_dataset.dims["points_physiques"]

    # Longitude:
    iim = (
        np.argwhere(
            (my_dataset.longitude[2:] == my_dataset.longitude[1]).values
        ).squeeze()[0]
        + 1
    )
    print("iim = ", iim)
    longitude = my_dataset.longitude.values[1 : iim + 2].copy()
    longitude[-1] = longitude[0] + 360

    # Latitude:
    jjm = (klon - 2) // iim + 1
    print("jjm = ", jjm)
    latitude = np.empty(jjm + 1)
    latitude[0] = my_dataset.variables["latitude"][0]
    latitude[1:] = my_dataset.variables["latitude"][1::iim]

    return longitude, latitude


def gr_fi_dyn(pfi, longitude, latitude):
    i = pfi.dims.index("points_physiques")
    klon = pfi.points_physiques.size
    iim = longitude.size - 1
    jjm = latitude.size - 1

    v = np.empty(pfi.shape[:i] + (jjm + 1, iim + 1) + pfi.shape[i + 1 :])

    # North pole:
    v[(slice(None),) * i + (0, slice(None, -1))] = pfi.values[
        (slice(None),) * i + (0, np.newaxis)
    ]

    # Points other than poles:
    v[(slice(None),) * i + (slice(1, jjm), slice(None, -1))] = pfi.values[
        (slice(None),) * i + (slice(1, klon - 1),)
    ].reshape(pfi.shape[:i] + (jjm - 1, iim))

    # South pole:
    v[(slice(None),) * i + (jjm, slice(None, -1))] = pfi.values[
        (slice(None),) * i + (klon - 1, np.newaxis)
    ]

    # Duplicated longitude:
    v[(slice(None),) * (i + 1) + (-1,)] = v[(slice(None),) * (i + 1) + (0,)]

    coords = dict(
        {d: pfi[d] for d in pfi.dims[:i] + pfi.dims[i + 1 :]},
        latitude=latitude,
        longitude=longitude,
    )
    return xr.DataArray(
        v,
        dims=pfi.dims[:i] + ("latitude", "longitude") + pfi.dims[i + 1 :],
        coords=coords,
    )


if __name__ == "__main__":
    # Note: we do not simply plot with xarray plot because it does not
    # plot cells at the poles in projections, and it does not plot
    # cells at the limit of the colorbar.

    import argparse

    import cartopy.crs as ccrs
    import jumble
    import matplotlib
    from matplotlib import pyplot as plt, colors, cm, ticker

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help = "limit or startphy NetCDF file")
    parser.add_argument("--output_file", "-o",
                        help = "create a new NetCDF file")
    args = parser.parse_args()

    if args.output_file:
        input_ds = xr.open_dataset(args.input_file)
        output_ds = xr.Dataset()
        longitude, latitude = get_lon_lat(input_ds)

        for var_name in input_ds:
            if var_name not in {"longitude", "latitude"}:
                my_var = gr_fi_dyn(input_ds[var_name], longitude, latitude)
                output_ds[var_name] = my_var

        output_ds.to_netcdf(args.output_file)
    else:
        matplotlib.interactive(True)
        cmap = cm.autumn
        src_crs = ccrs.PlateCarree()

        # Whole world:
        extents = None
        fig_list = [(ccrs.Robinson(), extents)]

        # Northern hemisphere:
        extents = (-180, 180, 50, 90)
        fig_list.append(
            (ccrs.LambertAzimuthalEqualArea(central_latitude=90), extents)
        )

        # Southern hemisphere:
        extents = (-180, 180, -90, -50)
        fig_list.append(
            (ccrs.LambertAzimuthalEqualArea(central_latitude=-90), extents)
        )

        with xr.open_dataset(args.input_file) as my_dataset:
            longitude, latitude = get_lon_lat(my_dataset)

            # Longitude bounds:
            long_edge = jumble.edge(longitude)
            long_edge[0] = -180
            long_edge[-1] = 180

            # Latitude bounds:
            lat_edge = jumble.edge(latitude)
            lat_edge[0] = 90
            lat_edge[-1] = -90

            long_edge_mesh, lat_edge_mesh = np.meshgrid(long_edge, lat_edge)

            while True:
                var_name = input("Name of NetCDF primary variable? ")
                if len(var_name) == 0 or var_name.isspace():
                    break

                try:
                    pfi = my_dataset[var_name]
                except KeyError:
                    print("Not found")
                    print("Variables are:")
                    print(list(my_dataset.variables))
                    continue

                for dim in my_dataset[var_name].dims:
                    if dim != "points_physiques":
                        l = input(f"Subscript of {dim} (0-based)? ")
                        l = int(l)
                        pfi = pfi[{dim: l}]

                my_var = gr_fi_dyn(pfi, longitude, latitude)

                # Colorbar levels:
                level_min = my_var.values.min()
                level_max = my_var.values.max()
                levels = ticker.MaxNLocator(nbins=5).tick_values(
                    level_min, level_max
                )
                norm = colors.BoundaryNorm(levels, cmap.N)

                for projection, extents in fig_list:
                    plt.figure()
                    ax = plt.axes(projection=projection)
                    pcolormesh_return = ax.pcolormesh(
                        long_edge_mesh,
                        lat_edge_mesh,
                        my_var,
                        transform=src_crs,
                        cmap=cmap,
                        norm=norm,
                    )
                    if extents:
                        ax.set_extent(extents, crs=src_crs)
                    plt.colorbar(
                        pcolormesh_return, orientation="horizontal", shrink=0.5
                    )
                    ax.coastlines()
                    ax.gridlines(draw_labels=True)
                    ax.set_title(var_name)
                    plt.suptitle(args.input_file)
                    plt.subplots_adjust(top=0.85, bottom=0)

        # No plt.show() since matplotlib.interactive(True)
