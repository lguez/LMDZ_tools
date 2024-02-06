""" Author: Lionel Guez """

import math

import numpy as np

def find_idx_nearest_val(array, value):
    """Find index of nearest value in a one-dimensional monotonous array."""

    n = len(array)

    if array[0] <= array[-1]:
        sorted_array = array
        index = range(n)
    else:
        sorted_array = array[::- 1]
        index = range(n - 1, - 1, - 1)

    i = np.searchsorted(sorted_array, value)

    if i == n:
        return index[n - 1]
    elif i == 0:
        return index[0]
    else:
        if abs(value - sorted_array[i - 1]) < abs(value - sorted_array[i]):
            return index[i-1]
        else:
            return index[i]

def read_line_array(filename):
    """Read one real array per line from a text file.

    We assume that the first two token of each line are the name of
    the array and the equal sign. The arrays on different lines do not
    have to have the same length. The function returns a dictionary of
    Numpy arrays. The keys are the names of the arrays.
    """

    my_dict = {}

    with open(filename) as f:
        for line in f:
            splitted_line = line.split()
            value_list = [float(token) for token in splitted_line[2:]]
            my_dict[splitted_line[0]] = np.array(value_list)

    return my_dict

def read_line_with_header(filename):
    """Read one real array per line from a text file.

    We assume that the first token of each line is a quoted header,
    which can contain white space. The arrays on different lines do
    not have to have the same length. The function returns a
    dictionary of Numpy arrays. The keys are the headers.

    """

    import csv

    my_dict = {}

    with open(filename) as f:
        reader = csv.reader(f, delimiter = " ", skipinitialspace = True, 
                            quoting = csv.QUOTE_NONNUMERIC)

        for line in reader:
            if line[-1] == "":
                # (white space at the end of the input line has
                # produced an empty string list item)
                del line[-1]

            my_dict[line[0]] = np.array(line[1:])

    return my_dict

def edge(x):
    """Returns an array with elements halfway between input values, plus
    edge values.

    x should be a numpy array. Useful for pcolormesh if available
    coordinates are at z points.

    """

    return np.insert((x[:-1] + x[1:]) / 2, (0, x.size - 1),
                     ((3 * x[0] - x[1]) / 2, (3 * x[- 1] - x[- 2]) / 2))

def lsR(o):
    try:
        for k, v in o.items():
           print(k)
           lsR(v)
    except AttributeError:
        print(type(o))
        try:
            print(o.dtype, o.shape)
            if not np.issctype(o.dtype):
                print("First element:")
                lsR(o.flat[0])
            else:
                print()
        except AttributeError:
            print()

def cartesian(longitude, latitude):
    """Compute cartesian coordinates on the unit sphere."""
    x = math.cos(latitude) * math.cos(longitude)
    y = math.cos(latitude) * math.sin(longitude)
    z = math.sin(latitude)
    return x, y, z
