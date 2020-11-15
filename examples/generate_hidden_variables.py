"""
This scripts generates two examples of the files containing the hidden variables
  required to run the generate_edgelist_from_modelSD program.

Copyright (C) 2020  Antoine Allard (antoineallard.info)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import numpy as np
import pandas as pd
import scipy.stats as st


def to_fwf(df, fname, cols=None, width=15, prec=7, vnames=True):
    """Custom method 'to_fwf' for Pandas

    Parameters
    ----------
    fname : str
        The path to the new file in which the hidden parameters will be written.
    cols : list or array of strings, optional
        A list or an array containing the name of the columns (as strings) to
        be written in the file (if None all columns will be written).
    width : integer, optional
        Width of the columns (default is 15 characters).
    prec : decimal precision of the floats written into the file (default is 7).
    vnames : bool, optional
        Boolean indicating whether the first column contains the names of the
        vertices and should therefore be treated as strings (default is True).
    """

    if cols is None:
        cols = df.columns

    with open(fname, 'w') as ofile:

        fmt_str = '{:>' + '{}'.format(width) + '} '
        string_to_write = ''.join(fmt_str.format(i) for i in cols) + '\n'
        ofile.write(string_to_write.replace(' ', '#', 1))

        fmt_str = []
        nb_cols = len(cols)
        if vnames:
            fmt_str.append('%{}s '.format(width))
            nb_cols = nb_cols - 1

        fmt_str.extend(['%{}.{}f '.format(width, prec) for i in range(nb_cols)])
        fmt_str = ''.join(fmt_str)

        np.savetxt(ofile, df[cols].values, fmt=fmt_str)

# Adds the custom method to Pandas.
pd.DataFrame.to_fwf = to_fwf


if __name__ == "__main__":
    # Number of vertices.
    nb_vertices = 250

    # Generates the names of the vertices (v + number).
    vnames = np.array(['v' + str(i).zfill(5) for i in range(nb_vertices)])

    # Generates the kappas.
    #   distributed according to an exponential prob. density function.
    average_kappa = 15
    kappa = st.expon.rvs(scale=average_kappa, size=nb_vertices)

    # Generates the azimutal angular positions (for both S1 and S2)
    #   uniformly distributes in [0, 2pi)
    theta = 2 * np.pi * np.random.random_sample(nb_vertices)

    # Generates the polar angular positions (for S2).
    #   uniformly distributes in [0, pi)
    phi = np.arccos(2 * np.random.random_sample(nb_vertices) - 1)

    # Creates a dataframe with everything in it.
    df = pd.DataFrame(data={'vertex': vnames, 'kappa': kappa,
                            'theta': theta, 'phi': phi})

    # Writes an example of hidden variables files for the S1 model.
    df.to_fwf('graph01_expo_S1_hidvar.dat', ['vertex', 'kappa', 'theta'])

    # Writes an example of hidden variables files for the S2 model.
    df.to_fwf('graph01_expo_S2_hidvar.dat', ['vertex', 'kappa', 'theta', 'phi'])
