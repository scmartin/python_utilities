#!/usr/bin/env python

import pandas as pd
import numpy as np

def create_dframe(datafile):
    """
    Creates pandas dataframes from LAMMPS log file.
    To generate a pandas dataframes from the data in the LAMMPS log file, we need to pick out the lines with the data.
    pandas can read data from a .csv and we can tell it to ingnore a number of lines of header and footer.
    This function returns a list of dataframes. The [-1] element of the returned list is the thermodynamic output of
    final LAMMPS run command.

    E.g.

    datafile = 'log.lammps'
    data = create_dframe(datafile)
    plotter(data[-1],'Temp',datafile)

    would plot the temperature from the production run in log.lammps.
    """

    count = 0
    header_list = []
    footer_list = []
    dflist = []
    with open(datafile) as f:  # this opens a file as a file object called f, which we can then loop over
        for line in f:
            footertest = line.startswith('Loop time') or line.startswith('slurmstep') or line.startswith('srun: Job step')
            if footertest:
                footer_list.append(count)
            elif line.startswith('Step') or line.startswith('Time'):
                header_list.append(count)
            count += 1
    if len(footer_list) < len(header_list):
        footer_list.append(count)
    footer_list = np.array([count for i in range(len(footer_list))]) - np.array(footer_list)
    # now that we know what to ignore, creating a dataframe is as easy as
    dflist.append(pd.read_csv(datafile, delim_whitespace=True, skiprows=header_list[0], engine='python', skipfooter=footer_list[0], index_col=0))
    for i in range(len(header_list)-1):
        dflist.append(pd.read_csv(datafile, delim_whitespace=True, skiprows=header_list[i+1], engine='python', skipfooter=footer_list[i+1], index_col=0))
    return dflist

def plotter(dframe, column, ax):
    """ Plot the <column> from <dframe>.
    <dframe> needs to be a pandas dataframe. Use lammps_log.create_dframe(datafile)
    where datafile is the path to the LAMMPS log file to create a list of dataframes.
    See the help for create_dframe() to learn more.
    """
    import matplotlib.pyplot as plt

    try:
        dframe.plot(y=column, ax=ax, marker='None', linestyle='-')
    except TypeError:
        print("The data in this column is not numeric data")


def savestats(dframe, column, datafile):
    """Print the basic statistics of <column> from <dframe>."""
    try:
        print(dframe[column].describe())
        # dframe[column].describe().to_csv(datafile+".stats")
    except TypeError:
        print("The data in this column is not numeric data")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, type=str, help="LAMMPS log file path")
    parser.add_argument('-f', '--func', type=str, nargs='*', default=['plot', 'stats'], choices=['plot', 'stats'], help="List of functions you want to perform. Default is all of them")
    parser.add_argument('-c', '--column', type=str, nargs='*', default=[], help="Column headers that you want to analyze. Default is all of them")
    parser.add_argument('-d','--dframe', type=int, default=-1, help="Choosing which data frame to analyze. Default is the production run")

    args = parser.parse_args()
    datafile = args.input
    dframe = args.dframe
    column = args.column
    func = args.func

    data = create_dframe(datafile)[dframe]
    if not column:
        column = data.columns.values

    func_dict = {'plot':plotter,'stats':savestats}


    for i in column:
        for j in func:
            function = func_dict[j]
            function(data, i, datafile)
