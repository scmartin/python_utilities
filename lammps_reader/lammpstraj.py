import numpy as np
import argparse
import pandas as pd
import timeit

def traj_reader(fname):
    step = 0
    atoms = 0
    header = []
    first = True
    columntypes = {'q':'float','x':'float','y':'float','z':'float',
                   'element':'str','vx':'float','vy':'float',
                   'vz':'float','c_st[1]':'float',
                   'c_st[2]':'float','c_st[3]':'float','c_st[4]':'float',
                   'c_st[5]':'float','c_st[6]':'float','fx':'float',
                   'fy':'float','fz':'float'}
    with open(fname) as f:
        count = 1
        test1 = False
        test2 = False
        test3 = False
        snapshot = []
        for line in f:
            if line == "ITEM: TIMESTEP\n":
                step = int(f.readline())
                test1 = True
            elif line == "ITEM: NUMBER OF ATOMS\n":
                atoms = int(f.readline())
                test2 = True
            elif line.startswith("ITEM: ATOMS "):
                header = line.split()[3:]
                test3 = True
            elif (test1 and test2 and test3):
                snapshot.append(line.split())
                count += 1
                if (count > atoms):
                    count = 1
                    test1, test2, test3 = [False, False, False]
                    index = [np.repeat(step, atoms),
                            np.array(snapshot)[:,0].astype('int')]
                    snapshot = np.array(snapshot)[:,1:]
                    if first:
                        traj = pd.DataFrame(snapshot, index=index,
                                            columns=header)
                        first = False
                    else:
                        tstep = pd.DataFrame(snapshot, index=index, 
                                            columns=header)
                        traj = pd.concat((traj,tstep))
                    snapshot = []

    return traj.astype(columntypes)

def pandas_reader(fname, *args, columns=None):
    
    count = 0
    timesteps = []
    headers = {}
    frames = {}
    box = {}
    atoms = {}
    with open(fname) as f:  
        for line in f:
            if line == "ITEM: TIMESTEP\n":
                step = int(f.readline())
                timesteps.append(step)
                junk = f.readline()
                atoms[step] = int(f.readline())
                junk = f.readline()
                x = np.array(f.readline().split(), dtype=float)
                y = np.array(f.readline().split(), dtype=float)
                z = np.array(f.readline().split(), dtype=float)
                box[step] = np.array([x,y,z])
                header = f.readline().split()[2:]
                count += 9
                headers[step] = count
            else:
                count += 1
    # now that we know what to ignore, creating a dataframe is as easy as
    if not columns:
        columns = header
    for step in timesteps:
        frame = pd.read_csv(fname, sep=' ', skiprows=headers[step], 
                                   engine='c', nrows=atoms[step], 
                                   index_col=False, names=header, 
                                   header=None, usecols=columns)
        frames[step] = frame.set_index('id', drop=True)
    frames = pd.concat(frames)
    return frames,box


if __name__ == '__main__':
    pass
