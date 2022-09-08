import numpy as np
import freud
import lammps
import argparse
from scipy.spatial import Delaunay, distance
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Analyze a LAMMPS snapshot to find surface atoms.")
parser.add_argument('-f','--file',type=str, required=True)
parser.add_argument('-c','--cutoff',type=float, default=1000, help='cutoff value of the dipole magnitude to be a surface atom')
parser.add_argument('-s','--sphere',type=float, default=6.0, help="sphere of neighbor list radius")
parser.add_argument('-w', '--weighted', type=bool, default=False,
                    help="Apply Gaussian weight to dipole. Default=False")

args = parser.parse_args()
traj = lammps.lammpstraj.pandas_reader(fname=args.file)

final = traj[0].index.levels[0][-1]
box = traj[1][final]
frame = traj[0].loc[final].copy()
box = freud.box.Box(box[0,1] - box[0,0],
                    box[1,1] - box[1,0],
                    box[2,1] - box[2,0])

# oxygen = frame[frame['type'] == 2].copy()
# silicon = frame[frame['type'] == 1].copy()
# 
# nlist = freud.locality.AABBQuery(box, frame.loc[:,['x','y','z']])
# onlist = freud.locality.AABBQuery(box, oxygen.loc[:,['x','y','z']])
# sinlist = freud.locality.AABBQuery(box, silicon.loc[:,['x','y','z']])
# 
# 
# RDF = freud.density.RDF(30,10)
# RDF.compute(system=nlist)
# osiRDF = freud.density.RDF(30,10)
# osiRDF.compute(system=sinlist, query_points=oxygen.loc[:,['x','y','z']].to_numpy())
# ooRDF = freud.density.RDF(30,10)
# ooRDF.compute(system=onlist)
# sisiRDF = freud.density.RDF(30,10)
# sisiRDF.compute(system=sinlist)
# 
# fig, ax = plt.subplots(figsize=(10,7))
# ax.plot(RDF.bin_centers, RDF.rdf, '-k', label="all")
# ax.plot(osiRDF.bin_centers, osiRDF.rdf, '-r', label="O-Si")
# ax.plot(ooRDF.bin_centers, ooRDF.rdf, '-b', label="O-O")
# ax.plot(sisiRDF.bin_centers, sisiRDF.rdf, '-m', label="Si-SI")
# ax.legend()
# plt.show()

def weight(points, mu, sigma):
    dist = np.sqrt(np.sum(points**2,axis=1))
    return np.tile(np.exp(-0.5*(dist-mu)**2/sigma**2), (3,1)).T

mu, sigma = input("enter mu and sigma: ").split()
weightParams = [float(mu), float(sigma)]
# simu, sisigma = input("enter Si mu and sigma: ").split()
# si_weightParams = [float(simu), float(sisigma)]
# paramDict = dict(zip([2,1],[o_weightParams, si_weightParams]))

dipoles = np.zeros(frame.shape[0])
zdipoles = np.zeros(frame.shape[0])

for i in range(frame.shape[0]):
    x = frame['x'].values
    y = frame['y'].values
    z = frame['z'].values
    atomtype = frame.loc[i+1,'type']
    # wrap the coordinates
    x = x - x[i]
    x[np.abs(x) > 0.5*box.Lx] -= np.sign(x[np.abs(x) > 0.5*box.Lx])*box.Lx
    y = y - y[i]
    y[np.abs(y) > 0.5*box.Ly] -= np.sign(y[np.abs(y) > 0.5*box.Ly])*box.Ly
    z = z - z[i]
    coords = np.vstack((x,y,z)).T
    neighs = coords#[frame['type'] != atomtype] 
    neighs = neighs[neighs[:,0]**2 + neighs[:,1]**2 + neighs[:,2]**2 < args.sphere**2,:]
    if args.weighted:
        dipole = np.sum(neighs*weight(neighs, *weightParams), axis=0)
    else:
        dipole = np.sum(neighs, axis=0)
    zdipoles[i] = np.dot([0,0,1],dipole)
    dipoles[i] = np.dot(dipole,dipole)

frame.loc[slice(None),'dipole'] = dipoles
frame.loc[slice(None),'zdipole'] = zdipoles
fig, ax = plt.subplots(2,1, figsize=(6,6))
frame[frame["type"] == 2].hist(column='dipole', bins=21, ax=ax[0],
                               color='r', label='O')
frame[frame["type"] == 1].hist(column='dipole', bins=21, ax=ax[0],
                               color='k', label="Si")
frame[frame["type"] == 2].hist(column='zdipole', bins=21, ax=ax[1], color='r')
frame[frame["type"] == 1].hist(column='zdipole', bins=21, ax=ax[1], color='k')
ax[0].legend()
plt.show()

f = open('dipoles.xyz', 'w')
f.write('576\n')
frame.to_csv(f, sep=' ')
f.close()

surface = frame.index.values[dipoles >= args.cutoff]

f = open('surface.xyz', 'w')
f.write(f'{surface.shape[0]}\n')
frame.loc[surface].to_csv('surface.xyz', sep=' ')
f.close()
