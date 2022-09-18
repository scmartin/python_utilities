import numpy as np
import freud
import lammps_reader as lammps
import argparse
from scipy.spatial import Delaunay, distance
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Analyze a LAMMPS snapshot to find surface atoms.")
parser.add_argument('-f','--file',type=str, required=True)
parser.add_argument('-c','--cutoff',type=float, default=1000, help='cutoff value of the dipole magnitude to be a surface atom')
parser.add_argument('-r','--radius',type=float, default=6.0, help="neighbor list radius")
parser.add_argument('-m','--mu',type=float, default=0.0, help='center of the gaussian weighting function')
parser.add_argument('-s','--sigma',type=float, default=2.5, help='standard deviation of gaussian weight function')

# Read a lammps trajectory
args = parser.parse_args()
mu = args.mu
sigma = args.sigma
traj = lammps.lammpstraj.pandas_reader(fname=args.file)

# Take the final frame from the lammps trajectory
final = traj[0].index.levels[0][-1]
box = traj[1][final]
frame = traj[0].loc[final].copy()
# Generate a Freud box from the lammps trajectory
box = freud.box.Box(box[0,1] - box[0,0],
                    box[1,1] - box[1,0],
                    box[2,1] - box[2,0])

# Generate neighborlists
nlist = freud.locality.AABBQuery(box, frame.loc[:,['x','y','z']])
query = nlist.query(nlist.points,dict(mode='ball',r_max=args.sphere))
newnlist = query.toNeighborList()

def weight(points, mu, sigma):
    """Generate gaussian weights for neighbors of a central point.
    
    Keyword arguments:
    points -- cartesian coordinates of neighbors. Coordinate origin is located at the central point.
    mu -- center of the gaussian weighting function. mu = 0 sets maximum weight at the central point.
    sigma -- width of the gaussian weighting function.
    """
    dist = np.sqrt(np.sum(points**2,axis=1))
    return np.tile(np.exp(-0.5*(dist-mu)**2/sigma**2), (3,1)).T

# mu, sigma = input("enter mu and sigma: ").split()
# mu, sigma = float(mu), float(sigma)

dipoles = np.zeros(nlist.points.shape[0])
zdipoles = np.zeros(nlist.points.shape[0])
print(nlist.box.Lx, nlist.box.Ly, nlist.box.Lz)
for i in range(nlist.points.shape[0]):
    points = nlist.points[newnlist.point_indices[newnlist.query_point_indices == i]] - nlist.points[i]
    xwrap = np.abs(points[:,0]) > nlist.box.Lx*0.5
    ywrap = np.abs(points[:,1]) > nlist.box.Ly*0.5
    points[xwrap,0] = (points[xwrap,0] - np.sign(points[xwrap,0])*nlist.box.Lx)
    points[ywrap,1] = (points[ywrap,1] - np.sign(points[ywrap,1])*nlist.box.Ly)
    dipole = np.sum(points*weight(points, mu, sigma), axis=0)
    zdipoles[i] = np.dot([0,0,1],dipole)
    dipoles[i] = np.dot(dipole,dipole)

frame.loc[slice(None),'dipole_magnitude'] = dipoles
frame.loc[slice(None),'zdipole'] = zdipoles
fig, ax = plt.subplots(2,1)
frame.hist(column='dipole_magnitude', bins=35, ax=ax[0])
frame.hist(column='zdipole', bins=25, ax=ax[1])
plt.show()

frame.to_csv('dipoles.xyz', sep=' ')
surface = frame.index.values[dipoles >= args.cutoff]

frame.loc[surface].to_csv('surffinal.lammpstrj', sep=' ')
