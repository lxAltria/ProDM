import parmed as pmd
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import rms
import numpy
import os
import numpy.linalg
import MDAnalysis
from MDAnalysis.analysis import distances
from numpy import *
from glob import glob
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

filename = '../r-100.bin'
arr = np.fromfile(filename, dtype=np.float32).reshape(46, 3, 1584039)

traj = "../md.1.mdcrd.nc"
top = pmd.load_file("../amber.prmtop")
u0 = mda.Universe(top, traj)

arr = arr.transpose((0, 2, 1))

# Update positions
# for i in range(0, 46):
#     u0.trajectory[i].positions = arr[i]
#     ref.trajectory[i].positions = arr[i]
    
with mda.Writer('traj_modified.dcd', u0.atoms.n_atoms) as W:
    i = 0
    for ts in u0.trajectory:
        a = u0.atoms
        a.positions = arr[i]
        W.write(a)
        i += 1

print("Update completed")

u = mda.Universe(top, 'traj_modified.dcd')

# ------ calculate map ------ #


timestep = 46
n_frames = 0
start = 45

# group_1 = u.select_atoms('resid 315 317 323 324 325 and name CA')
# group_2 = u.select_atoms('resid 315 317 323 324 325 and name CA')

group_1 = u.select_atoms('protein')
group_2 = u.select_atoms('protein')

n1 = len(group_1)
n2 = len(group_2)


contact_sum = numpy.zeros((n1, n2))

max_distance = 10.0

    
# group_1 = u.select_atoms('resid 315 317 323 324 325 and name CA')
# group_2 = u.select_atoms('resid 315 317 323 324 325 and name CA')
group_1 = u.select_atoms('protein')
group_2 = u.select_atoms('protein')


for ts in u.trajectory[start::timestep]:
    gp1 = group_1.positions
    gp2 = group_2.positions
    ts_dist = distances.distance_array(gp1, gp2, box=u.dimensions)
    ts_dist[ts_dist < max_distance] = 1
    ts_dist[ts_dist > max_distance] = 0
    contact_sum = ts_dist + contact_sum
    n_frames+=1

contact_ratio = contact_sum/n_frames

print(contact_ratio)

np.save('r100.npz', contact_ratio)
