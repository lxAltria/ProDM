'''
Calculate contact map between two groups of atoms (residues)
based on a specified cut-off distance using multiple MD simulation
trajectories (plotting is not included).   
'''

import numpy
import os
import parmed as pmd
import numpy.linalg
import MDAnalysis
from MDAnalysis.analysis import distances
from numpy import *
from glob import glob
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


inpath = '/home/cc/dataset/amber_138_f000000100000'
outpath = '/home/cc/dataset/amber_138_f000000100000/analysis'

# topo="../amber_nowat.prmtop"
# traj="../md.nowat.mdcrd"
topo="../amber.prmtop"
traj="../md.1.mdcrd.nc"
top = pmd.load_file(topo)
u = MDAnalysis.Universe(top,traj)

sims = []

# store target simulation directories in a list
# simulation directory names read from sim_list.dat file
with open('sim_list.dat') as fp:
    for line in fp:
        data = line.split()
        sims.append(data[0])


os.chdir(inpath)

timestep = 46
n_frames = 0
start = 0

# group_1 = u.select_atoms('resid 315 317 323 324 325 and name CA')
# group_2 = u.select_atoms('resid 315 317 323 324 325 and name CA')

group_1 = u.select_atoms('protein')
group_2 = u.select_atoms('protein')

n1 = len(group_1)
n2 = len(group_2)


contact_sum = numpy.zeros((n1, n2))

max_distance = 10.0


for i in range(0, len(sims)):
    p = sims[i]
    print(p)
    os.chdir(p)
    

    topo="amber.prmtop"
    traj= "md.1.mdcrd.nc"
    top = pmd.load_file(topo)
    u = MDAnalysis.Universe(top,traj)
    
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

    os.chdir('../')


contact_ratio = contact_sum/n_frames

print(contact_ratio)
np.save('conmap_orig.npz', contact_ratio)
