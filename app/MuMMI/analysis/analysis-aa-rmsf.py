''' 
Measure and plot root mean square fluctuations for each residue
in various parts of the protein using MD trajectories.
'''

import sys
import math
import os
import random
from os import path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg
import parmed as pmd
import MDAnalysis as mda
import MDAnalysis.transformations
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
from MDAnalysis.analysis.rms import RMSF


mainFolder="/home/cc/dataset/amber_138_f000000100000"
cdir=mainFolder+"/analysis"

sims_full = mainFolder #+"/simulation_directory" 

os.chdir(cdir)
os.getcwd()



name = "7mfd_braf_parts"

lstPatch = [f for f in os.listdir(sims_full) if f.startswith('min_7mfd_')]

num_sims = len(lstPatch)

nframes=4600


# test run rgyr for one simulation - not used
array_rmsf = np.empty([num_sims,738,9], dtype= float)  # the second axis is the number of residues BRAF all is the max with 583 residue but numbered up to 738
array_rmsf[:] = np.nan
array_rmsf_name = np.array(['BRAF-all', 'BRAF-RBD', 'BRAF-CRD', 'BRAF-loop1', 'BRAF-loop2', 'BRAF-KD', '14_3_3-S729', '14_3_3-S365'])

s = 0
for p in lstPatch:

    # topo = sims_full + "/" + p + "/amber_nowat.prmtop"
    # traj = sims_full + "/" + p + "/md.nowat.mdcrd"
    topo = sims_full + "/" + p + "/amber.prmtop"
    traj = sims_full + "/" + p + "/md.1.mdcrd.nc"
    print(f"Read i, top {topo}, trj {traj}")


    # if data alreaddy computed use that 
    if os.path.isfile(f"{name}_rmsf_run_{s}.npz"):
        load_data = np.load(f"{name}_rmsf_run_{s}.npz", allow_pickle=True)
        array_rmsf[s,:,:] = load_data['rmsf'][s,:,:]
    else:
        top = pmd.load_file(topo)
        cSim_mobile = mda.Universe(top,traj)
        cSim_ref = mda.Universe(top,traj)

        cSim_mobile.trajectory[0:nframes]
        atoms_cSim_braf_bb    = cSim_mobile.select_atoms("name CA and resid 0:582")
        atoms_cSim_braf_rbd   = cSim_mobile.select_atoms("name CA and resid 0:72")   # Residue 156 to 228
        atoms_cSim_braf_crd   = cSim_mobile.select_atoms("name CA and resid 81:119")   # Residue 237 to 275
        atoms_cSim_braf_loop1 = cSim_mobile.select_atoms("name CA and resid 126:203")   # Residue 282 to 359 
        atoms_cSim_braf_loop2 = cSim_mobile.select_atoms("name CA and resid 215:292")  # Residue 371 to 448 # note 8 residue loop 2 ext not included
        atoms_cSim_braf_kd    = cSim_mobile.select_atoms("name CA and resid 301:567") # Residue 457 to 723

        atoms_cSim_14_3_3_bb          = cSim_mobile.select_atoms("name CA and resid 586:1045")
        atoms_cSim_14_3_3_bb_S729_CR3 = cSim_mobile.select_atoms("name CA and resid 586:815") # Residue 1 to 230
        atoms_cSim_14_3_3_bb_S365_CR2 = cSim_mobile.select_atoms("name CA and resid 816:1045") # Residue 1 to 230 

        cSim_ref.trajectory[0] 
        atoms_ref_braf_bb    = cSim_ref.select_atoms("name CA and resid 0:582")
        atoms_ref_braf_rbd   = cSim_ref.select_atoms("name CA and resid 0:72")   # Residue 156 to 228
        atoms_ref_braf_crd   = cSim_ref.select_atoms("name CA and resid 81:119")   # Residue 237 to 275
        atoms_ref_braf_loop1 = cSim_ref.select_atoms("name CA and resid 126:203")   # Residue 282 to 359
        atoms_ref_braf_loop2 = cSim_ref.select_atoms("name CA and resid 215:292")  # Residue 371 to 448 # note 8 residue loop 2 ext not included
        atoms_ref_braf_kd    = cSim_ref.select_atoms("name CA and resid 301:567") # Residue 457 to 723

        atoms_ref_14_3_3_bb          = cSim_ref.select_atoms("name CA and resid 586:1045")
        atoms_ref_14_3_3_bb_S365_CR2 = cSim_ref.select_atoms("name CA and resid 586:815") # Residue 1 to 230
        atoms_ref_14_3_3_bb_S729_CR3 = cSim_ref.select_atoms("name CA and resid 816:1045") # Residue 1 to 230

        

        print(f" - {s} RMSF 1 BRAF full, loop1 and loop2")
        aligner = align.AlignTraj(cSim_mobile, cSim_ref, select="name CA and resid 0:582", in_memory=True).run() # align on BRAF full 
        rmsfer = RMSF(atoms_cSim_braf_bb, verbose=True).run()
        array_rmsf[s,atoms_cSim_braf_bb.resnums[0]:atoms_cSim_braf_bb.resnums[-1]+1,1] = rmsfer.rmsf
        rmsfer = RMSF(atoms_cSim_braf_loop1, verbose=True).run()
        array_rmsf[s,atoms_cSim_braf_loop1.resnums[0]:atoms_cSim_braf_loop1.resnums[-1]+1,4] = rmsfer.rmsf
        rmsfer = RMSF(atoms_cSim_braf_loop2, verbose=True).run()
        array_rmsf[s,atoms_cSim_braf_loop2.resnums[0]:atoms_cSim_braf_loop2.resnums[-1]+1,5] = rmsfer.rmsf
        
        print(f" - {s} RMSF 2 BRAF RBD")
        aligner = align.AlignTraj(cSim_mobile, cSim_ref, select="name CA and resid 0:72", in_memory=True).run() # align on BRAF RBD 
        rmsfer = RMSF(atoms_cSim_braf_rbd, verbose=True).run()
        array_rmsf[s,atoms_cSim_braf_rbd.resnums[0]:atoms_cSim_braf_rbd.resnums[-1]+1,2] = rmsfer.rmsf

        print(f" - {s} RMSF 3 BRAF CRD")
        aligner = align.AlignTraj(cSim_mobile, cSim_ref, select="name CA and resid 81:119", in_memory=True).run() # align on BRAF CRD 
        rmsfer = RMSF(atoms_cSim_braf_crd, verbose=True).run()
        array_rmsf[s,atoms_cSim_braf_crd.resnums[0]:atoms_cSim_braf_crd.resnums[-1]+1,3] = rmsfer.rmsf

        print(f" - {s} RMSF 6 BRAF KD")
        aligner = align.AlignTraj(cSim_mobile, cSim_ref, select="name CA and resid 301:567", in_memory=True).run() # align on BRAF KD 
        rmsfer = RMSF(atoms_cSim_braf_kd, verbose=True).run()
        array_rmsf[s,atoms_cSim_braf_kd.resnums[0]:atoms_cSim_braf_kd.resnums[-1]+1,6] = rmsfer.rmsf

        print(f" - {s} RMSF 7 14_3_3_S365_CR2")
        aligner = align.AlignTraj(cSim_mobile, cSim_ref, select="name CA and resid 586:815", in_memory=True).run() # align on 14_3_3_S365_CR2
        rmsfer = RMSF(atoms_cSim_14_3_3_bb_S365_CR2, verbose=True).run()
        #array_rmsf[s,atoms_cSim_14_3_3_bb_S365_CR2.resnums[0]-1:atoms_cSim_14_3_3_bb_S365_CR2.resnums[-1],7] = rmsfer.rmsf
        array_rmsf[s,0:230,7] = rmsfer.rmsf

        print(f" - {s} RMSF 8 14_3_3_S729_CR3")
        aligner = align.AlignTraj(cSim_mobile, cSim_ref, select="name CA and resid 816:1045", in_memory=True).run() # align on 14_3_3_S729_CR3
        rmsfer = RMSF(atoms_cSim_14_3_3_bb_S729_CR3, verbose=True).run()
        #array_rmsf[s,atoms_cSim_14_3_3_bb_S729_CR3.resnums[0]-1:atoms_cSim_14_3_3_bb_S729_CR3.resnums[-1],8] = rmsfer.rmsf
        array_rmsf[s,0:230,8] = rmsfer.rmsf

        np.savez(f"{name}_rmsf_run_{s}", rmsf=array_rmsf, rmsf_names=array_rmsf_name) # save in cwd 

    s=s+1

np.savez(f"{name}_rmsf_all", rmsf=array_rmsf, rmsf_names=array_rmsf_name) # save in cwd 


# plot combined RMSF
plt.figure(figsize=(15, 20))
fig, axs = plt.subplots(8, figsize=(12, 18))
colors = ['red', 'black', 'blue', 'green']
x = np.arange(0,738) 
for i in range(8):
    for s in range(num_sims):
        if i<6:
            axs[i].plot(x+156, array_rmsf[s,:,i+1], lw=1, color=colors[s]);
            axs[i].text(1.01, 0.5, array_rmsf_name[i], va="center", ha="left", size=12, transform=axs[i].transAxes)
            axs[i].set_xlim([-1, 739])
        else:
            axs[i].plot(x, array_rmsf[s,:,i+1], lw=1, color=colors[s]);
            axs[i].text(1.01, 0.5, array_rmsf_name[i], va="center", ha="left", size=12, transform=axs[i].transAxes)
            axs[i].set_xlim([-1, 739])

fig.savefig(f"{name}_rmsf_01_new.eps", bbox_inches='tight')
fig.savefig(f"{name}_rmsf_01_new.pdf", bbox_inches='tight')

      
# calculate avg and sd - not really meaningfull for RMSF but here anyway 
for i in range(8):
    print(f"RMSF - {array_rmsf_name[i]:<11s} avg {np.nanmean(array_rmsf[:,:,i+1]):5.2f} sd_avg {np.std(np.nanmean(array_rmsf[:,:,i+1], axis=1)):.2f} sd_all {np.nanstd(array_rmsf[:,:,i+1]):.2f} - avg over residue and all x4 repeats")
