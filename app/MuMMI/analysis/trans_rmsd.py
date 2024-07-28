import parmed as pmd
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import rms

filename = '../r-10.bin'
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
ref = mda.Universe(top, 'traj_modified.dcd')


# ------ calculate RMSD ------ #

# choose reference frame
ref_protein = ref.select_atoms('protein')
protein = u.select_atoms('protein')

# init RMSD analysis
R = rms.RMSD(protein, ref_protein)
R.run()

rmsd = R.results.rmsd.T  # transpose makes it easier for plotting

print(u.trajectory[0].positions[0:100])
# print("RMSD (Å): ", rmsd[2])
formatted_rmsd = " ".join(f"{value:.6f}" for value in rmsd[2])
print(f"RMSD (Å): [{formatted_rmsd}]")