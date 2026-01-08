import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import warnings
warnings.filterwarnings('ignore')

import MDAnalysis as mda
from MDAnalysis.tests.datafiles import GRO, TRR, PSF, DCD
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import distances, contacts, align
from MDAnalysis.analysis.distances import distance_array

import sys, os

u = mda.Universe('dna_0.pdb')
molecules_list="resname D*"
molecules = u.select_atoms(molecules_list)
print(len(u.trajectory))
frame_ndx = 0
for ts in u.trajectory:
    time = u.trajectory.time
    print(u.trajectory.frame)
    os.makedirs('./reference_from_pdb', exist_ok=True)
    molecules = u.select_atoms('resname D*')
    molecules.write('./reference_from_pdb/dna_0_'+str(frame_ndx)+'.gro')
    frame_ndx += 1
   

input_gro = 'out_0.gro'
input_trr = 'out_gro.trr'

scores = []
for frame_ndx in range(len(u.trajectory)):
    frame_name = './reference/dna_0'+'_'+str(frame_ndx)+'.gro'
    print(frame_name)
    ref    = mda.Universe(frame_name)
    ref_P  = ref.select_atoms("name P and resname D*")
    ref_P0 = ref_P.positions - ref_P.center_of_mass()

    u1 = mda.Universe(input_gro, input_trr)
    ts_ndx = 0
    ts_list = []
    for ts in u1.trajectory:
        tag   = u1.atoms
        tag_P = u1.select_atoms("name P")
        P_rmsd = rms.rmsd(ref_P.positions, tag_P.positions, superposition=True)
        tag_P0 = tag_P.positions - tag_P.center_of_mass()
        R, rmsd = align.rotation_matrix(tag_P0, ref_P0)

        tag.translate(-tag_P.center_of_mass())
        tag.rotate(R)
        tag.translate(ref_P.center_of_mass())

        #tag.write("tag_on_ref.pdb")

        tag_P1 = tag.select_atoms("name P")
        cutoffs = [2., 3., 4., 5.]

        print(tag_P1)
        dist_list = []
        for ref_atom, tag_atom in zip(ref_P, tag_P1):
            atom_p1 = ref_atom.position
            atom_p2 = tag_atom.position
            dist = ((atom_p1[0]-atom_p2[0])**2+(atom_p1[1]-atom_p2[1])**2+(atom_p1[2]-atom_p2[2])**2)**0.5
            dist_list.append(dist)

        dist_array = np.asarray(dist_list)

        gdts = []
        for cutoff in cutoffs:
            gdts.append((dist_array <= cutoff).sum()/float(len(dist_array)))
        #GDT = np.mean(gdts)
        GDT = gdts[1]
        print(frame_ndx, ts_ndx, ' '.join('{:0.4f}'.format(i) for i in gdts), '%.4f'%GDT)
        ts_list.append([ts_ndx, GDT])
        ts_ndx += 1
    scores.append(ts_list)
    scores_arr = np.array(scores)


plt.ion()
#colors = [color for color in mcolors.TABLEAU_COLORS]
colors = list(mcolors.CSS4_COLORS.values())

try:
    _gdt_log_fh = open('calculated_GDT.txt', 'w')
except Exception as _e:
    print(f"[WARN] Could not open calculated_GDT.txt for writing: {_e}")
    _gdt_log_fh = None
    
fig, ax = plt.subplots()
for ndx in range(len(scores)):
    label_txt = 'pdb_'+str(ndx)
    ts_frame = scores_arr[ndx][:,0]
    ts_score = scores_arr[ndx][:,1]
    #plt.scatter(ts_frame, ts_score, label=label_txt, c=colors[ndx], marker='o')
    ax.plot(ts_frame, ts_score, label=label_txt, c=colors[ndx])
    max_index = np.argmax(ts_score)
    max_score = np.max(ts_score)
    print('native: ', ndx, ' streamline: ', max_index, ' max score: %6.3f'%max_score)
    
    if _gdt_log_fh is not None:
        _gdt_log_fh.write(f"native:  {ndx}  streamline:  {max_index}  max score:  {max_score:6.3f}\n")

# Close the log file once after the loop
if _gdt_log_fh is not None:
    _gdt_log_fh.flush()
    _gdt_log_fh.close()
    
ax.set_ylabel('# GDT Score',fontsize=15)
ax.set_xlabel('AlphaFold3 Frame',fontsize=15)

ax.legend()

plt.xticks([-1, 0, 1, 2, 3, 4, 5])
plt.show()
input('enter any key to close ...')

plt.savefig('dna_0_scores.eps',format='eps',dpi=1200)
plt.savefig('dna_0_scores.pdf')
plt.savefig('dna_0_scores.png')

plt.close()














