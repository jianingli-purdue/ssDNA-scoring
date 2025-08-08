#edit the path and name of the gro, trr files
#dna_0 is the pdb file downloaded from PDB
#dna_0_# is the saved specific structure from dna_0. One entry from PDB could have multiple entries in one file. 
#Save each sructure from dna_0 separately as dna_0_0, dna_0_1, etc. 

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

import sys

#dna_sequence   = '6j2w'
#ligand_resname = 'B2U'

#output file
gdt_score_log = []

u = mda.Universe('dna_0.pdb')
molecules_list="resname D*"
molecules = u.select_atoms(molecules_list)
print(len(u.trajectory))
frame_ndx = 0
for ts in u.trajectory:
    time = u.trajectory.time
    print(u.trajectory.frame)
#change the file path    
    molecules.write('./references/dna_0'+'_'+str(frame_ndx)+'.gro')
    frame_ndx += 1
   
#change the file name to the correct name
input_gro = '1.gro'
input_trr = '1.trr'

scores = []
for frame_ndx in range(len(u.trajectory)):
#change the file path
    frame_name = './references/dna_0'+'_'+str(frame_ndx)+'.gro'
    print(frame_name)
    ref    = mda.Universe(frame_name)
    ref_P  = ref.select_atoms("name P and resname D*")
    ref_P0 = ref_P.positions - ref_P.center_of_mass()

    u1 = mda.Universe(input_gro, input_trr)
    ts_ndx = 0
    ts_list = []
    for ts in u1.trajectory[:-1]:
        tag   = u1.atoms
        
        #Yutong removed tag_P = u1.select_atoms("name P") and changed to:
        tag_P = u1.select_atoms("name P and resname D*")
        
        P_rmsd = rms.rmsd(ref_P.positions, tag_P.positions, superposition=True)
        tag_P0 = tag_P.positions - tag_P.center_of_mass()
        R, rmsd = align.rotation_matrix(tag_P0, ref_P0)

        tag.translate(-tag_P.center_of_mass())
        tag.rotate(R)
        tag.translate(ref_P.center_of_mass())

        #tag.write("tag_on_ref.pdb")

        #Yutong removed tag_P1 = tag.select_atoms("name P") and changed to:
        tag_P1 = tag.select_atoms("name P and resname D*")
        cutoffs = [2., 3., 4., 5.]

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
        #print(frame_ndx, ts_ndx, ' '.join('{:0.4f}'.format(i) for i in gdts), '%.4f'%GDT)
        ts_list.append([ts_ndx, GDT])
        ts_ndx += 1
    scores.append(ts_list)
    scores_arr = np.array(scores)


plt.ion()
colors = plt.cm.tab20.colors

fig, ax = plt.subplots()
for ndx in range(len(scores)):
    label_txt = 'pdb_'+str(ndx)
    ts_frame = scores_arr[ndx][:,0]
    ts_score = scores_arr[ndx][:,1]
    ax.plot(ts_frame, ts_score, label=label_txt, c=colors[ndx])
    max_index = np.argmax(ts_score)
    max_score = np.max(ts_score)
    
    #print('native: ', ndx, ' streamline: ', max_index, ' max score: %6.3f'%max_score)
    #native: the specific entry from the main dna_0 file. 
    #streamline: the frame that has the highest GDT score with the native. 
    log_line = f'native:  {ndx}  streamline:  {max_index}  max score:  {max_score:6.3f}'
    print(log_line)
    gdt_score_log.append(log_line)

ax.set_ylabel('# GDT Score',fontsize=15)
ax.set_xlabel('Time (ns)',fontsize=15)

ax.legend()

plt.show()
input('enter any key to close ...')

#output score file
with open('GDT_scores1.txt', 'w') as f:
    f.write('\n'.join(gdt_score_log) + '\n')
    
plt.savefig('dna_0_scores.eps',format='eps',dpi=1200)
plt.savefig('dna_0_scores.pdf')
plt.savefig('dna_0_scores.png')

plt.close()















