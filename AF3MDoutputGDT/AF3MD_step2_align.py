#edit the path and name of the gro, trr files
#reference and target frames are the detected best score combination from step 1. 
#reference: 0 (it is the first pdb entry) 
#target frame: the # of frame that has the best score compared to the reference

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

import MDAnalysis as mda
from MDAnalysis.tests.datafiles import GRO, TRR, PSF, DCD
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import distances, contacts, align
from MDAnalysis.analysis.distances import distance_array

import sys, os

cutoffs = [2., 3., 4., 5.]

#customize the file names
input_gro = '1.gro'
input_trr = '1.trr'

ref_ndx = input("please enter reference frame index: \n")
tag_ndx = input("please enter taget frame index: \n")
print(ref_ndx, tag_ndx)

#customize the path to the saved reference file
ref_file = './references/dna_0_'+str(ref_ndx)+'.gro'
ref = mda.Universe(ref_file)
ref_P = ref.select_atoms("name P and resname D*")

tag = mda.Universe(input_gro, input_trr)

ts_ndx = 0
for ts in tag.trajectory[:-1]:
    if ts.frame != int(tag_ndx):
        continue

    print(ts.frame)
    tag_P = tag.select_atoms("name P and resname D*")

    P_rmsd = rms.rmsd(ref_P.positions, tag_P.positions, superposition=True)
    tag_P0 = tag_P.positions - tag_P.center_of_mass()
    ref_P0 = ref_P.positions - ref_P.center_of_mass()
    R, rmsd = align.rotation_matrix(tag_P0, ref_P0)
    
    tag.atoms.translate(-tag_P.center_of_mass())
    tag.atoms.rotate(R)
    tag.atoms.translate(ref_P.center_of_mass())
#the aligned target frame to reference is saved as aligned.gro 
#file name can be changed here
    tag.atoms.write("aligned.gro")
    tag_P1 = tag.select_atoms("name P and resname D*")
    
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
    GDT = np.mean(gdts)
    print(ts_ndx, gdts, ' mean score: %.4f'%GDT)

    






