from schrodinger import structure
from schrodinger.structutils import rmsd, measure

from os import listdir
from os.path import isfile, join

import sys


def get_backbone_atoms(st):
    backbone_atoms=[]
    resnum = 0
    for res in st.residue:
        res_atoms = [None, None, None, None, None, None]   # P, O5', C5', C4', C3', O3'
        if res.chain=="A" and res.pdbres.strip() in ['DA', 'DC', 'DG', 'DT']:
            
  
            for at in res.atom:
                
                if at.pdbname.strip()=="P":
                    res_atoms[0] = at
                    #print(at)
                elif at.pdbname.strip()=="O5'":
                    res_atoms[1] = at
                elif at.pdbname.strip()=="C5'":
                    res_atoms[2] = at
                elif at.pdbname.strip()=="C4'":
                    res_atoms[3] = at  
                elif at.pdbname.strip()=="C3'":
                    res_atoms[4] = at
                elif at.pdbname.strip()=="O3'":
                    res_atoms[5] = at
            backbone_atoms.append(res_atoms)
            resnum = resnum + 1

    
    print("#ires   alpha    beta     gamma     delta    epsilon   zeta")
    for ires in range(resnum):
        if ires>0 and ires<resnum-1:
            alpha=measure.measure_dihedral_angle(backbone_atoms[ires-1][5], backbone_atoms[ires][0], backbone_atoms[ires][1], backbone_atoms[ires][2])
            beta=measure.measure_dihedral_angle(backbone_atoms[ires][0], backbone_atoms[ires][1], backbone_atoms[ires][2], backbone_atoms[ires][3])
            gamma=measure.measure_dihedral_angle(backbone_atoms[ires][1], backbone_atoms[ires][2], backbone_atoms[ires][3], backbone_atoms[ires][4])
            delta=measure.measure_dihedral_angle(backbone_atoms[ires][2], backbone_atoms[ires][3], backbone_atoms[ires][4], backbone_atoms[ires][5])
            epsilon=measure.measure_dihedral_angle(backbone_atoms[ires][3], backbone_atoms[ires][4], backbone_atoms[ires][5], backbone_atoms[ires+1][0])
            zeta = measure.measure_dihedral_angle(backbone_atoms[ires][4], backbone_atoms[ires][5], backbone_atoms[ires+1][0], backbone_atoms[ires+1][1])
            print("%4d %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f"%(ires, alpha, beta, gamma, delta, epsilon, zeta))
        elif ires == 0:
            gamma=measure.measure_dihedral_angle(backbone_atoms[ires][1], backbone_atoms[ires][2], backbone_atoms[ires][3], backbone_atoms[ires][4])
            delta=measure.measure_dihedral_angle(backbone_atoms[ires][2], backbone_atoms[ires][3], backbone_atoms[ires][4], backbone_atoms[ires][5])
            epsilon=measure.measure_dihedral_angle(backbone_atoms[ires][3], backbone_atoms[ires][4], backbone_atoms[ires][5], backbone_atoms[ires+1][0])
            zeta = measure.measure_dihedral_angle(backbone_atoms[ires][4], backbone_atoms[ires][5], backbone_atoms[ires+1][0], backbone_atoms[ires+1][1])
            print("%4d    NA       NA    %8.1f %8.1f %8.1f %8.1f"%(ires, gamma, delta, epsilon, zeta))
        elif ires==resnum-1:
            alpha=measure.measure_dihedral_angle(backbone_atoms[ires-1][5], backbone_atoms[ires][0], backbone_atoms[ires][1], backbone_atoms[ires][2])
            beta=measure.measure_dihedral_angle(backbone_atoms[ires][0], backbone_atoms[ires][1], backbone_atoms[ires][2], backbone_atoms[ires][3])
            gamma=measure.measure_dihedral_angle(backbone_atoms[ires][1], backbone_atoms[ires][2], backbone_atoms[ires][3], backbone_atoms[ires][4])
            delta=measure.measure_dihedral_angle(backbone_atoms[ires][2], backbone_atoms[ires][3], backbone_atoms[ires][4], backbone_atoms[ires][5])
            print("%4d %8.1f %8.1f %8.1f %8.1f    NA       NA   "%(ires, alpha, beta, gamma, delta))
            
    print()
    return(True)
    




    
#mypath='/Users/jli27/Dropbox/jianing/Research/Nano_materials/DNA_aptamer/PDB/Training_Set/'
#mypath='/Users/jli27/Dropbox/jianing/Research/Nano_materials/DNA_aptamer/PDB/B-DNA/'
#mypath='/Users/jli27/Dropbox/jianing/Research/Nano_materials/DNA_aptamer/PDB/A-DNA/'
#mypath='/Users/jli27/Dropbox/jianing/Research/Nano_materials/DNA_aptamer/PDB/Stem_loop/'
#mypath='/Users/jli27/Dropbox/jianing/Research/Nano_materials/DNA_aptamer/PDB/Z-DNA/'
#mypath='/Users/jli27/Dropbox/jianing/Research/Nano_materials/DNA_aptamer/PDB/G-quadruplex/'

mypath='/home/shi631/Public/Angle_dis_ssDNA/'

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

ifile = 0
for fname in onlyfiles:

   # if ifile < 10:
   #     ifile = ifile + 1
   # else:
   #     sys.exit(1)
    if ".ent" not in fname:
        continue
    
    filename=mypath+fname
    print("#Working on "+filename)
    st = structure.StructureReader.read(filename)

    get_backbone_atoms(st)
