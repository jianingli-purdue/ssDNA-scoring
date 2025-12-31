from schrodinger import structure
from schrodinger.structutils import rmsd, measure

from os import listdir
from os.path import isfile, join

import sys
import math


def get_backbone_atoms(st, chi_list=None):
    """
    Collect backbone atoms and compute backbone dihedrals (alpha–zeta)
    plus glycosidic chi angle for each residue in chain A.

    chi_list (optional): list to which all chi values (in degrees) will be appended.
    """
    backbone_atoms = []
    resnum = 0

    for res in st.residue:
        # res_atoms:
        # 0: P, 1: O5', 2: C5', 3: C4', 4: C3', 5: O3'
        # 6: O4', 7: C1', 8: base N (N9 or N1), 9: base C (C4 or C2)
        res_atoms = [None] * 10
        if res.chain == "A" and res.pdbres.strip() in ['DA', 'DC', 'DG', 'DT']:
            resname = res.pdbres.strip()

            # Decide which base atoms we need for chi depending on residue type
            if resname in ['DA', 'DG']:
                chi_N_name = "N9"
                chi_C_name = "C4"
            else:  # DC or DT
                chi_N_name = "N1"
                chi_C_name = "C2"

            for at in res.atom:
                aname = at.pdbname.strip()
                if aname == "P":
                    res_atoms[0] = at
                elif aname == "O5'":
                    res_atoms[1] = at
                elif aname == "C5'":
                    res_atoms[2] = at
                elif aname == "C4'":
                    res_atoms[3] = at
                elif aname == "C3'":
                    res_atoms[4] = at
                elif aname == "O3'":
                    res_atoms[5] = at
                elif aname == "O4'":
                    res_atoms[6] = at
                elif aname == "C1'":
                    res_atoms[7] = at
                elif aname == chi_N_name:
                    res_atoms[8] = at
                elif aname == chi_C_name:
                    res_atoms[9] = at

            backbone_atoms.append(res_atoms)
            resnum += 1

    print("#ires   alpha    beta     gamma     delta    epsilon   zeta      chi")
    for ires in range(resnum):
        # Compute chi for this residue if all four atoms are present
        chi = None
        o4p = backbone_atoms[ires][6]
        c1p = backbone_atoms[ires][7]
        nbase = backbone_atoms[ires][8]
        cbase = backbone_atoms[ires][9]
        if o4p is not None and c1p is not None and nbase is not None and cbase is not None:
            chi = measure.measure_dihedral_angle(o4p, c1p, nbase, cbase)
            if chi_list is not None:
                chi_list.append(chi)

        # Helper to format float or NA
        def fmt(angle):
            if angle is None:
                return "      NA"
            else:
                return f"{angle:8.1f}"

        if ires > 0 and ires < resnum - 1:
            alpha   = measure.measure_dihedral_angle(backbone_atoms[ires-1][5],
                                                     backbone_atoms[ires][0],
                                                     backbone_atoms[ires][1],
                                                     backbone_atoms[ires][2])
            beta    = measure.measure_dihedral_angle(backbone_atoms[ires][0],
                                                     backbone_atoms[ires][1],
                                                     backbone_atoms[ires][2],
                                                     backbone_atoms[ires][3])
            gamma   = measure.measure_dihedral_angle(backbone_atoms[ires][1],
                                                     backbone_atoms[ires][2],
                                                     backbone_atoms[ires][3],
                                                     backbone_atoms[ires][4])
            delta   = measure.measure_dihedral_angle(backbone_atoms[ires][2],
                                                     backbone_atoms[ires][3],
                                                     backbone_atoms[ires][4],
                                                     backbone_atoms[ires][5])
            epsilon = measure.measure_dihedral_angle(backbone_atoms[ires][3],
                                                     backbone_atoms[ires][4],
                                                     backbone_atoms[ires][5],
                                                     backbone_atoms[ires+1][0])
            zeta    = measure.measure_dihedral_angle(backbone_atoms[ires][4],
                                                     backbone_atoms[ires][5],
                                                     backbone_atoms[ires+1][0],
                                                     backbone_atoms[ires+1][1])
            print(f"{ires:4d} {alpha:8.1f} {beta:8.1f} {gamma:8.1f} {delta:8.1f} "
                  f"{epsilon:8.1f} {zeta:8.1f} {fmt(chi)}")
        elif ires == 0:
            gamma   = measure.measure_dihedral_angle(backbone_atoms[ires][1],
                                                     backbone_atoms[ires][2],
                                                     backbone_atoms[ires][3],
                                                     backbone_atoms[ires][4])
            delta   = measure.measure_dihedral_angle(backbone_atoms[ires][2],
                                                     backbone_atoms[ires][3],
                                                     backbone_atoms[ires][4],
                                                     backbone_atoms[ires][5])
            epsilon = measure.measure_dihedral_angle(backbone_atoms[ires][3],
                                                     backbone_atoms[ires][4],
                                                     backbone_atoms[ires][5],
                                                     backbone_atoms[ires+1][0])
            zeta    = measure.measure_dihedral_angle(backbone_atoms[ires][4],
                                                     backbone_atoms[ires][5],
                                                     backbone_atoms[ires+1][0],
                                                     backbone_atoms[ires+1][1])
            # alpha and beta not defined for first residue
            print(f"{ires:4d}    NA       NA {gamma:8.1f} {delta:8.1f} "
                  f"{epsilon:8.1f} {zeta:8.1f} {fmt(chi)}")
        elif ires == resnum - 1:
            alpha = measure.measure_dihedral_angle(backbone_atoms[ires-1][5],
                                                   backbone_atoms[ires][0],
                                                   backbone_atoms[ires][1],
                                                   backbone_atoms[ires][2])
            beta  = measure.measure_dihedral_angle(backbone_atoms[ires][0],
                                                   backbone_atoms[ires][1],
                                                   backbone_atoms[ires][2],
                                                   backbone_atoms[ires][3])
            gamma = measure.measure_dihedral_angle(backbone_atoms[ires][1],
                                                   backbone_atoms[ires][2],
                                                   backbone_atoms[ires][3],
                                                   backbone_atoms[ires][4])
            delta = measure.measure_dihedral_angle(backbone_atoms[ires][2],
                                                   backbone_atoms[ires][3],
                                                   backbone_atoms[ires][4],
                                                   backbone_atoms[ires][5])
            # epsilon and zeta not defined for last residue
            print(f"{ires:4d} {alpha:8.1f} {beta:8.1f} {gamma:8.1f} {delta:8.1f} "
                  f"   NA       NA {fmt(chi)}")

    print()
    return True


#mypath='/Users/jli27/Dropbox/jianing/Research/Nano_materials/DNA_aptamer/PDB/Training_Set/'
#mypath='/Users/jli27/Dropbox/jianing/Research/Nano_materials/DNA_aptamer/PDB/B-DNA/'
#mypath='/Users/jli27/Dropbox/jianing/Research/Nano_materials/DNA_aptamer/PDB/A-DNA/'
#mypath='/Users/jli27/Dropbox/jianing/Research/Nano_materials/DNA_aptamer/PDB/Stem_loop/'
#mypath='/Users/jli27/Dropbox/jianing/Research/Nano_materials/DNA_aptamer/PDB/Z-DNA/'
#mypath='/Users/jli27/Dropbox/jianing/Research/Nano_materials/DNA_aptamer/PDB/G-quadruplex/'

mypath = '/home/shi631/Public/Angle_dis_ssDNA/'

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

chi_all = []  # collect chi angles over all structures

ifile = 0
for fname in onlyfiles:
    # if ifile < 10:
    #     ifile = ifile + 1
    # else:
    #     sys.exit(1)
    if ".ent" not in fname:
        continue

    filename = mypath + fname
    print("#Working on " + filename)
    st = structure.StructureReader.read(filename)

    get_backbone_atoms(st, chi_all)

# Simple chi distribution summary
if chi_all:
    n = len(chi_all)
    mean_chi = sum(chi_all) / n
    var_chi = sum((x - mean_chi) ** 2 for x in chi_all) / n
    std_chi = math.sqrt(var_chi)

    print("#Chi angle distribution over all residues")
    print(f"# N = {n}")
    print(f"# mean = {mean_chi:8.2f} deg")
    print(f"# std  = {std_chi:8.2f} deg")

