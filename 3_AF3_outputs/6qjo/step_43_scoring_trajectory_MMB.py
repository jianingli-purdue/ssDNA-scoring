import os
import glob
import numpy as np
import warnings
warnings.filterwarnings('ignore')

import MDAnalysis as mda
from MDAnalysis.analysis import rms, align

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------
ref_pdb = "last.1.pdb"      # reference structure in working directory
ref_dir = "./reference"     # directory containing .gro files
output_log = "MMB_GDT.txt"  # output file for GDT scores
cutoffs = [2.0, 3.0, 4.0, 5.0]  # Å, same as original script

# ------------------------------------------------------------
# Load reference structure
# ------------------------------------------------------------
print(f"[INFO] Loading reference PDB: {ref_pdb}")
ref = mda.Universe(ref_pdb)

# Use phosphate atoms in DNA residues (same as original)
ref_P = ref.select_atoms("name P and resname D*")
if len(ref_P) == 0:
    raise ValueError("No atoms selected with 'name P and resname D*' in last.1.pdb")

# Centered coordinates for alignment
ref_P0 = ref_P.positions - ref_P.center_of_mass()

# ------------------------------------------------------------
# Collect all .gro files in reference directory
# ------------------------------------------------------------
gro_files = sorted(glob.glob(os.path.join(ref_dir, "*.gro")))
if not gro_files:
    raise FileNotFoundError(f"No .gro files found in {ref_dir}")

print("[INFO] Found the following .gro files:")
for g in gro_files:
    print("  ", g)

# ------------------------------------------------------------
# Open output log
# ------------------------------------------------------------
try:
    log_fh = open(output_log, "w")
    print(f"[INFO] Writing GDT scores to: {output_log}")
except Exception as e:
    print(f"[WARN] Could not open {output_log} for writing: {e}")
    log_fh = None

# Header line (optional)
header = (
    "# index  file  GDT_2A  GDT_3A  GDT_4A  GDT_5A  GDT_main(3A)\n"
)
if log_fh is not None:
    log_fh.write(header)

# ------------------------------------------------------------
# Main loop: compute GDT between last.1.pdb and each .gro
# ------------------------------------------------------------
for idx, gro_file in enumerate(gro_files):
    print(f"\n[INFO] Processing {idx}: {gro_file}")
    u = mda.Universe(gro_file)

    tag = u.atoms
    # In original script: tag_P = select_atoms("name P")
    tag_P = u.select_atoms("name P")

    if len(tag_P) != len(ref_P):
        print(
            f"[WARN] Number of P atoms in {gro_file} ({len(tag_P)}) "
            f"does not match reference ({len(ref_P)}); skipping."
        )
        continue

    # Superposition & alignment (same pattern as original)
    _P_rmsd = rms.rmsd(ref_P.positions, tag_P.positions, superposition=True)

    tag_P0 = tag_P.positions - tag_P.center_of_mass()
    R, _rmsd = align.rotation_matrix(tag_P0, ref_P0)

    tag.translate(-tag_P.center_of_mass())
    tag.rotate(R)
    tag.translate(ref_P.center_of_mass())

    # Re-select after transformation
    tag_P1 = tag.select_atoms("name P")

    # Distances between corresponding P atoms
    dist_list = []
    for ref_atom, tag_atom in zip(ref_P, tag_P1):
        p1 = ref_atom.position
        p2 = tag_atom.position
        dist = np.linalg.norm(p1 - p2)
        dist_list.append(dist)

    dist_array = np.asarray(dist_list)

    # GDT fractions at given cutoffs (same logic as original)
    gdts = []
    for cutoff in cutoffs:
        gdts.append((dist_array <= cutoff).sum() / float(len(dist_array)))

    # Use 3 Å cutoff (index 1) as main GDT score, same as original script
    GDT_main = gdts[1]

    print(
        f"  {os.path.basename(gro_file)}  "
        + " ".join(f"{g:0.4f}" for g in gdts)
        + f"  GDT(3Å) = {GDT_main:0.4f}"
    )

    # Write to log
    if log_fh is not None:
        line = (
            f"{idx:4d}  {os.path.basename(gro_file)}  "
            + "  ".join(f"{g:0.4f}" for g in gdts)
            + f"  {GDT_main:0.4f}\n"
        )
        log_fh.write(line)
        log_fh.flush()

# ------------------------------------------------------------
# Close log file
# ------------------------------------------------------------
if log_fh is not None:
    log_fh.close()
    print(f"\n[INFO] Finished. Results saved in {output_log}")















