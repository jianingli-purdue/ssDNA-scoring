import MDAnalysis as mda
from MDAnalysis.analysis import align, rms

# --- Load universes ---
u_ref = mda.Universe("8y0f_9.pdb")   # reference
u_mob = mda.Universe("aligned.pdb")    # model to align

# --- Selection: heavy atoms, no ions (no element keyword!) ---
ion_resnames = ["NA", "K", "CL", "MG", "CA", "ZN"]
ion_sel = " ".join(ion_resnames)

# no hydrogens, no ions
selection = f"not name H* and not resname {ion_sel}"

# If you also want to exclude water, uncomment this instead:
# selection = (
#     f"not name H* and not resname {ion_sel} "
#     "and not resname HOH TIP3 TIP3P SOL WAT"
# )

ref_atoms = u_ref.select_atoms(selection)
mob_atoms = u_mob.select_atoms(selection)

print(f"Ref atoms: {len(ref_atoms)}, Model atoms: {len(mob_atoms)}")

if len(ref_atoms) != len(mob_atoms):
    raise ValueError(
        "Atom count mismatch after selection.\n"
        "Check that native.pdb and model.pdb have the same atoms/order,\n"
        "and adjust ion_resnames if they use different names."
    )

# --- RMSD BEFORE alignment (raw coordinates) ---
rmsd_before = rms.rmsd(
    mob_atoms.positions,
    ref_atoms.positions,
    center=False,
    superposition=False,
)
print(f"RMSD BEFORE alignment (heavy, no ions): {rmsd_before:.3f} Å")

# --- Align model -> native using the same selection ---
# This modifies u_mob coordinates in place; returns (rmsd, rotation_matrix)
align.alignto(u_mob, u_ref, select=selection)

# Rebuild selections (optional but safe)
ref_atoms = u_ref.select_atoms(selection)
mob_atoms = u_mob.select_atoms(selection)

# --- RMSD AFTER alignment (no extra superposition) ---
rmsd_after = rms.rmsd(
    mob_atoms.positions,
    ref_atoms.positions,
    center=False,
    superposition=False,
)
print(f"RMSD AFTER alignment (heavy, no ions):  {rmsd_after:.3f} Å")




