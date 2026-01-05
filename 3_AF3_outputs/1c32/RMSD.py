import MDAnalysis as mda
from MDAnalysis.analysis import rms

# Load reference (native) and model structures
u_ref = mda.Universe("native.pdb")
u_mob = mda.Universe("model.pdb")

# Selection:
# - exclude hydrogens (name H*, element H)
# - exclude common ions (extend this list if needed)
ion_resnames = ["NA", "K", "CL", "MG", "CA", "ZN"]
ion_sel = " ".join(ion_resnames)

selection = "nucleic and not name H* and not element H"


ref_atoms = u_ref.select_atoms(selection)
mob_atoms = u_mob.select_atoms(selection)

print(f"Ref atoms: {len(ref_atoms)}, Model atoms: {len(mob_atoms)}")

if len(ref_atoms) != len(mob_atoms):
    raise ValueError("Atom count mismatch after selection – "
                     "check that both PDBs have the same atoms/order "
                     "and add/remove resnames in ion_resnames if needed.")

# RMSD with optimal superposition (Kabsch)
rmsd_value = rms.rmsd(mob_atoms.positions, ref_atoms.positions)

rmsd_no_align = rms.rmsd(mob_atoms.positions, ref_atoms.positions,
                         center=False, superposition=False)

rmsd_align = rms.rmsd(mob_atoms.positions, ref_atoms.positions)
# default: center=True, superposition=True

print("RMSD without alignment:", rmsd_no_align)
print("RMSD with optimal alignment:", rmsd_align)



