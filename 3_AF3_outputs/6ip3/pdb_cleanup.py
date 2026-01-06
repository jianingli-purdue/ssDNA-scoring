#!/usr/bin/env python3
import sys

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} in.pdb out_clean.pdb")
    sys.exit(1)

in_pdb = sys.argv[1]
out_pdb = sys.argv[2]

# Keep altloc blank or 'A'; adjust if you want only A
KEEP_ALT = (" ", "A")

seen_keys = set()

with open(in_pdb) as fin, open(out_pdb, "w") as fout:
    for line in fin:
        if line.startswith(("ATOM", "HETATM")):
            # PDB format:
            # col 17 (index 16) = altLoc
            altloc = line[16]

            # Build a key that identifies the "logical atom":
            # serial (6–11), atom name (12–16), residue info (17–26)
            key = (line[6:11], line[12:16], line[17:26])

            # Skip unwanted altlocs
            if altloc not in KEEP_ALT:
                continue

            # If we've already written this logical atom, skip duplicates
            if key in seen_keys:
                continue

            seen_keys.add(key)

        fout.write(line)

