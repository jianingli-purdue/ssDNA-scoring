#!/usr/bin/env bash
# Split a (multi-model) PDB into per-frame .gro files via VMD
# Downloads PDB (uppercase), saves as ./dna_0.pdb and ./reference/{PDBID}.pdb,
# then uses VMD to write frames to ./reference/dna_0_<n>.gro
# Usage: ./split_pdb_to_gro.sh <pdb_id>

# usage:

#chmod +x split_pdb_to_gro.sh
#./split_pdb_to_gro.sh 8or8
#



set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <pdb_id>"
  exit 1
fi

# Dependencies check
command -v wget >/dev/null 2>&1 || { echo "[ERROR] wget not found"; exit 1; }
command -v vmd  >/dev/null 2>&1 || { echo "[ERROR] VMD not found";  exit 1; }

# Normalize ID
pdb_lower="$(echo "$1" | tr '[:upper:]' '[:lower:]')"
pdb_upper="$(echo "$1" | tr '[:lower:]' '[:upper:]')"

# 1) Make reference dir in current directory
mkdir -p reference

# 2) Download PDB once to a temp file (uppercase in URL)
tmp_pdb="$(mktemp -t "${pdb_upper}_XXXXXX.pdb")"
echo "[INFO] Downloading ${pdb_upper}.pdb ..."
wget -q "https://files.rcsb.org/download/${pdb_upper}.pdb" -O "${tmp_pdb}"

# 3) Save to current directory as dna_0.pdb and to reference/{PDBID}.pdb
cp -f "${tmp_pdb}" "./dna_0.pdb"
cp -f "${tmp_pdb}" "reference/${pdb_upper}.pdb"
rm -f "${tmp_pdb}"
echo "[INFO] Saved ./dna_0.pdb and ./reference/${pdb_upper}.pdb"

# 4) Use VMD to split frames -> reference/dna_0_<frame>.gro
tcl_script="$(mktemp -t vmd_split_XXXXXX.tcl)"
trap 'rm -f "$tcl_script"' EXIT

cat > "$tcl_script" <<'TCL'
# VMD TCL script
# Args: <pdb_file>
if {$argc < 1} {
    puts "Usage: vmd -dispdev text -e script.tcl -args <pdb_file>"
    quit
}
set pdbfile [lindex $argv 0]
mol new "$pdbfile" type pdb waitfor all
set molid [molinfo top]
set nf    [molinfo $molid get numframes]
set sel   [atomselect $molid "all"]

file mkdir "reference"

for {set i 0} {$i < $nf} {incr i} {
    animate goto $i
    $sel frame $i
    set out [format "reference/dna_0_%d.gro" $i]
    animate write gro "$out" sel $sel
    puts "Wrote $out"
}
quit
TCL

echo "[INFO] Splitting frames from ./dna_0.pdb into GRO files under ./reference ..."
vmd -dispdev text -e "$tcl_script" -args "./dna_0.pdb" -eofexit >/dev/null

echo "[DONE] GRO files saved in ./reference (dna_0_0.gro, dna_0_1.gro, ...)"

