# strip_H.tcl
# Usage:
#   vmd -dispdev text -eofexit -e removeH.tcl -args native.pdb native_noH.pdb

# Get arguments
set in_file  [lindex $argv 0]
set out_file [lindex $argv 1]

# Load structure
mol new $in_file waitfor all

# Select all non-hydrogen atoms
set sel [atomselect top "all and noh and not (index 1 2 3)"]

# Write out H-stripped PDB
$sel writepdb $out_file

# Clean up
$sel delete
mol delete top

# Exit VMD
exit


