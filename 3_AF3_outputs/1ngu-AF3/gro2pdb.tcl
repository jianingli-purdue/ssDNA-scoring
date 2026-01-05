#usage:
#vmd -dispdev text -e gro2pdb.tcl


mol new native.gro
set sel [atomselect top "all"]
$sel writepdb native.pdb
exit

