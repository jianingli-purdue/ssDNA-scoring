#usage:
#vmd -dispdev text -e gro2pdb.tcl


mol new model.gro
set sel [atomselect top "all"]
$sel writepdb model.pdb
exit

