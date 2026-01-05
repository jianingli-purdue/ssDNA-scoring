#usage:
#vmd -dispdev text -e cif2pdb.tcl


mol new model.cif
set sel [atomselect top "all"]
$sel writepdb model.pdb
exit

