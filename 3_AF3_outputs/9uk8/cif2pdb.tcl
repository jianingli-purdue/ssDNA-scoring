#usage:
#vmd -dispdev text -e cif2pdb.tcl


mol new model.cif
set sel [atomselect top "all and not name K"]
$sel writepdb model.pdb
exit

