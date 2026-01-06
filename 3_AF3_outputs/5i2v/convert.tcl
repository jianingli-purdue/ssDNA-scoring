

set sel [atomselect top "all not (serial 3 4 5 6)"]
$sel writegro $argv

exit
