

set sel [atomselect top "all not (serial 4 5 6 7)"]
$sel writegro $argv

exit
