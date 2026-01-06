

set sel [atomselect top "all not (serial 5 6 7 8)"]
$sel writegro $argv

exit
