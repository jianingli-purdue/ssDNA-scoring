

set sel [atomselect top "all and not (serial 3 4 5 6)"]
$sel writegro $argv

exit
