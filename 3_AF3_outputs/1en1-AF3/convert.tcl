

set sel [atomselect top "all not (serial 1 2 3 4)"]
$sel writegro $argv

exit
