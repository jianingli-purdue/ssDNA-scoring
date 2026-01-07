

set sel [atomselect top "all not (serial 2 3 4 5)"]
$sel writegro $argv

exit
