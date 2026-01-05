

set sel [atomselect top "all and not (index 1 2 3 4)"]
$sel writegro $argv

exit
