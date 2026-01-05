

set sel [atomselect top "all and not (index 4 5 6 7)"]
$sel writegro $argv

exit
