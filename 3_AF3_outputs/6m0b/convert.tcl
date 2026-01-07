

set sel [atomselect top "all not (serial 8 9 10 11)"]
$sel writegro $argv

exit
