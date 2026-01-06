

set sel [atomselect top "all not (serial 10 11 12 13)"]
$sel writegro $argv

exit
