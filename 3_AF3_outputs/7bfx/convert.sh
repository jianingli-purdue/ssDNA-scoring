
ndx=0
for alpha3 in `ls fold_*.cif`
do
    mol2=out_$ndx.gro
    sed -e '/ATOM 1   O OP3/d' -e '/ATOM 2   P P/d' -e '/ATOM 3   O OP1/d' -e '/ATOM 4   O OP2/d' $alpha3 > alpha3_$ndx.cif
    vmd -dispdev text -m alpha3_$ndx.cif -args $mol2 < convert.tcl
    echo $ndx $mol $alpha3 $mol2
    ndx=$((ndx + 1))
    
done
