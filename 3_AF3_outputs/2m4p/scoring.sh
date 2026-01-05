for ndx in `seq 0 1 4`
do
echo $ndx
gmx trjconv -f out_$ndx.gro -o out_$ndx.xtc
gmx trjconv -f out_$ndx.xtc -o out_$ndx.trr
gmx trjconv -f out_$ndx.trr -o out_add_$ndx.trr -timestep 1
done
echo -e "1 \n 2 \n 3 \n 4 \n 5" |  gmx trjcat -f out_add_0.trr out_add_1.trr out_add_2.trr out_add_3.trr out_add_4.trr -o out_gro.trr -tu fs -settime

#cp ../4-gdt/references ../4-gdt/dna_0.pdb .
python3 step_43_scoring_trajectory.py
