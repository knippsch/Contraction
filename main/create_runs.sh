#!/bin/bash

outpath="/hiskp2/werner/LapH/correlators_test"
lattice="A40.20"

start_config=1300
end_config=1300
delta_config=6
stepping=1

number_of_max_mom=2
max_mom_in_one_dir=2
dirac_min=05
dirac_max=05
displ_min=0
displ_max=0


filenumber=1

mkdir $lattice
cd $lattice
mkdir infiles_$lattice
cd infiles_$lattice

cp ../../infiles/LapHs$lattice.in .

for i in $(seq $start_config $(($stepping * delta_config)) $end_config); do
  echo "creating infile $filenumber"

  perl -pe "s/.*/start_config = "$i" / if $. == 10" < LapHs$lattice.in > LapHs_temp.in
  perl -pe "s/.*/end_config = "$(($i + delta_config * (stepping - 1)))" / if $. == 11" < LapHs_temp.in > LapHs_$filenumber.in
  perl -pe "s/.*/number_of_max_mom = "$number_of_max_mom" / if $. == 26" < LapHs_$filenumber.in > LapHs_temp.in
  perl -pe "s/.*/max_mom_in_one_dir = "$max_mom_in_one_dir" / if $. == 27" < LapHs_temp.in > LapHs_$filenumber.in
  perl -pe "s/.*/dirac_min = "$dirac_min" / if $. == 29" < LapHs_$filenumber.in > LapHs_temp.in
  perl -pe "s/.*/dirac_max = "$dirac_max" / if $. == 30" < LapHs_temp.in > LapHs_$filenumber.in
  perl -pe "s/.*/displ_min = "$displ_min" / if $. == 32" < LapHs_$filenumber.in > LapHs_temp.in
  perl -pe "s/.*/displ_max = "$displ_max" / if $. == 33" < LapHs_temp.in > LapHs_$filenumber.in

  filenumber=$(($filenumber + 1))
done

cd ../../
mkdir $outpath
mkdir $outpath"/"$lattice
mkdir $outpath"/"$lattice"/dirac_"$dirac_min"_"$dirac_max"_p_0_"$number_of_max_mom"_displ_"$displ_min"_"$displ_max
mkdir $outpath"/"$lattice"/dirac_"$dirac_min"_"$dirac_max"_p_0_0_displ_"$displ_min"_"$displ_max
