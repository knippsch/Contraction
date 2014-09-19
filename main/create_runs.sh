#!/bin/bash

outpath="/hiskp2/werner/LapH/correlators"
lattice="A100_test"

start_config=1685
end_config=1837
delta_config=8
stepping=1

number_of_max_mom=0
max_mom_in_one_dir=0
dirac_min=00
dirac_max=05
displ_min=0
displ_max=0


filenumber=1

mkdir -p $lattice
cd $lattice
mkdir -p infiles_$lattice
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
mkdir -p $outpath
mkdir -p $outpath"/"$lattice

for dirac1 in $(seq $dirac_min 01 $dirac_max); do
  while [ ${#dirac1} -ne 2 ]; do
    dirac1="0"$dirac1;
  done
  for dirac2 in $(seq $dirac_min 01 $dirac_max); do
    while [ ${#dirac2} -ne 2 ]; do
      dirac2="0"$dirac2;
    done
    for p1 in $(seq 0 1 $(($number_of_max_mom* $number_of_max_mom))); do
      for p2 in $(seq 0 1 $(($number_of_max_mom* $number_of_max_mom))); do
        for displ1 in $(seq $displ_min 1 $displ_max); do
          for displ2 in $(seq $displ_min 1 $displ_max); do
            if [ "$p1" -eq "$p2" ]; then
              mkdir -p $outpath"/"$lattice"/dirac_"$dirac1"_"$dirac2"_p_"$p1"_"$p2"_displ_"$displ1"_"$displ2"_unsuppressed";
            fi
            mkdir -p $outpath"/"$lattice"/dirac_"$dirac1"_"$dirac2"_p_"$p1"_"$p2"_displ_"$displ1"_"$displ2;
          done;
        done;
      done;
    done;
  done;
done

