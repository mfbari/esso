#!/bin/bash

for s in {1..10}; do
  arr_rate=$( echo $s*0.01 | bc);
  algo='f'
  dataset=$(printf '13129_%4.2f_3600' "$arr_rate")
  eval 'python run_simulation.py -'$algo'ri'$dataset' ../data/'$dataset'/ &> '$algo'_ocf_output_'$dataset'.dat &'
  echo $!
done
