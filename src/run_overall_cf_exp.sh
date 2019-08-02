#!/bin/bash

for s in {1..10}; do
  std_dev=$( echo $s*0.01 | bc);
  algo='f'
  dataset=$(printf '13129_0.0001_3600_%4.2f' "$std_dev")
  eval 'python run_simulation.py -'$algo'ri'$dataset' ../data/'$dataset'/ &> '$algo'_ocf_run_time'$dataset'.dat &'
  echo $!
done
