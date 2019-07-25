#!/bin/bash

for i in {1..10}; do
  for j in {3600..3600..1800}; do
    arrival_rate=$( echo $i*0.01 | bc);
    algo='c'
    dataset=$(printf '13129_%4.2f_%d' "$arrival_rate" "$j")
    eval 'python run_simulation.py -'$algo'ri'$dataset'_3 ../data/'$dataset'/ &> '$algo'_running_time_'$dataset'.dat &'
    echo $!
  done
done
