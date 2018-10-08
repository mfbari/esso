#!/bin/bash

for i in {1..10}; do
  echo $i
  name=$( echo $i*0.003 | bc);
  eval 'python run_simulation.py -tri13129_0.5_0'$name' ../data/7170_0.5_0'$name'/ > heuristic_0'$name'.dat'
done
