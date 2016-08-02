#!/bin/bash
for j in {1..100}
do
 sbatch me.sh $j
done
