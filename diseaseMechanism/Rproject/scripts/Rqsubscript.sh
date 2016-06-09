#!/bin/bash

rep=0
i=0
name="ldl"
for ((complexity=1; complexity<2; complexity=complexity+1))
do

qsub -l mem=63G,vmem=63G,nodes=1:ppn=1,walltime=239:50:00 -N $name"_$complexity" /hpf/largeprojects/agoldenb/aziz/diseaseMechanism/Rproject/scripts/analyse_ldl.sh -F "$complexity $i $rep"

done

