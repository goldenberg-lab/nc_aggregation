#!/bin/bash

# $1 is complexity
# $2 is ratioSignal
# $3 is i
# $4 is rep

/hpf/largeprojects/agoldenb/aziz/diseaseMechanism/Rproject/analyse_data.r /hpf/largeprojects/agoldenb/aziz/diseaseMechanism/ /hpf/largeprojects/agoldenb/aziz/FHS/processed/LDL/ /hpf/largeprojects/agoldenb/aziz/FHS/processed/LDL/results_$1/ 0 0 patients.txt variants1.txt genotype.txt 0 $1 $2 $3 $4

