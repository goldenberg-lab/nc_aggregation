rep=0
i=0
for ((complexity=0; complexity<4; complexity=complexity+1))
do
name="eomi"
echo "/hpf/largeprojects/agoldenb/aziz/diseaseMechanism/Rproject/analyse_data.r /hpf/largeprojects/agoldenb/aziz/diseaseMechanism/ /hpf/largeprojects/agoldenb/aziz/FHS/processed/EOMI/ /hpf/largeprojects/agoldenb/aziz/FHS/processed/EOMI/results_$complexity/ 0 0 patients.txt variants1.txt genotype.txt 0 $complexity 0.9 $i $rep" | qsub -l mem=63G,vmem=63G,nodes=1:ppn=1,walltime=239:50:00 -N $name"_$complexity"
done
