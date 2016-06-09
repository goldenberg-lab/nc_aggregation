rep=0
i=0
for ((complexity=0; complexity<3; complexity=complexity+1))
do
name="myoc"
echo "/hpf/largeprojects/agoldenb/aziz/diseaseMechanism/Rproject/analyse_data.r /hpf/largeprojects/agoldenb/aziz/diseaseMechanism/ /hpf/largeprojects/agoldenb/aziz/myocardial-infarction/ /hpf/largeprojects/agoldenb/aziz/myocardial-infarction/results_$complexity/ 0 0 patients.txt variants1.txt genotype1.txt 1 $complexity 0.9 $i $rep" | qsub -l mem=63G,vmem=63G,nodes=1:ppn=1,walltime=71:50:00 -N $name"_$complexity"
done
