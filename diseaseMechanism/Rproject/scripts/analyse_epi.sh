rep=0
i=0
for ((complexity=0; complexity<4; complexity=complexity+1))
do
name="epi"
echo "/hpf/largeprojects/agoldenb/aziz/diseaseMechanism/Rproject/analyse_data.r /hpf/largeprojects/agoldenb/aziz/diseaseMechanism/ /hpf/largeprojects/agoldenb/aziz/epi4k/ /hpf/largeprojects/agoldenb/aziz/epi4k/results_$complexity/ 0 0 patients.txt variants.txt genotype.txt 1 $complexity $i $rep" | qsub -l mem=63G,vmem=63G,nodes=1:ppn=1,walltime=71:50:00 -N $name"_$complexity"
done
