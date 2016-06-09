rep=0
i=0
name="bipolar"
for ((complexity=1; complexity<4; complexity=complexity+1))
do
echo "/hpf/largeprojects/agoldenb/aziz/diseaseMechanism/Rproject/analyse_data.r /hpf/largeprojects/agoldenb/aziz/diseaseMechanism/ /hpf/largeprojects/agoldenb/aziz/bipolar/processed/ /hpf/largeprojects/agoldenb/aziz/bipolar/processed/results_old-$complexity/ 0 0 patients_training.txt variants.txt genotype.txt 0 $complexity 0.9 $i $rep patients_validation.txt" | qsub -l mem=63G,vmem=63G,nodes=1:ppn=1,walltime=239:50:00 -N $name"_old_$complexity"
done

