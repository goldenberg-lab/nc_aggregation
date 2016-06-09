#This code suppose a balanced set of cases/controls with 0/1 assignments
#file="/hpf/largeprojects/agoldenb/aziz/bipolar/4-Bipolar_exomes_training.txt"
#output="/hpf/largeprojects/agoldenb/aziz/bipolar/processed/"
#file="/hpf/largeprojects/agoldenb/aziz/epi4k/patients.txt"
#output="/hpf/largeprojects/agoldenb/aziz/epi4k/"
#file="/hpf/largeprojects/agoldenb/aziz/FHS/processed/LDL_High/patients1.txt"
#output="/hpf/largeprojects/agoldenb/aziz/FHS/processed/LDL_High/"
#file="/hpf/largeprojects/agoldenb/aziz/diseaseMechanism/simul/simul_10_0_200/rep1/patients.txt"
#output="/hpf/largeprojects/agoldenb/aziz/diseaseMechanism/simul/simul_10_0_200/rep1/"
file="/hpf/largeprojects/agoldenb/aziz/brca/processed/LumB-Basal/patients_LumB-Basal.txt"
output="/hpf/largeprojects/agoldenb/aziz/brca/processed/LumB-Basal/"

nfold=5
/hpf/largeprojects/agoldenb/aziz/diseaseMechanism/Rproject/scripts/cross_validation_gen.r "$file" "$output" "$nfold" 0 0  #header was 1 for bipolar, paired is 1 for epi4k
name="brca_LumB-Basal"
rep=0
i=0
for ((fold=1; fold<=$nfold; fold=fold+1))
do
for ((complexity=1; complexity<4; complexity=complexity+1))
do
for ((ratioSignal=10; ratioSignal<95; ratioSignal=ratioSignal+10))
do
echo "/hpf/largeprojects/agoldenb/aziz/diseaseMechanism/Rproject/analyse_data.r /hpf/largeprojects/agoldenb/aziz/diseaseMechanism/ $output $output/results_$fold-$complexity-$ratioSignal/ 0 0 patients_training$fold.txt variants.txt genotype.txt 0 $complexity $ratioSignal $i $rep patients_validation$fold.txt" | qsub -l mem=40G,vmem=40G,nodes=1:ppn=1,walltime=239:50:00 -N $name"_$fold""_$complexity""_$ratioSignal"
done
done
done

