rep=0
i=0
for ((complexity=1; complexity<5; complexity=complexity+1))
do
name="brca_LumB-Basal"
echo "/hpf/largeprojects/agoldenb/aziz/diseaseMechanism/Rproject/analyse_data.r /hpf/largeprojects/agoldenb/aziz/diseaseMechanism/ /hpf/largeprojects/agoldenb/aziz/brca/processed/LumB-Basal/ /hpf/largeprojects/agoldenb/aziz/brca/processed/LumB-Basal/results_$complexity/ gene_expression_hareem.txt 0 patients.txt variants.txt genotype.txt 1 $complexity 0.9 $i $rep" | qsub -l mem=63G,vmem=63G,nodes=1:ppn=1,walltime=239:50:00 -N $name"_$complexity"
done
