dirname="/hpf/largeprojects/agoldenb/aziz/diseaseMechanism/simul/simul_10_0_100/rep"
for ((i=1; i<=20; i=i+1))
do
for modnet in -0.5 -0.2 0.2 0.5
do
echo "/hpf/largeprojects/agoldenb/aziz/diseaseMechanism/Rproject/scripts/test_network_robustness.r $dirname""$i"/" $modnet" | qsub -l mem=30G,vmem=30G,nodes=1:ppn=1,walltime=239:50:00 -N "net_""$i""_""$modnet"
done
done


