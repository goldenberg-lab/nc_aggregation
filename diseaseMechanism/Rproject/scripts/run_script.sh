rep=1
name="disMec"
n=50

for m in 10 20 50
do
for ((i=1; i<21; i=i+$rep ))
do
echo "/hpf/largeprojects/agoldenb/aziz/diseaseMechanism/Rproject/main_Simul_Pop\
.r /hpf/largeprojects/agoldenb/aziz/diseaseMechanism/ $m 0 $n 1 $i $rep" | q\
sub -l vmem=16G,nodes=1:ppn=1,walltime=23:50:00 -N $name"_$n""_$m""_$i"
done
done

