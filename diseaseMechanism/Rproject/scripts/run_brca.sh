rep=0
name="brca"
for ((a=1; a<5; a=a+1))
do
for ((b=1; b<5; b=b+1))
do
for ((t=1; t<4; t=t+1))
do
for ((i=0; i<1; i=i+1))
do
if [ $a -ne $b ]
then
echo "/hpf/largeprojects/agoldenb/aziz/diseaseMechanism/Rproject/scripts/analyse_brca_new.r /hpf/largeprojects/agoldenb/aziz/diseaseMechanism/ /hpf/largeprojects/agoldenb/aziz/brca/processed/ /hpf/largeprojects/agoldenb/aziz/brca/results/ $a $b $t $i $rep" | qsub -l mem=63G,vmem=63G,nodes=1:ppn=1 -N $name"_$a""_$b""_$t""_$i"
fi
done
done
done
done
