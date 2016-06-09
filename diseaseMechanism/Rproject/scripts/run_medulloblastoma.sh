rep=10
name="medullo"

for ((a=1; a<5; a=a+1))
do
for ((b=0; b<3; b=b+1))
do
for ((i=1; i<11; i=i+1))
do
echo "/hpf/largeprojects/agoldenb/aziz/diseaseMechanism/Rproject/medulloblastoma.r /hpf/largeprojects/agoldenb/aziz/diseaseMechanism/ $a $b $i $rep" | qsub -l vmem=20G,nodes=1:ppn=1 -N $name"_$a""_$b""_$i"
done
done
done
