for ((i=1; i<4; i=i+1)) do
for ((j=0; j<4; j=j+1)) do
cat "medulloblastoma/results/result"$i"_"$j"_"*"/lik.txt" > "medulloblastoma/results/final/lik"$i"_"$j".txt"
cat "medulloblastoma/results/result"$i"_"$j"_"*"/pred.txt" > "medulloblastoma/results/final/pred"$i"_"$j".txt"
done
done

