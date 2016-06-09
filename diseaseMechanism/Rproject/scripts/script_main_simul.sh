#PBS -N diseaseMechanism
#PBS -l vmem=12G,nodes=1:ppn=10

/hpf/largeprojects/agoldenb/aziz/diseaseMechanism/Rproject/main_Simul_Pop.r /hpf/largeprojects/agoldenb/aziz/diseaseMechanism/ 10 1 200 1 1 1
