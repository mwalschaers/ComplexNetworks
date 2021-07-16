#!/bin/bash
#SBATCH --job-name=WS-Graph-State
#SBATCH -p normal
#SBATCH -t 5-23:59:59           <---- time in hours:minutes:seconds
#SBATCH --mem 12gb	      <---- memory
hostname

nodes=1000
distance=5
rewireProb=0.2
dim=1
squeezing=15
repetitions=10
random=0
pathway=/users/jussieu/walschaers/Graph-States/WS_n${nodes}_nei${distance}_dim${dim}_p${rewireProb}_sq${squeezing}_randomSub${random}_ID$SLURM_JOBID

source activate mattia_graph_env

mkdir -p $pathway
cd $pathway

python3.7 /users/jussieu/walschaers/Graph-States/IntensityCorrWS_Cluster_Stat.py $SLURM_JOBID $nodes $distance $rewireProb $dim $squeezing $repetitions $random
