#!/bin/bash
#SBATCH --job-name=BA-Graph-State
#SBATCH -p normal
#SBATCH -t 13-23:59:59        <---- time in hours:minutes:seconds
#SBATCH --mem 1gb	      <---- memory
hostname

nodes=20
parameter=1
squeezing=3
repetitions=1
random=0
id=1234
pathway=/Users/mattiauser/ownCloud/ComplexNetworks/BA_n${nodes}_m${parameter}_sq${squeezing}_randomSub${random}_ID${id}
#pathway=/users/jussieu/walschaers/Graph-States/BA_n${nodes}_m${parameter}_sq${squeezing}_randomSub${random}_ID$SLURM_JOBID

#source activate mattia_graph_env

mkdir -p $pathway
cd $pathway

python3 /Users/mattiauser/ownCloud/ComplexNetworks/IntensityCorrBA_Cluster_Stat.py $id $nodes $parameter $squeezing $repetitions $random
