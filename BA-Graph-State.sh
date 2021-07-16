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
pathway=[INSERT LOCAL DIRECTORY PATH]BA_n${nodes}_m${parameter}_sq${squeezing}_randomSub${random}_ID${id}

mkdir -p $pathway
cd $pathway

python3 [INSERT LOCAL DIRECTORY PATH]IntensityCorrBA_Cluster_Stat.py $id $nodes $parameter $squeezing $repetitions $random
