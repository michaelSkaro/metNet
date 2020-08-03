#PBS -S /bin/bash
#PBS -q batch
#PBS -N mfs_multiqc
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00
#PBS -l mem=20g
#PBS -M mfs61958@uga.edu
#PBS -m abe

cd $PBS_O_WORKDIR

echo
echo "Job ID: $PBS_JOBID"
echo "Queue:  $PBS_QUEUE"
echo "Cores:  $PBS_NP"
echo "Nodes:  $(cat $PBS_NODEFILE | sort -u | tr '\n' ' ')"
echo

module load MultiQC/1.5-foss-2016b-Python-2.7.14


multiqc .



# Author: Michael Skaro
# Date: 8/3/2020
# Purpose: Use the multiQC module to complete RNA sequencing data multiQC report
# Done




