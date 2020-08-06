#PBS -S /bin/bash
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=480:00:00
#PBS -l mem=80gb
#PBS -M mfs61958@uga.edu
#PBS -m abe

cd $PBS_O_WORKDIR

echo
echo "Job ID: $PBS_JOBID"
echo "Queue:  $PBS_QUEUE"
echo "Cores:  $PBS_NP"
echo "Nodes:  $(cat $PBS_NODEFILE | sort -u | tr '\n' ' ')"
echo

module load FastQC/0.11.8-Java-1.8.0_144


time fastqc -f fastq *.fastq.gz



# Author: Michael Skaro
# Date: 8/3/2020
# Purpose: Use the fastQC module to complete RNA sequencing data QC
# Done




