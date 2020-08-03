#PBS -S /bin/bash
#PBS -N mfs_star
#PBS -q batch
#PBS -l nodes=1:ppn=8:AMD
#PBS -l walltime=480:00:00
#PBS -l mem=80gb

cd $PBS_O_WORKDIR

echo
echo "Job ID: $PBS_JOBID"
echo "Queue:  $PBS_QUEUE"
echo "Cores:  $PBS_NP"
echo "Nodes:  $(cat $PBS_NODEFILE | sort -u | tr '\n' ' ')"
echo


ml STAR/2.6.1c-foss-2016b


### Step 2: Alignment 1st Pass.

STAR
--genomeDir /work/jaalab/mskaro1/rcc
--readFilesIn PM1010_Z3_1_Case_RNASeq_S15_L003_R1.fastq.gz, PM1010_Z3_1_Case_RNASeq_S15_L003_R2.fastq.gz, PM1116_Z3_1_Case_RNASeq_S46_L006_R1_001.fastq.gz, PM1116_Z3_1_Case_RNASeq_S46_L006_R2_001.fastq.gz, PM1220_Z5_1_Case_RNASeq_S17_L003_R1_001.fastq.gz, PM1220_Z5_1_Case_RNASeq_S17_L003_R2_001.fastq.gz, PM1306_Z3_1_Case_RNASeq_S21_L003_R1_001.fastq.gz, PM1306_Z3_1_Case_RNASeq_S21_L003_R2_001.fastq.gz, PM1321_Z1_1_Case_RNASeq_S23_L003_R1_001.fastq.gz, PM1321_Z1_1_Case_RNASeq_S23_L003_R2_001.fastq.gz, PM492_Z5_1_Case_RNASeq_S6_L004_R1_001.fastq.gz, PM492_Z5_1_Case_RNASeq_S6_L004_R2_001.fastq.gz, PM492_Z6_1_Case_RNASeq_S7_L004_R1_001.fastq.gz, PM492_Z6_1_Case_RNASeq_S7_L004_R2_001.fastq.gz, PM492_Z7_1_Case_RNASeq_S8_L004_R1_001.fastq.gz, PM492_Z7_1_Case_RNASeq_S8_L004_R2_001.fastq.gz, Sample_PM1035_X1_1_Case_RNASeq_R1.fastq.gz, Sample_PM1035_X1_1_Case_RNASeq_R2.fastq.gz, Sample_PM1041_Z2_1_Case_RNASeq_R1.fastq.gz, Sample_PM1041_Z2_1_Case_RNASeq_R2.fastq.gz, Sample_PM157_Z7_1_RNA_R1.fastq.gz, Sample_PM157_Z7_1_RNA_R2.fastq.gz, Sample_PM215_Z4_1_Case_RNASeq_R1.fastq.gz, Sample_PM215_Z4_1_Case_RNASeq_R2.fastq.gz, Sample_PM255_Z1_1_Case_RNASeq_R1.fastq.gz, Sample_PM255_Z1_1_Case_RNASeq_R2.fastq.gz, Sample_PM269_Z1_1_Case_RNASeq_R1.fastq.gz, Sample_PM269_Z1_1_Case_RNASeq_R2.fastq.gz, Sample_PM316_Z17_1_Case_RNASeq_R1.fastq.gz, Sample_PM316_Z17_1_Case_RNASeq_R2.fastq.gz, Sample_PM316_Z18_1_Case_RNASeq_R1.fastq.gz, Sample_PM316_Z18_1_Case_RNASeq_R2.fastq.gz, Sample_PM316_Z19_1_Case_RNASeq_R1.fastq.gz, Sample_PM316_Z19_1_Case_RNASeq_R2.fastq.gz, Sample_PM316_Z22_1_Case_RNASeq_R1.fastq.gz, Sample_PM316_Z22_1_Case_RNASeq_R2.fastq.gz, Sample_PM316_Z27_1_Case_RNASeq_R1.fastq.gz, Sample_PM316_Z27_1_Case_RNASeq_R2.fastq.gz, Sample_PM316_Z8_1_Case_RNASeq_R1.fastq.gz, Sample_PM316_Z8_1_Case_RNASeq_R2.fastq.gz, Sample_PM351_Z6_1_Case_RNASeq_R1.fastq.gz, Sample_PM351_Z6_1_Case_RNASeq_R2.fastq.gz, Sample_PM594_Z7_1_Case_RNASeq_R1.fastq.gz, Sample_PM594_Z7_1_Case_RNASeq_R2.fastq.gz, Sample_PM685_Z1_1_Case_RNASeq_R1.fastq.gz, Sample_PM685_Z1_1_Case_RNASeq_R2.fastq.gz, Sample_PM869_X1_1_Case_RNASeq_R1.fastq.gz, Sample_PM869_X1_1_Case_RNASeq_R2.fastq.gz, Sample_PM876_Z1_1_Case_RNASeq_R1.fastq.gz, Sample_PM876_Z1_1_Case_RNASeq_R2.fastq.gz, Sample_PM937_Z1_1_Case_RNASeq_R1.fastq.gz, Sample_PM937_Z1_1_Case_RNASeq_R2.fastq.gz, Sample_PM976_X1_1_Case_RNASeq_R1.fastq.gz, Sample_PM976_X1_1_Case_RNASeq_R2.fastq.gz
--runThreadN 8
--outFilterMultimapScoreRange 1
--outFilterMultimapNmax 20
--outFilterMismatchNmax 10
--alignIntronMax 500000
--alignMatesGapMax 1000000
--sjdbScore 2
--alignSJDBoverhangMin 1
--genomeLoad NoSharedMemory
--readFilesCommand zcat
--outFilterMatchNminOverLread 0.33
--outFilterScoreMinOverLread 0.33
--sjdbOverhang 100
--outSAMstrandField intronMotif
--outSAMtype None
--outSAMmode None


### Step 3: Intermediate Index Generation.

STAR
--runMode genomeGenerate
--genomeDir /work/jaalab/mskaro1/rcc
--genomeFastaFiles GRCh38.d1.vd1.fa.tar.gz
--sjdbOverhang 100
--runThreadN 8
--sjdbFileChrStartEnd <SJ.out.tab from previous step>



### Step 4: Alignment 2nd Pass.

STAR
--genomeDir <output_path from previous step>
--readFilesIn PM1010_Z3_1_Case_RNASeq_S15_L003_R1.fastq.gz, PM1010_Z3_1_Case_RNASeq_S15_L003_R2.fastq.gz, PM1116_Z3_1_Case_RNASeq_S46_L006_R1_001.fastq.gz, PM1116_Z3_1_Case_RNASeq_S46_L006_R2_001.fastq.gz, PM1220_Z5_1_Case_RNASeq_S17_L003_R1_001.fastq.gz, PM1220_Z5_1_Case_RNASeq_S17_L003_R2_001.fastq.gz, PM1306_Z3_1_Case_RNASeq_S21_L003_R1_001.fastq.gz, PM1306_Z3_1_Case_RNASeq_S21_L003_R2_001.fastq.gz, PM1321_Z1_1_Case_RNASeq_S23_L003_R1_001.fastq.gz, PM1321_Z1_1_Case_RNASeq_S23_L003_R2_001.fastq.gz, PM492_Z5_1_Case_RNASeq_S6_L004_R1_001.fastq.gz, PM492_Z5_1_Case_RNASeq_S6_L004_R2_001.fastq.gz, PM492_Z6_1_Case_RNASeq_S7_L004_R1_001.fastq.gz, PM492_Z6_1_Case_RNASeq_S7_L004_R2_001.fastq.gz, PM492_Z7_1_Case_RNASeq_S8_L004_R1_001.fastq.gz, PM492_Z7_1_Case_RNASeq_S8_L004_R2_001.fastq.gz, Sample_PM1035_X1_1_Case_RNASeq_R1.fastq.gz, Sample_PM1035_X1_1_Case_RNASeq_R2.fastq.gz, Sample_PM1041_Z2_1_Case_RNASeq_R1.fastq.gz, Sample_PM1041_Z2_1_Case_RNASeq_R2.fastq.gz, Sample_PM157_Z7_1_RNA_R1.fastq.gz, Sample_PM157_Z7_1_RNA_R2.fastq.gz, Sample_PM215_Z4_1_Case_RNASeq_R1.fastq.gz, Sample_PM215_Z4_1_Case_RNASeq_R2.fastq.gz, Sample_PM255_Z1_1_Case_RNASeq_R1.fastq.gz, Sample_PM255_Z1_1_Case_RNASeq_R2.fastq.gz, Sample_PM269_Z1_1_Case_RNASeq_R1.fastq.gz, Sample_PM269_Z1_1_Case_RNASeq_R2.fastq.gz, Sample_PM316_Z17_1_Case_RNASeq_R1.fastq.gz, Sample_PM316_Z17_1_Case_RNASeq_R2.fastq.gz, Sample_PM316_Z18_1_Case_RNASeq_R1.fastq.gz, Sample_PM316_Z18_1_Case_RNASeq_R2.fastq.gz, Sample_PM316_Z19_1_Case_RNASeq_R1.fastq.gz, Sample_PM316_Z19_1_Case_RNASeq_R2.fastq.gz, Sample_PM316_Z22_1_Case_RNASeq_R1.fastq.gz, Sample_PM316_Z22_1_Case_RNASeq_R2.fastq.gz, Sample_PM316_Z27_1_Case_RNASeq_R1.fastq.gz, Sample_PM316_Z27_1_Case_RNASeq_R2.fastq.gz, Sample_PM316_Z8_1_Case_RNASeq_R1.fastq.gz, Sample_PM316_Z8_1_Case_RNASeq_R2.fastq.gz, Sample_PM351_Z6_1_Case_RNASeq_R1.fastq.gz, Sample_PM351_Z6_1_Case_RNASeq_R2.fastq.gz, Sample_PM594_Z7_1_Case_RNASeq_R1.fastq.gz, Sample_PM594_Z7_1_Case_RNASeq_R2.fastq.gz, Sample_PM685_Z1_1_Case_RNASeq_R1.fastq.gz, Sample_PM685_Z1_1_Case_RNASeq_R2.fastq.gz, Sample_PM869_X1_1_Case_RNASeq_R1.fastq.gz, Sample_PM869_X1_1_Case_RNASeq_R2.fastq.gz, Sample_PM876_Z1_1_Case_RNASeq_R1.fastq.gz, Sample_PM876_Z1_1_Case_RNASeq_R2.fastq.gz, Sample_PM937_Z1_1_Case_RNASeq_R1.fastq.gz, Sample_PM937_Z1_1_Case_RNASeq_R2.fastq.gz, Sample_PM976_X1_1_Case_RNASeq_R1.fastq.gz, Sample_PM976_X1_1_Case_RNASeq_R2.fastq.gz
--runThreadN 8
--outFilterMultimapScoreRange 1
--outFilterMultimapNmax 20
--outFilterMismatchNmax 10
--alignIntronMax 500000
--alignMatesGapMax 1000000
--sjdbScore 2
--alignSJDBoverhangMin 1
--genomeLoad NoSharedMemory
--limitBAMsortRAM 0
--readFilesCommand zcat
--outFilterMatchNminOverLread 0.33
--outFilterScoreMinOverLread 0.33
--sjdbOverhang 100
--outSAMstrandField intronMotif
--outSAMattributes NH HI NM MD AS XS
--outSAMunmapped Within
--outSAMtype BAM SortedByCoordinate
--outSAMheaderHD @HD VN:1.4
--outSAMattrRGline <formatted RG line provided by wrapper> # gotta figure this one out?


# Author: Michael Skaro
# Date: 8/3/2020
# Purpose: Use the star module to complete alignment of the RNA sequencing data 
# Done


