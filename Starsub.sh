#PBS -S /bin/bash
#PBS -N mfs_star_all_samples_1_run
#PBS -q highmem_q
#PBS -l nodes=1:ppn=16
#PBS -l walltime=72:00:00
#PBS -l mem=500gb
#PBS -M mfs61958@uga.edu 
#PBS -m ae

cd $PBS_O_WORKDIR

echo
echo "Job ID: $PBS_JOBID"
echo "Queue:  $PBS_QUEUE"
echo "Cores:  $PBS_NP"
echo "Nodes:  $(cat $PBS_NODEFILE | sort -u | tr '\n' ' ')"
echo


# Load STAR, RSEM, picard,

ml STAR/2.6.1c-foss-2016b
module load RSEM/1.3.0-foss-2016b
module load picard/2.16.0-Java-1.8.0_144

### Step 1: Index Generation.

STAR --runMode genomeGenerate \
--genomeDir /work/jaalab/mskaro1/rcc/fq/Star_indicies/ \
--genomeFastaFiles GRCh38.d1.vd1.fa \
--sjdbGTFfile /work/jaalab/mskaro1/rcc/fq/gencode.v22.annotation.gtf \
--sjdbOverhang 74 --runThreadN 16


### Step 1: Two Pass alignment and quantification.
STAR --runMode alignReads --genomeDir /work/jaalab/mskaro1/rcc/fq/Star_indicies/ --twopassMode Basic   --runThreadN 16 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --readFilesCommand zcat --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbGTFfile /work/jaalab/mskaro1/rcc/fq/gencode.v22.annotation.gtf --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 --outSAMattrRGline ID:PM1010_Z3_1 , ID:PM1116_Z3_1 , ID:PM1220_Z5_1 , ID:PM1306_Z3_1 , ID:PM1321_Z1_1 , ID:PM492_Z5_1 , ID:PM492_Z6_1 , ID:PM492_Z7_1 , ID:PM1035_X1_1 , ID:PM1041_Z2_1 , ID:PM157_Z7_1 , ID:PM215_Z4_1 , ID:PM255_Z1_1 , ID:PM269_Z1_1 , ID:PM316_Z17_1 , ID:PM316_Z18_1 , ID:PM316_Z19_1 , ID:PM316_Z22_1 , ID:PM316_Z27_1 , ID:PM316_Z8_1 , ID:PM351_Z6_1 , ID:PM594_Z7_1 , ID:PM685_Z1_1 , ID:PM869_X1_1 , ID:PM876_Z1_1 , ID:PM937_Z1_1 , ID:PM976_X1_1 --quantMode TranscriptomeSAM GeneCounts



### Step 2: Two Pass alignment and quantification.
STAR --runMode alignReads \
--runThreadN 16 \
--genomeDir /work/jaalab/mskaro1/rcc/fq/Star_indicies/ \
--twopassMode Basic \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outFilterType BySJout \
--outFilterScoreMinOverLread 0.33 \
--outFilterMatchNminOverLread 0.33 \
--limitSjdbInsertNsj 1200000 \
--readFilesIn PM1010_Z3_1_Case_RNASeq_S15_L003_R1.fastq.gz,PM1116_Z3_1_Case_RNASeq_S46_L006_R1_001.fastq.gz,PM1220_Z5_1_Case_RNASeq_S17_L003_R1_001.fastq.gz,PM1306_Z3_1_Case_RNASeq_S21_L003_R1_001.fastq.gz,PM1321_Z1_1_Case_RNASeq_S23_L003_R1_001.fastq.gz,PM492_Z5_1_Case_RNASeq_S6_L004_R1_001.fastq.gz,PM492_Z6_1_Case_RNASeq_S7_L004_R1_001.fastq.gz,PM492_Z7_1_Case_RNASeq_S8_L004_R1_001.fastq.gz,Sample_PM1035_X1_1_Case_RNASeq_R1.fastq.gz,Sample_PM1041_Z2_1_Case_RNASeq_R1.fastq.gz,Sample_PM157_Z7_1_RNA_R1.fastq.gz,Sample_PM215_Z4_1_Case_RNASeq_R1.fastq.gz,Sample_PM255_Z1_1_Case_RNASeq_R1.fastq.gz,Sample_PM269_Z1_1_Case_RNASeq_R1.fastq.gz,Sample_PM316_Z17_1_Case_RNASeq_R1.fastq.gz,Sample_PM316_Z18_1_Case_RNASeq_R1.fastq.gz,Sample_PM316_Z19_1_Case_RNASeq_R1.fastq.gz,Sample_PM316_Z22_1_Case_RNASeq_R1.fastq.gz,Sample_PM316_Z27_1_Case_RNASeq_R1.fastq.gz,Sample_PM316_Z8_1_Case_RNASeq_R1.fastq.gz,Sample_PM351_Z6_1_Case_RNASeq_R1.fastq.gz,Sample_PM594_Z7_1_Case_RNASeq_R1.fastq.gz,Sample_PM685_Z1_1_Case_RNASeq_R1.fastq.gz,Sample_PM869_X1_1_Case_RNASeq_R1.fastq.gz,Sample_PM876_Z1_1_Case_RNASeq_R1.fastq.gz,Sample_PM937_Z1_1_Case_RNASeq_R1.fastq.gz,Sample_PM976_X1_1_Case_RNASeq_R1.fastq.gz PM1010_Z3_1_Case_RNASeq_S15_L003_R2.fastq.gz,PM1116_Z3_1_Case_RNASeq_S46_L006_R2_001.fastq.gz,PM1220_Z5_1_Case_RNASeq_S17_L003_R2_001.fastq.gz,PM1306_Z3_1_Case_RNASeq_S21_L003_R2_001.fastq.gz,PM1321_Z1_1_Case_RNASeq_S23_L003_R2_001.fastq.gz,PM492_Z5_1_Case_RNASeq_S6_L004_R2_001.fastq.gz,PM492_Z6_1_Case_RNASeq_S7_L004_R2_001.fastq.gz,PM492_Z7_1_Case_RNASeq_S8_L004_R2_001.fastq.gz,Sample_PM1035_X1_1_Case_RNASeq_R2.fastq.gz,Sample_PM1041_Z2_1_Case_RNASeq_R2.fastq.gz,Sample_PM157_Z7_1_RNA_R2.fastq.gz,Sample_PM215_Z4_1_Case_RNASeq_R2.fastq.gz,Sample_PM255_Z1_1_Case_RNASeq_R2.fastq.gz,Sample_PM269_Z1_1_Case_RNASeq_R2.fastq.gz,Sample_PM316_Z17_1_Case_RNASeq_R2.fastq.gz,Sample_PM316_Z18_1_Case_RNASeq_R2.fastq.gz,Sample_PM316_Z19_1_Case_RNASeq_R2.fastq.gz,Sample_PM316_Z22_1_Case_RNASeq_R2.fastq.gz,Sample_PM316_Z27_1_Case_RNASeq_R2.fastq.gz,Sample_PM316_Z8_1_Case_RNASeq_R2.fastq.gz,Sample_PM351_Z6_1_Case_RNASeq_R2.fastq.gz,Sample_PM594_Z7_1_Case_RNASeq_R2.fastq.gz,Sample_PM685_Z1_1_Case_RNASeq_R2.fastq.gz,Sample_PM869_X1_1_Case_RNASeq_R2.fastq.gz,Sample_PM876_Z1_1_Case_RNASeq_R2.fastq.gz,Sample_PM937_Z1_1_Case_RNASeq_R2.fastq.gz,Sample_PM976_X1_1_Case_RNASeq_R2.fastq.gz \
--readFilesCommand zcat \
--outSAMstrandField intronMotif \
--outFilterIntronMotifs None \
--alignSoftClipAtReferenceEnds Yes \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped None \
--genomeLoad NoSharedMemory \
--chimSegmentMin 15 \
--chimJunctionOverhangMin 15 \
--chimOutType Junctions WithinBAM SoftClip \
--chimMainSegmentMultNmax 1 \
--outSAMattributes NH HI AS nM NM ch \
--outSAMattrRGline ID:PM1010_Z3_1 , ID:PM1116_Z3_1 , ID:PM1220_Z5_1 , ID:PM1306_Z3_1 , ID:PM1321_Z1_1 , ID:PM492_Z5_1 , ID:PM492_Z6_1 , ID:PM492_Z7_1 , ID:PM1035_X1_1 , ID:PM1041_Z2_1 , ID:PM157_Z7_1 , ID:PM215_Z4_1 , ID:PM255_Z1_1 , ID:PM269_Z1_1 , ID:PM316_Z17_1 , ID:PM316_Z18_1 , ID:PM316_Z19_1 , ID:PM316_Z22_1 , ID:PM316_Z27_1 , ID:PM316_Z8_1 , ID:PM351_Z6_1 , ID:PM594_Z7_1 , ID:PM685_Z1_1 , ID:PM869_X1_1 , ID:PM876_Z1_1 , ID:PM937_Z1_1 , ID:PM976_X1_1



