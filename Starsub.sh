#PBS -S /bin/bash
#PBS -N mfs_star
#PBS -q batch
#PBS -l nodes=1:ppn=8
#PBS -l walltime=64:00:00
#PBS -l mem=100gb
#PBS -M mfs61958@uga.edu 
#PBS -m ae

cd $PBS_O_WORKDIR

echo
echo "Job ID: $PBS_JOBID"
echo "Queue:  $PBS_QUEUE"
echo "Cores:  $PBS_NP"
echo "Nodes:  $(cat $PBS_NODEFILE | sort -u | tr '\n' ' ')"
echo

# Step 1: Load STAR 

ml STAR/2.6.1c-foss-2016b

### Step 1: Index Generation.

STAR --runMode genomeGenerate --genomeDir /work/jaalab/mskaro1/rcc/fq/Star_indicies/ --genomeFastaFiles GRCh38.d1.vd1.fa --sjdbOverhang 74 --runThreadN 8 

### Step 2: Alignment 1st Pass with pre-made indicies

STAR --genomeDir /work/jaalab/mskaro1/rcc/fq/Star_indicies/ --readFilesManifest /work/jaalab/mskaro1/rcc/fq/manifest.tsv  --runThreadN 8 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --readFilesCommand zcat --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbGTFfile /work/jaalab/mskaro1/rcc/fq/gencode.v22.annotation.gtf --sjdbOverhang 74 --outSAMstrandField intronMotif --outSAMtype None --outSAMmode None

### Step 3: Intermediate Index Generation.

STAR --runMode genomeGenerate --genomeDir /work/jaalab/mskaro1/rcc/fq/Star_indicies/ --genomeFastaFiles GRCh38.d1.vd1.fa --sjdbOverhang 74 --runThreadN 8 --sjdbFileChrStartEnd /work/jaalab/mskaro1/rcc/fq/Star_indicies/sjdbList.out.tab

### Step 4: Alignment 2nd Pass.
STAR --genomeDir /work/jaalab/mskaro1/rcc/fq/Star_indicies/ ---readFilesManifest /work/jaalab/mskaro1/rcc/fq/manifest.tsv --runThreadN 8 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --readFilesCommand zcat --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbGTFfile /work/jaalab/mskaro1/rcc/fq/gencode.v22.annotation.gtf --sjdbOverhang 74 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 --quantMode GeneCounts








