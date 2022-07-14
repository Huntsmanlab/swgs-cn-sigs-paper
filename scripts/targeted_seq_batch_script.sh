#!/bin/bash

# targeted_seq_batch_script.sh
# Description: Batch script for targeted seq. data processing

# All paths are absolute for now
# Important main paths to set
PROJECT_DIR='cn_signatures_shallowWGS'
TMP_PATH=/projects/molonc/scratch/madouglas/cn_signatures
INPUT_PATH=/projects/molonc/scratch/madouglas/transfers/targeted_panel_sequencing/20220519/fastqs
OUTPUT_PATH=${TMP_PATH}/targeted_panel_seq
REF_GENOME=/projects/molonc/huntsman_lab/madouglas/bin/reference_genomes/dlp_refdata/human/GRCh37-lite.fa
MUTECT_RESOURCES=/projects/molonc/huntsman_lab/share/hg19_genome/forMutect2
FUNCOTATOR_RESOURCES=/projects/molonc/aparicio_lab/ezaikova/ref/hg38/funcotator_dataSources.v1.7.20200521s

SAMPLE_DIRS=$(ls $INPUT_PATH)

for i in $SAMPLE_DIRS
do
	SAMPLE=$(ls ${INPUT_PATH}/$i | egrep -Eo '^[^_]+' | uniq)

	RUNSCRIPTS=${TMP_PATH}/runscripts/${SAMPLE}
	FASTQC_RAW=${OUTPUT_PATH}/fastQC/${SAMPLE}/raw
	FASTQC_PROC=/projects/molonc/huntsman_lab/madouglas/cn_signatures_shallowWGS/targetted_panel_seq/fastQC/${SAMPLE}/post_alignment
	LOGS=${OUTPUT_PATH}/logs/${SAMPLE}
	ALIGNED=/projects/molonc/huntsman_lab/madouglas/cn_signatures_shallowWGS/targetted_panel_seq/aligned/${SAMPLE}
	FASTQ_FILES=(${INPUT_PATH}/$i/*)
	
	mkdir -p $LOGS $FASTQC_RAW $FASTQC_PROC
	mkdir -p $ALIGNED $RUNSCRIPTS

	echo "#!/bin/bash
#SBATCH --partition=upgrade
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=${LOGS}/${SAMPLE}.std_out.2.txt
#SBATCH --job-name=targeted_panel_seq
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32
#SBATCH --workdir=/projects/molonc/archive/madouglas

# fastqc -t 4 ${INPUT_PATH}/$i/*.gz -o $FASTQC_RAW/

# echo "${FASTQ_FILES[0]}"
# gunzip -c ${FASTQ_FILES[0]} | awk '{if(NR%4==2) print length($1)}' > $ALIGNED/$SAMPLE.input.readslength.txt
# textHistogram -binSize=10 $ALIGNED/$SAMPLE.input.readslength.txt
# echo "${FASTQ_FILES[1]}"
# gunzip -c ${FASTQ_FILES[1]} | awk '{if(NR%4==2) print length($1)}' > $ALIGNED/$SAMPLE.input.readslength.txt
# textHistogram -binSize=10 $ALIGNED/$SAMPLE.input.readslength.txt

########
########
#### Paired-End Alignment and Pre-processing

# bwa-mem2 mem -M -t 24 $REF_GENOME ${FASTQ_FILES[0]} ${FASTQ_FILES[1]} > $ALIGNED/$SAMPLE.pe.bwa.sam
# samtools view -h -b -S $ALIGNED/$SAMPLE.pe.bwa.sam > $ALIGNED/$SAMPLE.pe.bwa.bam
java -jar /projects/molonc/huntsman_lab/madouglas/bin/picard.jar AddOrReplaceReadGroups \
	I=$ALIGNED/$SAMPLE.pe.bwa.bam \
	O=$ALIGNED/$SAMPLE.pe.bwa.RG.bam \
	SORT_ORDER=queryname RGID=1 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 \
	RGSM=$SAMPLE \
	RGCN=BC-GSC
java -jar /projects/molonc/huntsman_lab/madouglas/bin/picard.jar FixMateInformation \
    I=$ALIGNED/$SAMPLE.pe.bwa.RG.bam \
    O=$ALIGNED/$SAMPLE.pe.bwa.fixmate.bam \
    ADD_MATE_CIGAR=true
samtools sort -@ 24 -o $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.bam $ALIGNED/$SAMPLE.pe.bwa.fixmate.bam
samtools index -@ 24 $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.bam
# java -jar /projects/molonc/huntsman_lab/madouglas/bin/picard.jar MarkDuplicates \
# 		I=$ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.bam \
# 		O=$ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.bam \
# 		M=$ALIGNED/$SAMPLE.pe.marked_dup_metrics.txt

java -jar /projects/molonc/huntsman_lab/madouglas/bin/picard.jar UmiAwareMarkDuplicatesWithMateCigar \
    I=$ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.bam \
    O=$ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.bam \
    M=$ALIGNED/$SAMPLE.pe.marked_dup_metrics.txt \
    UMI_METRICS=$ALIGNED/$SAMPLE.pe.dup_rm.umi_metrics.txt

gatk BaseRecalibrator \
    -I $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.bam \
    -R $REF_GENOME \
    --known-sites $MUTECT_RESOURCES/1000G_phase1.snps.high_confidence.b37.vcf.gz \
    --known-sites $MUTECT_RESOURCES/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
    --known-sites $MUTECT_RESOURCES/dbsnp_138.b37.excluding_sites_after_129.vcf.gz \
    -O $ALIGNED/$SAMPLE.recal_data.table
   
gatk ApplyBQSR \
    -R $REF_GENOME \
    -I $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.bam \
    --bqsr-recal-file $ALIGNED/$SAMPLE.recal_data.table \
    -O $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.recal.bam
    
gatk BaseRecalibrator \
    -I $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.recal.bam \
    -R $REF_GENOME \
    --known-sites $MUTECT_RESOURCES/1000G_phase1.snps.high_confidence.b37.vcf.gz \
    --known-sites $MUTECT_RESOURCES/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
    --known-sites $MUTECT_RESOURCES/dbsnp_138.b37.excluding_sites_after_129.vcf.gz \
    -O $ALIGNED/$SAMPLE.recal_data.post.table

gatk AnalyzeCovariates \
    -before $ALIGNED/$SAMPLE.recal_data.table \
    -after $ALIGNED/$SAMPLE.recal_data.post.table \
    -csv $ALIGNED/$SAMPLE.AnalyzeCovariates.csv
   
Rscript_4.1 /projects/molonc/huntsman_lab/madouglas/bin/BQSR.R \
    $ALIGNED/$SAMPLE.AnalyzeCovariates.csv $ALIGNED/$SAMPLE.recal_data.table $ALIGNED/$SAMPLE.AnalyzeCovariates.pdf

samtools index -@ 24 $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.recal.bam

java -jar /projects/molonc/huntsman_lab/madouglas/bin/picard.jar ValidateSamFile \
    IGNORE_WARNINGS=true MODE=VERBOSE \
    I=$ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.recal.bam


rm $ALIGNED/$SAMPLE.pe.bwa.sam
rm $ALIGNED/$SAMPLE.pe.bwa.RG.bam
rm $ALIGNED/$SAMPLE.pe.bwa.fixmate.bam
rm $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.bam
rm $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.bam

fastqc $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.recal.bam -o $FASTQC_PROC/

########
########
#### Somatic short mutation detection using Mutect2

gatk Mutect2 \
    -R $REF_GENOME \
    -I $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.recal.bam \
    --germline-resource $MUTECT_RESOURCES/af-only-gnomad.raw.sites.vcf \
    --panel-of-normals $MUTECT_RESOURCES/Mutect2-WGS-panel-b37.vcf \
    --tmp-dir $TMP_PATH \
    --f1r2-tar-gz $ALIGNED/$SAMPLE.f1r2.tar.gz \
    -O $ALIGNED/$SAMPLE.vcf.gz

gatk LearnReadOrientationModel -I $ALIGNED/$SAMPLE.f1r2.tar.gz \
    -O $ALIGNED/$SAMPLE.read-orientation-model.tar.gz \
    --num-em-iterations 50

gatk GetPileupSummaries \
    -I $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.recal.bam \
    -V $MUTECT_RESOURCES/small_exac_common_3.vcf \
    -L $MUTECT_RESOURCES/small_exac_common_3.vcf \
    -O $ALIGNED/$SAMPLE.getpileupsummaries.table

gatk CalculateContamination \
    -I $ALIGNED/$SAMPLE.getpileupsummaries.table \
    -tumor-segmentation $ALIGNED/$SAMPLE.segments.table \
    -O $ALIGNED/$SAMPLE.contamination.table

gatk FilterMutectCalls -R $REF_GENOME \
    -V $ALIGNED/$SAMPLE.vcf.gz \
    --tumor-segmentation $ALIGNED/$SAMPLE.segments.table \
    --contamination-table $ALIGNED/$SAMPLE.contamination.table \
    --ob-priors $ALIGNED/$SAMPLE.read-orientation-model.tar.gz \
    -O $ALIGNED/$SAMPLE.filtered.vcf

gatk Funcotator \
    --variant $ALIGNED/$SAMPLE.filtered.vcf \
    --reference $REF_GENOME \
    --ref-version hg19 \
    --data-sources-path $FUNCOTATOR_RESOURCES \
    --output $ALIGNED/$SAMPLE.filtered.funcotated.vcf \
    --output-file-format VCF    

" > $RUNSCRIPTS/${SAMPLE}.targetedseq.run.job
	echo "Submitting $RUNSCRIPTS/${SAMPLE}.targetedseq.run.job to the cluster"
	# sbatch $RUNSCRIPTS/${SAMPLE}.targetedseq.run.job
done
