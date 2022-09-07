#!/bin/bash

# targeted_seq_batch_script.sh
# Description: Batch script for targeted seq. data processing

# All paths are absolute for now
# Important main paths to set
PROJECT_DIR='cn_signatures_shallowWGS'
TMP_PATH=/projects/molonc/scratch/madouglas/cn_signatures
INPUT_PATH=/projects/molonc/scratch/madouglas/transfers/targeted_panel_sequencing/all
OUTPUT_PATH=${TMP_PATH}/targeted_panel_seq
REF_GENOME=/projects/molonc/huntsman_lab/madouglas/bin/reference_genomes/dlp_refdata/human/GRCh37-lite.fa
MUTECT_RESOURCES=/projects/molonc/huntsman_lab/share/hg19_genome/forMutect2
FUNCOTATOR_RESOURCES=/projects/molonc/aparicio_lab/ezaikova/ref/hg38/funcotator_dataSources.v1.7.20200521s

SAMPLE_DIRS=$(ls $INPUT_PATH)

for i in $SAMPLE_DIRS
do
	SAMPLE=${i}

	RUNSCRIPTS=${TMP_PATH}/runscripts/${SAMPLE}
	FASTQC_RAW=${OUTPUT_PATH}/fastQC/${SAMPLE}/raw
	FASTQC_PROC=/projects/molonc/huntsman_lab/madouglas/cn_signatures_shallowWGS/targeted_panel_seq/fastQC/${SAMPLE}/post_alignment
	LOGS=${OUTPUT_PATH}/logs/${SAMPLE}
	ALIGNED=/projects/molonc/scratch/madouglas/cn_signatures/targeted_panel_seq/aligned/${SAMPLE}
	FASTQ_FILES=(${INPUT_PATH}/${i}/*)

	mkdir -p /projects/molonc/scratch/madouglas/cn_signatures/targeted_panel_seq/aligned/${SAMPLE}/mutect2
	# mkdir -p $ALIGNED $RUNSCRIPTS

	echo "#!/bin/bash
#SBATCH --partition=upgrade
#SBATCH --time=32:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=${LOGS}/${SAMPLE}.std_out.2.txt
#SBATCH --job-name=targeted_panel_seq
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32
#SBATCH --workdir=/projects/molonc/archive/madouglas

# fastqc -t 4 ${INPUT_PATH}/$i/*.gz -o $FASTQC_RAW/
# source /home/madouglas/.bash
# conda activate py37
# fastp -i ${FASTQ_FILES[0]} -I ${FASTQ_FILES[1]} -o $INPUT_PATH/${SAMPLE}/${SAMPLE}.R1.fastp.fastq.gz -O $INPUT_PATH/${SAMPLE}/${SAMPLE}.R2.fastp.fastq.gz -U --umi_loc=read2 --umi_len=12 -F 11 -w 16 -h $ALIGNED/$SAMPLE.fastp.html
#
# echo "${INPUT_PATH}/${SAMPLE}/${SAMPLE}.R1.fastp.fastq.gz"
# gunzip -c ${INPUT_PATH}/${SAMPLE}/${SAMPLE}.R1.fastp.fastq.gz | awk '{if(NR%4==2) print length($1)}' > $ALIGNED/$SAMPLE.input.readslength.txt
# textHistogram -binSize=10 $ALIGNED/$SAMPLE.input.readslength.txt
# echo "${INPUT_PATH}/${SAMPLE}/${SAMPLE}.R2.fastp.fastq.gz"
# gunzip -c ${INPUT_PATH}/${SAMPLE}/${SAMPLE}.R2.fastp.fastq.gz | awk '{if(NR%4==2) print length($1)}' > $ALIGNED/$SAMPLE.input.readslength.txt
# textHistogram -binSize=10 $ALIGNED/$SAMPLE.input.readslength.txt
#
# ########
# ########
# #### Paired-End Alignment and Pre-processing
#
# bwa-mem2 mem -M -t 24 $REF_GENOME $INPUT_PATH/${SAMPLE}/${SAMPLE}.R1.fastp.fastq.gz $INPUT_PATH/${SAMPLE}/${SAMPLE}.R2.fastp.fastq.gz > $ALIGNED/$SAMPLE.pe.bwa.sam
# samtools view -h -b -S $ALIGNED/$SAMPLE.pe.bwa.sam > $ALIGNED/$SAMPLE.pe.bwa.bam
# java -jar /projects/molonc/huntsman_lab/madouglas/bin/picard.jar AddOrReplaceReadGroups \
# 	I=$ALIGNED/$SAMPLE.pe.bwa.bam \
# 	O=$ALIGNED/$SAMPLE.pe.bwa.RG.bam \
# 	SORT_ORDER=queryname RGID=1 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 \
# 	RGSM=$SAMPLE \
# 	RGCN=BC-GSC
# java -jar /projects/molonc/huntsman_lab/madouglas/bin/picard.jar FixMateInformation \
#     I=$ALIGNED/$SAMPLE.pe.bwa.RG.bam \
#     O=$ALIGNED/$SAMPLE.pe.bwa.fixmate.bam \
#     ADD_MATE_CIGAR=true
# samtools sort -@ 24 -o $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.bam $ALIGNED/$SAMPLE.pe.bwa.fixmate.bam
# samtools index -@ 24 $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.bam

# java -jar /projects/molonc/huntsman_lab/madouglas/bin/picard.jar UmiAwareMarkDuplicatesWithMateCigar \
#     I=$ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.bam \
#     O=$ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.bam \
#     M=$ALIGNED/$SAMPLE.pe.marked_dup_metrics.txt \
#     UMI_METRICS=$ALIGNED/$SAMPLE.pe.dup_rm.umi_metrics.txt
#
# gatk BaseRecalibrator \
#     -I $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.bam \
#     -R $REF_GENOME \
#     --known-sites $MUTECT_RESOURCES/1000G_phase1.snps.high_confidence.b37.vcf.gz \
#     --known-sites $MUTECT_RESOURCES/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
#     --known-sites $MUTECT_RESOURCES/dbsnp_138.b37.excluding_sites_after_129.vcf.gz \
#     -O $ALIGNED/$SAMPLE.recal_data.table
#
# gatk ApplyBQSR \
#     -R $REF_GENOME \
#     -I $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.bam \
#     --bqsr-recal-file $ALIGNED/$SAMPLE.recal_data.table \
#     -O $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.recal.bam
#
# gatk BaseRecalibrator \
#     -I $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.recal.bam \
#     -R $REF_GENOME \
#     --known-sites $MUTECT_RESOURCES/1000G_phase1.snps.high_confidence.b37.vcf.gz \
#     --known-sites $MUTECT_RESOURCES/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
#     --known-sites $MUTECT_RESOURCES/dbsnp_138.b37.excluding_sites_after_129.vcf.gz \
#     -O $ALIGNED/$SAMPLE.recal_data.post.table
#
# gatk AnalyzeCovariates \
#     -before $ALIGNED/$SAMPLE.recal_data.table \
#     -after $ALIGNED/$SAMPLE.recal_data.post.table \
#     -csv $ALIGNED/$SAMPLE.AnalyzeCovariates.csv
#
# Rscript_4.1 /projects/molonc/huntsman_lab/madouglas/bin/BQSR.R \
#     $ALIGNED/$SAMPLE.AnalyzeCovariates.csv $ALIGNED/$SAMPLE.recal_data.table $ALIGNED/$SAMPLE.AnalyzeCovariates.pdf
#
# samtools index -@ 24 $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.recal.bam
#
# java -jar /projects/molonc/huntsman_lab/madouglas/bin/picard.jar ValidateSamFile \
#     IGNORE_WARNINGS=true MODE=VERBOSE \
#     I=$ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.recal.bam

##################
##################
######### QIAseq bams
#########
samtools view -H $ALIGNED/$SAMPLE.qiaseq.bam > $ALIGNED/$SAMPLE.qiaseq.old.header.sam

Rscript_4.1 /projects/molonc/huntsman_lab/madouglas/cn_signatures_shallowWGS/scripts/reheadering_qiaseq_bams.R -i $ALIGNED/$SAMPLE.qiaseq.old.header.sam -o $ALIGNED/$SAMPLE.qiaseq.new.header.sam

samtools reheader $ALIGNED/$SAMPLE.qiaseq.new.header.sam $ALIGNED/$SAMPLE.qiaseq.bam > $ALIGNED/$SAMPLE.qiaseq.reheadered.bam

java -jar /projects/molonc/huntsman_lab/madouglas/bin/picard.jar AddOrReplaceReadGroups \
		I=$ALIGNED/$SAMPLE.qiaseq.reheadered.bam \
		O=$ALIGNED/$SAMPLE.qiaseq.reheadered.RG.bam \
		SORT_ORDER=queryname RGID=1 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 \
		RGSM=$SAMPLE \
		RGCN=BC-GSC
java -jar /projects/molonc/huntsman_lab/madouglas/bin/picard.jar FixMateInformation \
    I=$ALIGNED/$SAMPLE.qiaseq.reheadered.RG.bam \
    O=$ALIGNED/$SAMPLE.qiaseq.reheadered.RG.fixmate.bam \
    ADD_MATE_CIGAR=true
samtools sort -@ 24 -o $ALIGNED/$SAMPLE.qiaseq.reheadered.RG.fixmate.sorted.bam $ALIGNED/$SAMPLE.qiaseq.reheadered.RG.fixmate.bam
samtools index -@ 24 $ALIGNED/$SAMPLE.qiaseq.reheadered.RG.fixmate.sorted.bam

gatk BaseRecalibrator \
    -I $ALIGNED/$SAMPLE.qiaseq.reheadered.RG.fixmate.sorted.bam \
    -R $REF_GENOME \
    --known-sites $MUTECT_RESOURCES/1000G_phase1.snps.high_confidence.b37.vcf.gz \
    --known-sites $MUTECT_RESOURCES/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
    --known-sites $MUTECT_RESOURCES/dbsnp_138.b37.excluding_sites_after_129.vcf.gz \
    -O $ALIGNED/$SAMPLE.qiaseq.recal_data.table

gatk ApplyBQSR \
    -R $REF_GENOME \
    -I $ALIGNED/$SAMPLE.qiaseq.reheadered.RG.fixmate.sorted.bam \
    --bqsr-recal-file $ALIGNED/$SAMPLE.qiaseq.recal_data.table \
    -O $ALIGNED/$SAMPLE.qiaseq.reheadered.RG.fixmate.sorted.recal.bam

gatk BaseRecalibrator \
    -I $ALIGNED/$SAMPLE.qiaseq.reheadered.RG.fixmate.sorted.recal.bam \
    -R $REF_GENOME \
    --known-sites $MUTECT_RESOURCES/1000G_phase1.snps.high_confidence.b37.vcf.gz \
    --known-sites $MUTECT_RESOURCES/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
    --known-sites $MUTECT_RESOURCES/dbsnp_138.b37.excluding_sites_after_129.vcf.gz \
    -O $ALIGNED/$SAMPLE.qiaseq.recal_data.post.table

gatk AnalyzeCovariates \
    -before $ALIGNED/$SAMPLE.qiaseq.recal_data.table \
    -after $ALIGNED/$SAMPLE.qiaseq.recal_data.post.table \
    -csv $ALIGNED/$SAMPLE.qiaseq.AnalyzeCovariates.csv

Rscript_4.1 /projects/molonc/huntsman_lab/madouglas/bin/BQSR.R \
    $ALIGNED/$SAMPLE.qiaseq.AnalyzeCovariates.csv $ALIGNED/$SAMPLE.qiaseq.recal_data.table $ALIGNED/$SAMPLE.qiaseq.AnalyzeCovariates.pdf

samtools index -@ 24 $ALIGNED/$SAMPLE.qiaseq.reheadered.RG.fixmate.sorted.recal.bam

java -jar /projects/molonc/huntsman_lab/madouglas/bin/picard.jar ValidateSamFile \
    IGNORE_WARNINGS=true MODE=VERBOSE \
    I=$ALIGNED/$SAMPLE.qiaseq.reheadered.RG.fixmate.sorted.recal.bam

rm $ALIGNED/$SAMPLE.qiaseq.reheadered.bam
rm $ALIGNED/$SAMPLE.qiaseq.reheadered.RG.bam
rm $ALIGNED/$SAMPLE.qiaseq.reheadered.RG.fixmate.bam
rm $ALIGNED/$SAMPLE.qiaseq.reheadered.RG.fixmate.sorted.bam.bai
rm $ALIGNED/$SAMPLE.qiaseq.reheadered.RG.fixmate.sorted.bam
rm $ALIGNED/$SAMPLE.qiaseq.reheadered.RG.fixmate.sorted.recal.bai

###
### Somatic short mutation detection using Mutect2

gatk Mutect2 \
    -R $REF_GENOME \
    -I $ALIGNED/$SAMPLE.qiaseq.reheadered.RG.fixmate.sorted.recal.bam \
    --germline-resource $MUTECT_RESOURCES/af-only-gnomad.raw.sites.vcf \
    --panel-of-normals $MUTECT_RESOURCES/Mutect2-WGS-panel-b37.vcf \
    --tmp-dir $TMP_PATH \
    --f1r2-tar-gz $ALIGNED/mutect2/$SAMPLE.qiaseq.f1r2.tar.gz \
    -O $ALIGNED/mutect2/$SAMPLE.qiaseq.vcf.gz

gatk LearnReadOrientationModel -I $ALIGNED/mutect2/$SAMPLE.qiaseq.f1r2.tar.gz \
    -O $ALIGNED/mutect2/$SAMPLE.qiaseq.read-orientation-model.tar.gz \
    --num-em-iterations 50

gatk GetPileupSummaries \
    -I $ALIGNED/$SAMPLE.qiaseq.reheadered.RG.fixmate.sorted.recal.bam \
    -V $MUTECT_RESOURCES/small_exac_common_3.vcf \
    -L $MUTECT_RESOURCES/small_exac_common_3.vcf \
    -O $ALIGNED/mutect2/$SAMPLE.qiaseq.getpileupsummaries.table

gatk CalculateContamination \
    -I $ALIGNED/mutect2/$SAMPLE.qiaseq.getpileupsummaries.table \
    -tumor-segmentation $ALIGNED/mutect2/$SAMPLE.qiaseq.segments.table \
    -O $ALIGNED/mutect2/$SAMPLE.qiaseq.contamination.table

gatk FilterMutectCalls -R $REF_GENOME \
    -V $ALIGNED/mutect2/$SAMPLE.qiaseq.vcf.gz \
    --tumor-segmentation $ALIGNED/mutect2/$SAMPLE.qiaseq.segments.table \
    --contamination-table $ALIGNED/mutect2/$SAMPLE.qiaseq.contamination.table \
    --ob-priors $ALIGNED/mutect2/$SAMPLE.qiaseq.read-orientation-model.tar.gz \
    -O $ALIGNED/mutect2/$SAMPLE.qiaseq.filtered.vcf

# gatk Funcotator \
#     --variant $ALIGNED/$SAMPLE.qiaseq.filtered.vcf \
#     --reference $REF_GENOME \
#     --ref-version hg19 \
#     --data-sources-path $FUNCOTATOR_RESOURCES \
#     --output $ALIGNED/$SAMPLE.qiaseq.filtered.funcotated.vcf \
#     --output-file-format VCF \
#     --force-b37-to-hg19-reference-contig-conversion

gatk VariantsToTable \
    -V $ALIGNED/mutect2/$SAMPLE.qiaseq.filtered.vcf \
    -F CHROM -F POS -F REF -F ALT -F FILTER \
    -F GERMQ -F MBQ -F MFRL -F MMQ -F ROQ -F TLOD \
    -GF GT -GF AD -GF AF -GF DP -GF SB \
    -O $ALIGNED/mutect2/$SAMPLE.qiaseq.filtered.tsv

rm $ALIGNED/mutect2/$SAMPLE.qiaseq.vcf.gz
rm $ALIGNED/mutect2/$SAMPLE.qiaseq.vcf.gz.tbi
rm $ALIGNED/mutect2/$SAMPLE.qiaseq.vcf.gz.stats
##################
##################
##################


#
#
# rm $ALIGNED/$SAMPLE.pe.bwa.sam
# rm $ALIGNED/$SAMPLE.pe.bwa.RG.bam
# rm $ALIGNED/$SAMPLE.pe.bwa.fixmate.bam
# rm $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.bam
# rm $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.bam
#
# fastqc $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.recal.bam -o $FASTQC_PROC/
#
# ########
# ########
# #### Somatic short mutation detection using Mutect2
#
# gatk Mutect2 \
#     -R $REF_GENOME \
#     -I $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.recal.bam \
#     --germline-resource $MUTECT_RESOURCES/af-only-gnomad.raw.sites.vcf \
#     --panel-of-normals $MUTECT_RESOURCES/Mutect2-WGS-panel-b37.vcf \
#     --tmp-dir $TMP_PATH \
#     --f1r2-tar-gz $ALIGNED/$SAMPLE.f1r2.tar.gz \
#     -O $ALIGNED/$SAMPLE.vcf.gz
#
# gatk LearnReadOrientationModel -I $ALIGNED/$SAMPLE.f1r2.tar.gz \
#     -O $ALIGNED/$SAMPLE.read-orientation-model.tar.gz \
#     --num-em-iterations 50
#
# gatk GetPileupSummaries \
#     -I $ALIGNED/$SAMPLE.pe.bwa.sorted.fixmate.dup_rm.recal.bam \
#     -V $MUTECT_RESOURCES/small_exac_common_3.vcf \
#     -L $MUTECT_RESOURCES/small_exac_common_3.vcf \
#     -O $ALIGNED/$SAMPLE.getpileupsummaries.table
#
# gatk CalculateContamination \
#     -I $ALIGNED/$SAMPLE.getpileupsummaries.table \
#     -tumor-segmentation $ALIGNED/$SAMPLE.segments.table \
#     -O $ALIGNED/$SAMPLE.contamination.table
#
# gatk FilterMutectCalls -R $REF_GENOME \
#     -V $ALIGNED/$SAMPLE.vcf.gz \
#     --tumor-segmentation $ALIGNED/$SAMPLE.segments.table \
#     --contamination-table $ALIGNED/$SAMPLE.contamination.table \
#     --ob-priors $ALIGNED/$SAMPLE.read-orientation-model.tar.gz \
#     -O $ALIGNED/$SAMPLE.filtered.vcf
#
# gatk Funcotator \
#     --variant $ALIGNED/$SAMPLE.filtered.vcf \
#     --reference $REF_GENOME \
#     --ref-version hg19 \
#     --data-sources-path $FUNCOTATOR_RESOURCES \
#     --output $ALIGNED/$SAMPLE.filtered.funcotated.vcf \
#     --output-file-format VCF

" > $RUNSCRIPTS/${SAMPLE}.targetedseq.run.job
	echo "Submitting $RUNSCRIPTS/${SAMPLE}.targetedseq.run.job to the cluster"
	sbatch $RUNSCRIPTS/${SAMPLE}.targetedseq.run.job
done
