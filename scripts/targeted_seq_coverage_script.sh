#!/bin/bash

# targeted_seq_coverage_script.sh
# Description: Batch script for targeted seq. data processing

# All paths are absolute for now
# Important main paths to set
PROJECT_DIR='cn_signatures_shallowWGS'
TMP_PATH=/projects/molonc/scratch/madouglas/cn_signatures
INPUT_PATH=/projects/molonc/scratch/madouglas/transfers/targeted_panel_sequencing/20220519/fastqs
OUTPUT_PATH=${TMP_PATH}/targeted_panel_seq
REF_GENOME=/projects/molonc/huntsman_lab/madouglas/bin/reference_genomes/dlp_refdata/human/GRCh37-lite.fa
TARGETED_PANEL=/projects/molonc/scratch/madouglas/transfers/targeted_panel_sequencing/20220519/DHS-3501Z.roi.renamedchromosomes.bed

SAMPLE_DIRS=$(ls $INPUT_PATH)

for i in $SAMPLE_DIRS
do
	SAMPLE=$(ls ${INPUT_PATH}/$i | egrep -Eo '^[^_]+' | uniq)

	RUNSCRIPTS=${TMP_PATH}/runscripts/${SAMPLE}
	LOGS=${OUTPUT_PATH}/logs/${SAMPLE}
	ALIGNED=/projects/molonc/huntsman_lab/madouglas/cn_signatures_shallowWGS/targetted_panel_seq/aligned/${SAMPLE}
	
	mkdir -p $LOGS $ALIGNED $RUNSCRIPTS

	echo "#!/bin/bash
#SBATCH --partition=upgrade
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=${LOGS}/${SAMPLE}.std_out.toDEL.txt
#SBATCH --job-name=targeted_panel_seq
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --workdir=/projects/molonc/archive/madouglas

bedtools coverage -b $ALIGNED/$SAMPLE.pe.bwa.sorted.dup_rm.bam \
                  -a $TARGETED_PANEL > $ALIGNED/${SAMPLE}_bedtools_coverage_1.txt
bedtools coverage -mean -b $ALIGNED/$SAMPLE.pe.bwa.sorted.dup_rm.bam \
                  -a $TARGETED_PANEL > $ALIGNED/${SAMPLE}_bedtools_coverage_2.txt
bedtools coverage -hist -b $ALIGNED/$SAMPLE.pe.bwa.sorted.dup_rm.bam \
                  -a $TARGETED_PANEL > $ALIGNED/${SAMPLE}_bedtools_coverage_3.txt

" > $RUNSCRIPTS/${SAMPLE}.targetedseq.coverage.run.job
	echo "Submitting $RUNSCRIPTS/${SAMPLE}.targetedseq.coverage.run.job to the cluster"
	sbatch $RUNSCRIPTS/${SAMPLE}.targetedseq.run.job
done