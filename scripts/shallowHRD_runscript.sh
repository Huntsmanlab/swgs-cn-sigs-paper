#!/bin/bash

# shallowHRD_runscript.sh
# Description: Batch script for running shallowHRD on sWGS data

# All paths are absolute for now
# Important main paths to set
PROJECT_DIR='cn_signatures_shallowWGS'
TMP_PATH=/projects/molonc/scratch/madouglas/cn_signatures
PROJECT_PATH=/projects/molonc/huntsman_lab/madouglas/${PROJECT_DIR}


IFS=$'\t'
while read sample_id batch repeat_run tissue sample_type status hist cancer_type grade sequence library fastq_path_1 fastq_path_2
do
	# if [[ $batch == 1 ]] || [[ $batch == 2 ]] || [[ $batch == 3 ]] || \
	# 	[[ $batch == 4 ]] || [[ $batch == 5 ]] || [[ $batch == 6 ]] || \
	# 	[[ $batch == 7 ]] || [[ $batch == 8 ]] || [[ $batch == 9 ]] || \
	# 	[[ $batch == 10 ]] || [[ $batch == 11 ]];
	if [[ $repeat_run == 1 ]];
	then
  		# Make pipeline directories
  		RUNSCRIPTS=${TMP_PATH}/runscripts/${sample_id}
  		LOGS=${TMP_PATH}/logs/${sample_id}
			INTERIM_OUTPUT=${TMP_PATH}/shallowHRD/${sample_id}
  		mkdir -p $RUNSCRIPTS $LOGS $INTERIM_OUTPUT

			echo "#!/bin/bash
echo 'Running shallowHRD on ${sample_id}'
Rscript_4.1 ${PROJECT_PATH}/scripts/shallowHRD_hg19_1.13_QDNAseq_chrX.R ${TMP_PATH}/shallowHRDinputfiles/${sample_id}.bam_ratio.txt ${INTERIM_OUTPUT} ${PROJECT_PATH}/data/cytoband_adapted_hg19.csv
" > $RUNSCRIPTS/analysis2.bash

			echo "#!/bin/bash
#SBATCH --partition=upgrade
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=${LOGS}/${sample_id}.std_out2.txt
#SBATCH --job-name=shallowHRD_run
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32
#SBATCH --workdir=${TMP_PATH}
sh $RUNSCRIPTS/analysis2.bash" > $RUNSCRIPTS/run2.job
			echo "Submitting $RUNSCRIPTS/run2.job to the cluster"
			sbatch $RUNSCRIPTS/run2.job
  fi
done < "/projects/molonc/huntsman_lab/madouglas/cn_signatures_shallowWGS/metadata/metadata_unix.tsv"
