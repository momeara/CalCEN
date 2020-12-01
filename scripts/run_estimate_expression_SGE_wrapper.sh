#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N estimate_expression
#$ -e log.error
#$ -o log.output

### This is a SGE submission script to estimate gene expression from RNA-seq run data from SRA
###
###   Directories and files used by this script
###     <base_dir>
###        scripts/run_estiamte_expression.R    
###        intermediate_data/
###          estimate_expression_<TAG>/
###     /scratch/$(whoami)/
###        <TAG>/estimate_expression_<SGE_TASK_ID>
###
###   To submit job
###     cd <base_dir>
###     qsub -t 1-<n_runs> scripts/run_estimate_expression-wrapper.sh . <tag>
###     



set -e 1

PERSIST=/mnt/nfs/work/momeara/collaborations/CalCEN
TAG=180611

RUNS_FNAME=${PERSIST}/intermediate_data/todo_runs_180611.tsv
RESULTS_DIR=${PERSIST}/intermediate_data/estimated_expression_${TAG}
mkdir -p ${RESULTS_DIR}

LOGS_DIR=${PERSIST}/intermediate_data/estimated_expression_${TAG}/logs
mkdir -p ${LOGS_DIR}

TASK_NAME="estimated_expression_${SGE_TASK_ID}"
SCRATCH_DIR=/scratch
if [ ! -d $SCRATCH_DIR ]; then
    SCRATCH_DIR=/tmp
fi
TASK_DIR=$SCRATCH_DIR/$( whoami )/estimated_expression_${TAG}/${TASK_NAME}
mkdir -p ${TASK_DIR}


echo "PERSIST: ${PERSIST}"
echo "TAG: ${TAG}"
echo "SGE_TASK_ID: ${SGE_TASK_ID}"

cmd="/mnt/nfs/home/momeara/opt/bin/Rscript ${PERSIST}/scripts/run_estimate_expression.R \
  --runs_fname ${RUNS_FNAME} \
  --run_id ${SGE_TASK_ID} \
  --results_dir ${RESULTS_DIR} \
  --logs_dir ${LOGS_DIR} \
  --work_dir ${TASK_DIR} "
$cmd

rm -rf ${TASK_DIR}
