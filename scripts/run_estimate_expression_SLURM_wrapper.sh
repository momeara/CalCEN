#!/bin/bash
#SBATCH --job-name=estimate_ca_coexp_expression
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000m 
#SBATCH --time=100:00

### This is a SLURM submission script to estimate gene expression from RNA-seq run data from SRA
###
###   Directories and files used by this script
###     <base_dir>
###        scripts/run_estiamte_expression.R    
###        intermediate_data/
###          todo_runs_<TAG>.tsv
###          estimate_expression_<TAG>/
###     <scratch_dir>
###        <TAG>/estimate_expression_<SLURM_ARRAY_TASK_ID>
###     <
###
###   To submit job
###     cd <base_dir>
###     sbatch \
###       --array 1-<n_runs> \
###       --export=TAG='<YYYYMMDD>',BASE_DIR='$(pwd)',JOB_DIR='<job_dir>' \
###       scripts/run_estimate_expression_SLURM_wrapper.sh . <tag>
###     
###
###
###  A demon of running from the shell
###
###     export BASE_DIR=/home/maom/opt/ca_coexp
###     export JOB_DIR=/scratch/maom_root/maom99/maom/ca_coexp/estimate_expression_20201007
###     export SLURM_ARRAY_TASK_ID=1
###     /home/maom/opt/ca_coexp/scripts/run_estimate_expression_SLURM_wrapper.sh
###
###  this will set:
###
###     TAG: 20201007
###     BASE_DIR: /home/maom/opt/ca_coexp
###     SLURM_ARRAY_TASK_ID: 1
###     Rscript: /home/maom/opt/bin/Rscript
###     RUNS_FNAME: /home/maom/opt/ca_coexp/intermediate_data/todo_runs_20201007.tsv
###     LOGS_DIR: /home/maom/opt/ca_coexp/intermediate_data/estimated_expression_20201007/logs
###     TASK_NAME: estimated_expression_1
###     JOB_DIR: /scratch/maom_root/maom99/maom/ca_coexp/estimated_expression_20201007
###     TASK_DIR: /scratch/maom_root/maom99/maom/ca_coexp/estimated_expression_20201007/estimated_expression_1
###
### and then call
###
###     Rscript /home/maom/opt/ca_coexp/scripts/run_estimate_expression.R \
###         --runs_fname /home/maom/opt/ca_coexp/intermediate_data/todo_runs_20201007.tsv \
###         --run_id 1 \
###         --results_dir /home/maom/opt/ca_coexp/intermediate_data/estimated_expression_20201007 \
###         --logs_dir /home/maom/opt/ca_coexp/intermediate_data/estimated_expression_20201007/logs \
###         --work_dir /scratch/maom_root/maom99/maom/ca_coexp/estimated_expression_20201007/estimated_expression_1


set -e 1

# values passed in as parameters
# TAG
# BASE_DIR
# JOB_DIR



RUNS_FNAME=${BASE_DIR}/intermediate_data/todo_runs_${TAG}.tsv
RESULTS_DIR=${BASE_DIR}/intermediate_data/estimated_expression_${TAG}
mkdir -p ${RESULTS_DIR}

LOGS_DIR=${BASE_DIR}/intermediate_data/estimated_expression_${TAG}/logs
mkdir -p ${LOGS_DIR}

TASK_NAME="estimated_expression_${SLURM_ARRAY_TASK_ID}"
TASK_DIR=${JOB_DIR}/${TASK_NAME}
mkdir -p ${TASK_DIR}

echo "TAG: ${TAG}"
echo "BASE_DIR: ${BASE_DIR}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Rscript: $(which Rscript)"
echo "RUNS_FNAME: ${RUNS_FNAME}"
echo "LOGS_DIR: ${LOGS_DIR}"
echo "JOB_DIR: ${JOB_DIR}"
echo "TASK_NAME: ${TASK_NAME}"
echo "TASK_DIR: ${TASK_DIR}"

cmd="Rscript ${BASE_DIR}/scripts/run_estimate_expression.R \
  --runs_fname ${RUNS_FNAME} \
  --run_id ${SLURM_ARRAY_TASK_ID} \
  --results_dir ${RESULTS_DIR} \
  --logs_dir ${LOGS_DIR} \
  --work_dir ${TASK_DIR} "
echo $cmd

pushd ${BASE_DIR}
$cmd
popd

#rm -rf ${TASK_DIR}
