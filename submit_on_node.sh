#!/bin/bash
#PBS -P xf3
#PBS -q normal
#PBS -l walltime=12:00:00
#PBS -l ncpus=2
#PBS -l mem=1G
#PBS -l jobfs=4G
#PBS -l storage=scratch/xf3+gdata/xf3
#PBS -l wd

source /home/106/ht5438/snakemake_workflows/kmerGWAS_snake_v2/gadimod.sh
source /g/data/xf3/miniconda/etc/profile.d/conda.sh

conda activate snakemake

export LD_LIBRARY_PATH=/g/data/xf3/miniconda/envs/kmerGWAS/lib


set -ueo pipefail
logdir=log
mkdir -p $logdir
# mkdir -p data/log/
export TMPDIR=${PBS_JOBFS:-$TMPDIR}
TARGET=${TARGET:-all}

QSUB="qsub -q {cluster.queue} -l ncpus={cluster.threads} -l jobfs={cluster.jobfs}"
QSUB="$QSUB -l walltime={cluster.time} -l mem={cluster.mem} -N {cluster.name} -l storage=scratch/xf3+gdata/xf3"
QSUB="$QSUB -l wd -j oe -o $logdir -P {cluster.project}" 


snakemake																	\
	-j 1000																	\
	--max-jobs-per-second 2													\
	--use-conda \
	--cluster-config /home/106/ht5438/snakemake_workflows/kmerGWAS_snake_v2/cluster.yaml			\
	--local-cores ${PBS_NCPUS:-1}											\
	--js /home/106/ht5438/snakemake_workflows/kmerGWAS_snake_v2/jobscript.sh						\
	--nolock																\
	--keep-going															\
	--rerun-incomplete														\
	--cluster "$QSUB"														\
	"$TARGET"
	
