#!/bin/bash
#PBS -P xf3
#PBS -q normal
#PBS -l walltime=06:00:00
#PBS -l ncpus=2
#PBS -l mem=50G
#PBS -l jobfs=4G
#PBS -l storage=scratch/xf3+gdata/xf3
#PBS -l wd

source /home/106/ht5438/kmerGWAS_perm_snake/gadimod.sh

conda activate snakemake

export LD_LIBRARY_PATH=/g/data/xf3/ht5438/miniconda3/envs/kmerGWAS/lib


set -ueo pipefail
logdir=log
mkdir -p $logdir
export TMPDIR=${PBS_JOBFS:-$TMPDIR}
TARGET=${TARGET:-all}

QSUB="qsub -q {cluster.queue} -l ncpus={cluster.threads} -l jobfs={cluster.jobfs}"
QSUB="$QSUB -l walltime={cluster.time} -l mem={cluster.mem} -N {cluster.name} -l storage=scratch/xf3+gdata/xf3"
QSUB="$QSUB -l wd -j oe -o $logdir -P {cluster.project}" 


snakemake																	\
	-j 1000																	\
	--max-jobs-per-second 2													\
	--cluster-config /home/106/ht5438/kmerGWAS_perm_snake/cluster.yaml			\
	--local-cores ${PBS_NCPUS:-1}											\
	--js /home/106/ht5438/kmerGWAS_perm_snake/jobscript.sh						\
	--nolock																\
	--keep-going															\
	--rerun-incomplete														\
	--use-envmodules														\
	--latency-wait 500														\
	--cluster "$QSUB"														\
	"$TARGET"
	
