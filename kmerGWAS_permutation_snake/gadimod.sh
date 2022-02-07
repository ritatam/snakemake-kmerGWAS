module purge
export TMPDIR=${PBS_JOBFS:-/tmp}

function useconda() {
	                eval "$(/g/data/xf3/ht5438/miniconda3/bin/conda shell.zsh hook)"
			        }


useconda
conda activate kmerGWAS
export LD_LIBRARY_PATH=/g/data/xf3/ht5438/miniconda3/envs/kmerGWAS/lib
