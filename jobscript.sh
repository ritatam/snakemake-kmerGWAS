#!/bin/bash

source /home/106/ht5438/snakemake_workflows/kmerGWAS_snake_v2/gadimod.sh

export TMPDIR=$PBS_JOBFS

set -ueo pipefail
{exec_job}
