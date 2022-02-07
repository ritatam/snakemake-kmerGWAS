#!/bin/bash

source /home/106/ht5438/kmerGWAS_perm_snake/gadimod.sh

export TMPDIR=$PBS_JOBFS

set -ueo pipefail
{exec_job}