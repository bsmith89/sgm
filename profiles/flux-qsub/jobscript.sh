#!/bin/sh
# properties = {properties}

source ${{PBS_O_WORKDIR}}/profiles/flux-qsub/activate

{exec_job}
