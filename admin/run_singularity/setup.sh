#!/bin/bash

set -u

# +
DATA_BASE_DIR="/work/mic2021"
WORKSPACE_BASE_DIR="/work/$USER"

SINGULARITY_DIR="$HOME/container_images"
SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-$SINGULARITY_DIR/cachedir}"
SINGULARITY_PULLFOLDER="${SINGULARITY_PULLFOLDER:-$SINGULARITY_DIR/pulldir}"
SINGULARITY_IMAGE="/hpc/group/chsi/container_images/mic_course/jupyter-mic-2021_latest.sif"

CONTAINER_SCRATCH="${HOME}/scratch"

# -------------------------------------------------
# DO NOT MODIFY ANYTHING BELOW HERE!
# -------------------------------------------------

export SINGULARITY_CACHEDIR
export SINGULARITY_PULLFOLDER
# export SINGULARITY_IMAGE="${1:-$SINGULARITY_IMAGE}"

DATA="$DATA_BASE_DIR/rawdata"
WORKSPACE="$WORKSPACE_BASE_DIR/mic2021/workspace"
mkdir -p $WORKSPACE


if [ -d "${DATA}" ]; then
    # BIND_ARGS="--bind ${DATA}:/data:ro"
    BIND_ARGS="--bind ${DATA}:/data"
else
    echo "Make sure DATA exists: $DATA"
    echo "The location of DATA can be controlled by setting DATA_BASE_DIR before running this or in your .bashrc file"
    exit 1
fi

if [ -d "${WORKSPACE}" ]; then
    BIND_ARGS="$BIND_ARGS --bind ${WORKSPACE}:${CONTAINER_SCRATCH}"
else
    echo "Make sure WORKSPACE exists: $WORKSPACE"
    echo "The location of WORKSPACE can be controlled by setting WORKSPACE_BASE_DIR before running this or in your .bashrc file"
    exit 1

fi

#-----------------
SINGTMP="${WORKSPACE_BASE_DIR}/.tmp/`date +%s`_tmp"
mkdir -p $SINGTMP
echo "Binding /tmp to $SINGTMP"

trap "{ rm -rf $SINGTMP; }" EXIT

BIND_ARGS="$BIND_ARGS --bind $SINGTMP:/tmp"
#-----------------

mkdir -p $SINGULARITY_PULLFOLDER $SINGULARITY_CACHEDIR
