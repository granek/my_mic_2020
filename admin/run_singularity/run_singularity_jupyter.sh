#!/bin/bash

set -u
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
source ${SCRIPT_DIR}/setup.sh

singularity exec $BIND_ARGS $SINGULARITY_IMAGE /usr/local/bin/start.sh jupyter notebook --ip=0.0.0.0 --no-browser
