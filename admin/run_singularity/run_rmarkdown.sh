#!/bin/bash

# +
set -u
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
source ${SCRIPT_DIR}/setup.sh

RMARKDOWN="${1}"

echo "Singularity image: $SINGULARITY_IMAGE"
echo "Rmarkdown: $RMARKDOWN"
# -

KNIT_OUT='/workspace/knit_output_html'

singularity exec $BIND_ARGS $SINGULARITY_IMAGE Rscript -e "rmarkdown::render('${RMARKDOWN}', output_dir='${KNIT_OUT}')"
