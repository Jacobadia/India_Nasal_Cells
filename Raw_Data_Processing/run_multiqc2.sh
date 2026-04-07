#!/usr/bin/env bash

conda activate multiqc

OUTPUT_DIR="/India_Nasal_cells/Part_1/DOCUMENTS"

cd "/India_Nasal_cells/Part_1/DATA/fastqc_posttrim"

multiqc . \
-o $OUTPUT_DIR

conda deactivate