#!/usr/bin/env bash

cd /grphome/grp_tb/processing_scripts/results/untrimmed
multiqc . -o /grphome/grp_tb/processing_scripts/results/documents/untrimmed_qc

cd /grphome/grp_tb/processing_scripts/results/trimmed
multiqc . -o /grphome/grp_tb/processing_scripts/results/documents/trimmed_qc