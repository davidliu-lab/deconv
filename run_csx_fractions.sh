#!/bin/bash

export CSX_INPUT_DIR="/home/jupyter/csx/input"
export CSX_OUTPUT_DIR="/home/jupyter/csx/output"
mkdir $CSX_INPUT_DIR
mkdir $CSX_OUTPUT_DIR

export MIXTURE_FILE="Fig2b-WholeBlood_RNAseq.txt"
export REFSAMPLE_FILE="Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt"

rsync -v $(find /mnt/liulab/ -name "$MIXTURE_FILE") $CSX_INPUT_DIR/mixture.txt
rsync -v $(find /mnt/liulab/ -name "$REFSAMPLE_FILE") $CSX_INPUT_DIR/refsample.txt

ls -hl $CSX_INPUT_DIR

docker run \
    --rm \
    -v $CSX_INPUT_DIR:/src/data \
    -v $CSX_OUTPUT_DIR:/src/outdir \
    --user "$(id -u):$(id -g)" \
    cibersortx/fractions:latest \
    --username lyronctk@stanford.edu \
    --token dfeba2c8b9d61daebee5fa87026b8e56 \
    --single_cell TRUE \
    --refsample refsample.txt \
    --mixture mixture.txt \
    --rmbatchSmode TRUE \
    --verbose TRUE

#     --perm 10 \
#     --fraction 0 \
#     --sourceGEPs signature_matrix.txt
