#!/bin/bash

# variables
MIXTURE_FILE="Fig2b-WholeBlood_RNAseq.txt"
REFSAMPLE_FILE="Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt"
CSX_DIR="/home/jupyter/csx"

# setup
MIXTURE_FILE_PATH=$(find /mnt/liulab/ -name "$MIXTURE_FILE" | head -n 1)
REFSAMPLE_FILE_PATH=$(find /mnt/liulab/ -name "$REFSAMPLE_FILE" | head -n 1)
rm -r $CSX_DIR/*
mkdir $CSX_DIR/in $CSX_DIR/out
rsync -v $MIXTURE_FILE_PATH $CSX_DIR/in/mixture.txt
rsync -v $REFSAMPLE_FILE_PATH $CSX_DIR/in/refsample.txt

## verify setup
tree -h $CSX_DIR

# run csx
docker run \
    --rm \
    -v $CSX_DIR/in:/src/data \
    -v $CSX_DIR/out:/src/outdir \
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

tree -h $CSX_DIR
