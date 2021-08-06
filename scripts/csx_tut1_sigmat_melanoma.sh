#!/bin/bash

# variables
REFSAMPLE_FILE="/mnt/buckets/liulab/csx_example_files/Single_Cell_RNA-Seq_Melanoma_SuppFig_3b-d/scRNA-Seq_reference_melanoma_Tirosh_SuppFig_3b-d.txt"

# setup
CSX_DIR="/mnt/buckets/liulab/csx-runs/$(date '+%Y%m%d_%H%M%S')"
echo "setting up CIBERSORTx in directory:"
echo $CSX_DIR
mkdir -p $CSX_DIR/in $CSX_DIR/out
cp $REFSAMPLE_FILE $CSX_DIR/in/myrefsample.txt

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
    --refsample myrefsample.txt \
    --replicates 5 \
    --sampling 0.5 \
    --fraction 0.75 \
    --k.max 999 \
    --q.value 0.01 \
    --G.min 300 \
    --G.max 500 \
    --filter FALSE \
    --verbose TRUE \
    --QN FALSE

#     --perm 10 \
#     --fraction 0 \
#     --sourceGEPs signature_matrix.txt

tree -h $CSX_DIR
