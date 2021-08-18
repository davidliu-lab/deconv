#!/bin/bash

# variables
MIXTURE_FILE="/mnt/buckets/liulab/csx_example_files/Fig2ab-NSCLC_PBMCs/Fig2b-WholeBlood_RNAseq.txt"
REFSAMPLE_FILE="/mnt/buckets/liulab/csx_example_files/Fig2ab-NSCLC_PBMCs/Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt"

# setup
CSX_DIR="/mnt/buckets/liulab/csx-runs/$(date '+%Y%m%d_%H%M%S')"
echo "setting up CIBERSORTx in directory:"
echo $CSX_DIR
mkdir -p $CSX_DIR/in $CSX_DIR/out
cp $MIXTURE_FILE $CSX_DIR/in/mymixture.txt
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
    --mixture mymixture.txt \
    --rmbatchSmode TRUE \
    --verbose TRUE

#     --perm 10 \
#     --fraction 0 \
#     --sourceGEPs signature_matrix.txt

tree -h $CSX_DIR


# docker run \
#     --rm \
#     -v $CSX_DIR.2/in:/src/data \
#     -v $CSX_DIR.2/out:/src/outdir \
#     --user "$(id -u):$(id -g)" \
#     cibersortx/fractions:latest \
#     --username lyronctk@stanford.edu \
#     --token dfeba2c8b9d61daebee5fa87026b8e56 \
#     ----sigmatrix signature_matrix.txt \
#     --mixture mixture.txt \
#     --verbose TRUE
