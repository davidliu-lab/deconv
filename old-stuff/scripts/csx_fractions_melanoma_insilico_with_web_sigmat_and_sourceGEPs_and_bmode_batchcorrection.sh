#!/bin/bash

# variables
MIXTURE_FILE="/mnt/buckets/liulab/derek/simulations/experiments/cibersortx_sim_0.5sd.txt"
SIGMATRIX_FILE=""
SOURCE_GEPS_FILE=""

# setup
CSX_DIR="/mnt/buckets/liulab/csx-runs/$(date '+%Y%m%d_%H%M%S')"
echo "setting up CIBERSORTx in directory:"
echo $CSX_DIR
mkdir -p $CSX_DIR/in $CSX_DIR/out
pushd $CSX_DIR/in
ln -s $MIXTURE_FILE symlink-mymixture.txt
ln -s $SIGMATRIX_FILE symlink-mysigmatrix.txt
ln -s $SOURCE_GEPS_FILE symlink-mysourceGEPs.txt
cp $MIXTURE_FILE mymixture.txt
cp $SIGMATRIX_FILE mysigmatrix.txt
cp $SOURCE_GEPS_FILE mysourceGEPs.txt
popd

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
    --sigmatrix mysigmatrix.txt \
    --mixture mymixture.txt \
    --sourceGEPs mysourceGEPs.txt \
    --rmbatchBmode TRUE \
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
