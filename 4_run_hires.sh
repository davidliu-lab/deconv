#!/bin/bash

# variables
MIXTURE="gs://liulab/csx_experiments/varying_parameters/with_bmode/CIBERSORTx_Mixtures_Adjusted.txt"
SIGMATRIX="gs://liulab/csx_experiments/varying_parameters/with_bmode/CIBERSORTx_screfsampletirosh_inferred_phenoclasses.CIBERSORTx_screfsampletirosh_inferred_refsample.bm.K999.txt"
CIBRESULTS="gs://liulab/csx_experiments/varying_parameters/with_bmode/CIBERSORTx_Adjusted.txt"

# setup
CSX_DIR=$(mktemp -d)
echo "setting up CIBERSORTx in directory:"
echo $CSX_DIR
mkdir -p $CSX_DIR/in $CSX_DIR/out
gsutil cp $MIXTURE $CSX_DIR/in/mymixture.txt
gsutil cp $SIGMATRIX $CSX_DIR/in/mysigmatrix.txt
# gsutil cp $CIBRESULTS $CSX_DIR/out/mycibresults.txt

## verify setup
tree -h $CSX_DIR

# run csx
docker run \
    --rm \
    -it \
    -v $CSX_DIR/in:/src/data \
    -v $CSX_DIR/out:/src/outdir \
    cibersortx/hires:latest \
    --username lyronctk@stanford.edu \
    --token dfeba2c8b9d61daebee5fa87026b8e56 \
    --mixture mymixture.txt \
    --sigmatrix mysigmatrix.txt \
    --variableonly TRUE

#     --cibresults mycibresults.txt \
#     --user "$(id -u):$(id -g)" \
#     --rmbatchBmode TRUE \
#     --rmbatchSmode TRUE \
#     --sourceGEPs signature_matrix.txt

tree -h $CSX_DIR

sudo chown -R $(id -u):$(id -g) $CSX_DIR

# rsync -adv $CSX_DIR ./4_results/
gsutil rsync -d -r $CSX_DIR gs://liulab/csx_experiments/imputation
