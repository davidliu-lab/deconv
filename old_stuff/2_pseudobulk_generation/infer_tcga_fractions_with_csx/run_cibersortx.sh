MIXTURE_FILE="gs://liulab/tmp/mixtures_tcga_skcm_hg19_tpm.txt"
REFSAMPLE_FILE="gs://liulab/tmp/refsample_jerby_arnon.txt"

# setup
CSX_DIR=$(mktemp -d)
echo "setting up CIBERSORTx in directory:"
echo $CSX_DIR
mkdir -p $CSX_DIR/in $CSX_DIR/out
gsutil cp $MIXTURE_FILE $CSX_DIR/in/mymixture.txt
gsutil cp $REFSAMPLE_FILE $CSX_DIR/in/myrefsample.txt

tree -h $CSX_DIR

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
    --rmbatchBmode TRUE \
    --verbose TRUE

tree -h $CSX_DIR

gsutil rsync -r $CSX_DIR gs://liulab/tmp/csx_results
