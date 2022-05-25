MIXTURE_FILE="gs://liulab/data/from_amy/tpm_df_combatseq/pivoted_by_william_for_csx.tsv"
REFSAMPLE_FILE="gs://liulab/csx_example_files/Single_Cell_RNA-Seq_Melanoma_SuppFig_3b-d/scRNA-Seq_reference_melanoma_Tirosh_SuppFig_3b-d.txt"

# setup
CSX_DIR=$(mktemp -d)
echo "setting up CIBERSORTx in directory:"
echo $CSX_DIR
mkdir -p $CSX_DIR/in $CSX_DIR/out
gsutil cp $MIXTURE_FILE $CSX_DIR/in/mymixture.txt
gsutil cp $REFSAMPLE_FILE $CSX_DIR/in/myrefsample.txt

## verify setup
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
    --verbose TRUE \
    --absolute TRUE

gsutil rsync -r $CSX_DIR gs://liulab/data/from_amy/tpm_df_combatseq/csx_run_2_absolute
