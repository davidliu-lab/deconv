#!/bin/bash

docker build --progress=plain --tag grisaitis/bayesprism .

# run this from python:
docker run --rm -v $(pwd):/bayesprism grisaitis/bayesprism \
  --sc_rnaseq_uri 
  --sc_rnaseq_cell_types_uri 
  --sc_rnaseq_cell_states_uri 
  --bulk_rnaseq_uri 

exit()
python -m helpers.bayesprism.example


# list files in current working directory
docker run \
    --rm \
    --volume $(pwd):/bayesprism:ro \
    --workdir /bayesprism \
    grisaitis/bayesprism


    # --entrypoint ls \
    # -la


docker run \
    --rm \
    --volume $(cwd):/bayesprism:ro \
    --workdir /bayesprism \
    grisaitis/bayesprism \
    --reference_scrnaseq_uri gs://liulab/ftp/GSE115978/GSE115978_tpm.csv \
    --reference_scrnaseq_annotations_uri gs://liulab/ \
    --bulkrnaseq_uri gs://liulab/

# docker push grisaitis/bayesprism
