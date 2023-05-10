#!/bin/bash

docker build --progress=plain --tag grisaitis/bayesprism .

python -m helpers.bayesprism.example

exit(0)

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
