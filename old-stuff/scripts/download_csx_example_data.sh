#!/bin/bash

pushd /mnt/buckets/liulab/csx_example_files/

export BASE_URL="https://cibersortx.stanford.edu/inc/inc.download.page.handler.php"
# curl -O -J -L {$BASE_URL}?file=NSCLC_PBMCs_Single_Cell_RNA-Seq_Fig2ab.zip
# unzip NSCLC_PBMCs_Single_Cell_RNA-Seq_Fig2ab.zip
# curl -O -J -L {$BASE_URL}?file=RNA-Seq_mixture_melanoma_Tirosh_Fig2b-d.txt

tree -h

popd

