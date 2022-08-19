# In silico generation of bulk RNA-seq data

Deconvolution methods take bulk RNA-seq expression profiles as input and produce cell type proportion and gene expression as output. In order to evaluate the accuracy of a deconvolution method, one would ideally have paired input and output data from a cohort of samples taken from donors with this disease. To our knowledge, such data does not exist in the literature.



 needs input data (bulk RNA-seq) and known output data (cell type fractional composition and gene expression) to compare with the method's inferences. To our knowledge, such data with paired, known cell type-level data does not exist, at least publicly, for this disease. 

Therefore, we create synthetic, in silico bulk RNA-seq profiles from available single cell RNA-seq of metatastatic melanoma from (#todo citation). The approach we choose 



To evaluate the accuracy of a deconvolution method, ideally one would compare the outputs of a method with ground truth labels, namely cell type composition and gene expression. Unfortunately, we are not aware of bulk RNA-seq data from metastatic melanoma samples where these quantities are known. 

Therefore, we compute synthetic bulk transcriptomes of metastatic melanoma by linearly combining single cell gene expression profiles from a cohort of scRNA-seq of metastatic melanoma.

 one would have known  label data in addition to input data, or, in this context, known cell type composition and gene expression in addition to bulk RNA-seq
