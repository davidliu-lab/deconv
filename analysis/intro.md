# Overview

## Background

Understanding and predicting the efficacy of immunotherapies in treating metastatic melanoma is an active, unresolved area of investigation. Several investigations have identified genomic and transcriptomic correlates of response at a bulk sample level, including suggestions of biological mechanisms driving response. Furthermore, investigations of single cell transcriptomics have revealed cell type-level features characterizing the development of metastatic melanoma. 

Transcroptomic studies of metasrtatic melanoma and treatment with immunotherapies have identified

Cell-type specific causes and correlates of response to immunotherapy in metastatic melanoma are not well understood. Analyses of bulk transcriptomes of clinical samples have revealed diffentially expressed genes and gene sets associated with clinical outcomes, but the compositional reasons for these differences are not known. That is, it is not known to what extent these sample-level differences are due to (1) differences in cell type proportions in samples and/or (2) differences in gene expression in one or more cell types existing in the samples.

Furthermore, single cell transcriptomes of clinicial samples exist, and have elucidated intra and intertumoral heterogeneity of metastatic melanoma, but such data are not numerous and are mostly not of patients who received immunotherapy. So, while differential expression analyses of bulk transcriptomes of patients with metastatic melanoma have improved our understanding of the causes and correlates of immunotherapies, not much is known with much certainty about the cell type-level phenomena - compositionality and type-specific express - driving these sample-level differences. 

## Introduction

Given the limitations of previous investigations in identifying cell type-level differentiators of the impact of immunotherapies, in this work I investigate the use of deconvolution methods such as CIBERSORTx in inferring cell type composition and type-specific gene expression in bulk RNA-seq of clinical metastatic melanoma samples. With this, we might hope to improve our understanding of the cell type-level processes responsible for, or at least predictive of, the success or failure of immunotherapy in treating this disease. 

# Methods

## Generation of in silico bulk RNA-seq profiles of metastatic melanoma

Deconvolution methods take bulk RNA-seq expression profiles as input and produce cell type proportion and gene expression as output. In order to evaluate the accuracy of a deconvolution method, one would ideally have paired input and output data from a cohort of samples taken from donors with this disease. To our knowledge, such data does not exist in the literature.



 needs input data (bulk RNA-seq) and known output data (cell type fractional composition and gene expression) to compare with the method's inferences. To our knowledge, such data with paired, known cell type-level data does not exist, at least publicly, for this disease. 

Therefore, we create synthetic, in silico bulk RNA-seq profiles from available single cell RNA-seq of metatastatic melanoma from (#todo citation). The approach we choose 




To evaluate the accuracy of a deconvolution method, ideally one would compare the outputs of a method with ground truth labels, namely cell type composition and gene expression. Unfortunately, we are not aware of bulk RNA-seq data from metastatic melanoma samples where these quantities are known. 

Therefore, we compute synthetic bulk transcriptomes of metastatic melanoma by linearly combining single cell gene expression profiles from a cohort of scRNA-seq of metastatic melanoma.

 one would have known  label data in addition to input data, or, in this context, known cell type composition and gene expression in addition to bulk RNA-seq
