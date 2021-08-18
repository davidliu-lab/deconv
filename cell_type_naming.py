nice_to_weirds = {
    'Malignant': ['Malignant.cell', 'Mal'],
    'Endothelial': ['Endothelial.cell', 'Endothelial cells', 'Endo.'],
    'CAF': [],
    'T CD8': ['T.CD8', 'T cells CD8'],
    'NK': ['NK cells'],
    'Macrophage': ['Macrophages'],
    'T CD4': ['T.CD4', 'T cells CD4'],
    'B': ['B.cell', 'B cells']
}

weird_to_nice = {
    weird: nice
    for nice in nice_to_weirds
    for weird in nice_to_weirds[nice]
}
