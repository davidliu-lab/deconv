from contextlib import closing
import gzip
import os
import shutil
import urllib.request as request


def download_gz_from_ftp(ftp_url, base_dir):
    # base_dir = './data/downloaded_manually'
    # ftp_url = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127813/suppl/GSE127813_Whole_blood_gene_expression_matrix.txt.gz'

    filename = os.path.basename(ftp_url)
    filepath = os.path.join(base_dir, filename)

    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    with open(filepath, "wb") as f, closing(request.urlopen(ftp_url)) as r:
        shutil.copyfileobj(r, f)

    with gzip.open(filepath, "rb") as f_in:
        filepath_txt = os.path.splitext(filepath)[0]  # remove .gz extension
        with open(filepath_txt, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
