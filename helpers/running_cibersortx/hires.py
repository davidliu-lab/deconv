import tempfile


def run_and_upload(uri_bulk_rnaseq, uri_fractions, uri_refsample):
    with tempfile.TemporaryDirectory() as tmp_dir:
        set_up_inputs_in_dir(tmp_dir, uri_bulk_rnaseq, uri_fractions, uri_refsample)
        # run()
        # upload()
        pass
    pass


def set_up_inputs_in_dir(files_and_uris):
    """
    sets up...
    - bulk rna-seq
    - fractions
    -
    """
    # for file_name, uri in files_and_uris.items():
    pass


def run():
    command_arguments = " ".join(
        [
            "--mixture bulkrnaseq.tsv",
            "--sigmatrix sigmatrix.txt",
            "--cibresults cibresults.txt",
        ]
    )
    command_arguments = """
    --mixture       <file_name>  Mixture matrix [required]
    --sigmatrix     <file_name>  Signature matrix [required]
    --classes       <file_name>  Cell type groupings [optional; default: none]
    --cibresults    <file_name>  Previous CIBERSORTx cell fractions [default: run CIBERSORT]
"""
    pass


if __name__ == "__main__":
    # create_csx_fractions_tsv()
    # run()
    pass
