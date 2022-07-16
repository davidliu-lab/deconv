sudo apt update
sudo apt install -y \
    pandoc \
    tree \
    texlive-xetex \
    texlive-fonts-recommended \
    texlive-plain-generic

curl https://raw.githubusercontent.com/GitAlias/gitalias/master/gitalias.txt -o ~/.gitalias 
git config --global include.path ~/.gitalias

# adding stuff to base env
conda install -y -n base -c conda-forge mamba
mamba install -y -n base -c conda-forge \
    black-jupyter=22.6.0 \
    flit \
    isort \
    jupyter-book \
    jupyterlab_code_formatter \
    jupyterlab_miami_nights \
    jupyterlab-lsp \
    nb_conda_kernels \
    nbconvert \
    nbresuse \
    pre-commit \
    python-lsp-server

# to create deconv env:
mamba env create --file conda-env.yml

# to update deconv env:
mamba env update --file conda-env.yml
