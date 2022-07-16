sudo apt update
sudo apt install -y \
    pandoc \
    tree \
    texlive-xetex \
    texlive-fonts-recommended \
    texlive-plain-generic

curl https://raw.githubusercontent.com/GitAlias/gitalias/master/gitalias.txt -o ~/.gitalias 
git config --global include.path ~/.gitalias

conda install -y -n base -c conda-forge mamba
mamba install -y -n base -c conda-forge \
    black \
    black-jupyter \
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

mamba env create --file conda-env.yml
