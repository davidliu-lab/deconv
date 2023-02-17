# installation and machine setup

## install visual studio code

see https://code.visualstudio.com/docs/setup/linux#_installation

```shell
curl -L https://go.microsoft.com/fwlink/?LinkID=760868 --output visual_studio_code.deb
sudo apt install ./visual_studio_code.deb
```

and develop from a remote host (https://code.visualstudio.com/docs/remote/vscode-server)

```shell
code tunnel --accept-server-license-terms
```

## os dependencies

### a conda installation

this project assumes you have a conda installation. for a good one, check out [conda-forge/miniforge](https://github.com/conda-forge/miniforge).

### add git aliases

```shell
curl https://raw.githubusercontent.com/GitAlias/gitalias/master/gitalias.txt -o ~/.gitalias
git config --global include.path ~/.gitalias
```

### add some jupyter-related dependencies to the base conda env

```shell
conda install -y -n base -c conda-forge mamba

mamba env update --file conda-env-base-extras.yml

mamba update -n base -c conda-forge --update-all
```

### others

```shell
sudo apt update

sudo apt install -y \
    pandoc \
    tree \
    texlive-xetex \
    texlive-fonts-recommended \
    texlive-plain-generic
```

## conda env for project code

```shell
mamba env create --file conda-env.yml

mamba env update --file conda-env.yml

mamba activate deconv
```

### add `deconv` kernel to jupyter

```shell
# in the base env, in which jupyter lab runs
python -m ipykernel install --user --name=deconv
# or, to mimic nb_conda_kernels
python -m ipykernel install --user --name=conda-env-deconv-py
```

### to install `helpers`

```shell
# in the deconv env
pip install --verbose --no-build-isolation --editable .
```

### configuring access to google cloud storage with `gcloud auth`

```shell
gcloud auth application-default login
# and maybe also...
gcloud auth login
gcloud config set project keen-dispatch-316219
```

## restart jupyter lab server

```shell
sudo service jupyter status
sudo service jupyter restart
```

Source: https://cloud.google.com/vertex-ai/docs/general/troubleshooting-workbench#restart_the_jupyter_service
