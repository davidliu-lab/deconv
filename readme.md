### os dependencies

nice things

```shell
sudo apt update
sudo apt install -y \
    pandoc \
    tree \
    texlive-xetex \
    texlive-fonts-recommended \
    texlive-plain-generic
```

mamba

```shell
conda install -y -n base -c conda-forge mamba
```

### conda env

create:

```shell
mamba env create --file conda-env.yml
```

update:

```shell
mamba env update --name deconv --file conda-env.yml
```

activate:

```shell
mamba activate deconv
```

### to install `helpers`

to install the `helpers` package in this repo:

```shell
pip install --verbose --no-build-isolation --editable .
```

### configuring access to google cloud storage with `gcloud auth`

```
gcloud auth application-default login
```

and maybe also

```
gcloud auth login
gcloud config set project keen-dispatch-316219
```
