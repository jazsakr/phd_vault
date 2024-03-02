---
aliases: miniconda, conda
tags: env
type: 
---

# Install miniconda


> [!note] 
> For detailed installation instructions with images, see [[Miniconda Installation Instructions w Images For Beginners]]


To get the link to download [Miniconda](https://docs.anaconda.com/free/miniconda/) installer by right clicking [Miniconda3 Linux 64-bit](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh) under the [Latest Miniconda Installer Links](https://docs.anaconda.com/free/miniconda/#latest-miniconda-installer-links) section and copying the link.

```bash
# download miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# run script to insall conda
bash Miniconda3-latest-Linux-x86_64.sh

# review the license agreement by pressing the ENTER key
# type 'yes' when it asks: Do you accept the license terms? [yes|no]
# check location miniconda will be installed, press ENTER key
# type 'yes' when it asks to automatically initialize conda

source ~/.bashrc
# when you use `conda` for the very first time, `(base)` will appear in the command line prompt.

# update conda to latest version
conda update conda
```

If you get `conda: command not found` but you did download miniconda, run `source ~/anaconda3/bin/activate` to initialize conda. [source](https://askubuntu.com/questions/1143142/conda-init-gives-me-no-action-taken)

## configure miniconda

[To use bioconda channel](https://bioconda.github.io) enter the following lines of commands in the following order (order here matters!):
```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# then update conda again
conda update conda
```

---

# Usage

Refer to the [`conda` cheat sheet](https://docs.conda.io/projects/conda/en/stable/user-guide/cheatsheet.html) for comprehensive list of quick commands or the [Managing Environments page](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#managing-environments) for more information.

## conda create
```bash
# create env
conda create --name myenv

# creat env with specific python version and packages
conda create --name myenv python=3.8 package1=2.1.7 package2

# create env with packages from specific channels
conda create --name myenv bioconda::package2

# create env with lots of specifics
conda create --name myenv python=3.8 package1 bioconda::package2 package3=6.3
```

## importing & exporting env

[Creating an environment from an environment.yml file](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file)
```bash
conda env create -f environment.yml
```

[Example of yml file](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-file-manually):
```python
name: stats2
channels:
  - javascript
dependencies:
  - python=3.9
  - bokeh=2.4.2
  - conda-forge::numpy=1.21.*
  - nodejs=16.13.*
  - flask
  - pip
  - pip:
    - Flask-Testing
```

[Export your active environment to a yml file](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#sharing-an-environment)
```bash
conda env export > environment.yml
```

---

# symlink packages to conda

> [!warning]

When you install packages without using conda, you need to symlink it so that you can call it using conda otherwise you will get this error: `/usr/bin/bash: package_name: command not found`.

To use it globally, symlink it to `base` env AND <font color="red">EVERY</font> environment you want to use it because each environment uses its own bin folder. 
