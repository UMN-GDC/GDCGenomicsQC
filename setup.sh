#!/bin/bash

# setup the conda env
# conda env export > full_environment.yml

conda env create -f environment.yml

# install popvae separately
pip install git+https://github.com/coffm049/popvae/tree/master
