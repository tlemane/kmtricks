## Build package

```bash
conda create -p ./kmtricks-env
conda activate kmtricks-env
conda install conda-build
conda build ./kmtricks
```

## Install

```bash
PACKAGE_PATH=$(conda-build kmtricks --output)
conda install ${PACKAGE_PATH}
```