## Install from anaconda cloud
```bash
conda create -p ./kmtricks-env
conda activate ./kmtricks-env
conda install -c tlemane kmtricks
```
## Build package from source
```bash
git clone --recursive https://github.com/tlemane/kmtricks
conda create -p ./kmtricks-env
conda activate kmtricks-env
conda install conda-build
conda-build ./kmtricks/conda/kmtricks
```

## Install

```bash
PACKAGE_PATH=$(conda-build kmtricks --output)
conda install ${PACKAGE_PATH}
```