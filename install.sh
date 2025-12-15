#!/bin/bash

check_conda(){
    which conda &> /dev/null
    if [ $? -ne 0 ]
    then
        echo "$1$reset is required with -e <env>"
        exit 1
    fi
}

function kmtricks_build ()
{
  echo "Type=${1}, Modules=${2}, Socks=${3}, HowDe=${4}, Tests=${5}, Static=${6}, K=${7}, C=${8}, Native=${11}, Run=${10}, Plugin=${12}, j=${9}"

  mkdir build
  cd build

  cmake .. -DCMAKE_BUILD_TYPE=${1} \
           -DWITH_MODULES=${2} \
           -DCOMPILE_TESTS=${5} \
           -DSTATIC=${6} \
           -DKMER_LIST="${7}" \
           -DMAX_C=${8} \
           -DNATIVE=${11} \
           -DWITH_PLUGIN=${12}

  if [[ ${13} == 1 ]]; then
    make plugins
  else
    make -j${9}
  fi

  if [[ ${10} == 1 ]]; then
    ctest --verbose
  fi
}

function kmtricks_conda_build ()
{
  conda create -p ./km_conda_build
  conda activate ./km_conda_build
  if [ "$(uname)" == "Darwin" ]; then
    conda install -y clangxx_osx-64=11.1.0 cmake zlib -c conda-forge
    export CC=./km_conda_build/bin/clang
    export CXX=./km_conda_build/bin/clang++

  else
    conda install -y gxx_linux-64=9.3.0 cmake zlib -c conda-forge
    export CC=./km_conda_build/bin/x86_64-conda_cos6-linux-gnu-gcc
    export CXX=./km_conda_build/bin/x86_64-conda_cos6-linux-gnu-g++
  fi

  mkdir build_conda
  cd build_conda

  cmake .. -DCMAKE_BUILD_TYPE=${1} \
           -DWITH_MODULES=${2} \
           -DCOMPILE_TESTS=${5} \
           -DSTATIC=${6} \
           -DKMER_LIST="${7}" \
           -DMAX_C=${8} \
           -DNATIVE=${11} \
           -DCONDA_BUILD=ON \
           -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} \
           -DCMAKE_LIBRARY_PATH=${CONDA_PREFIX} \
           -DCMAKE_INCLUDE_PATH=${CODNA_PREFIX}

  if [[ ${13} == 1 ]]; then
    make plugins
  else
    make -j${9}
  fi

  if [[ ${10} == 1 ]]; then
    ctest --verbose
  fi
}

function usage ()
{
  echo "kmtricks build script - v1.5.1."
  echo "Usage: "
  echo "  ./install.sh [-r str] [-k LIST[int]] [-t int] [-c int] [-j int] [-m] [-s] [-n]
  [-e] [-p] [-q] [-h]"
  echo "Options: "
  echo "  -r <Release|Debug> -> build type {Release}."
  echo "  -k <LIST[INT]>     -> k-mer size {\"32 64 96 128\"}."
  echo "  -t <0|1|2>         -> tests: 0 = disabled, 1 = compile, 2 = compile and run {2}."
  echo "  -c <1|2|4>         -> byte per count {4}."
  echo "  -j <INT>           -> nb threads {8}."
  echo "  -m                 -> build all modules {disabled}."
  echo "  -s                 -> static build {disabled}."
  echo "  -n                 -> disable -march=native {enabled}."
  echo "  -e                 -> use conda to install compilers/dependencies {disabled}."
  echo "  -p                 -> with plugin support {disabled}."
  echo "  -q                 -> build plugins only {disabled}."
  echo "  -h                 -> show help."
  exit 1
}

rel="Release"
deb="Debug"

mode="Release"
ks="32 64 96 128"
count=4
dev="OFF"
static="OFF"
tests=2
tests_str="ON"
tests_run=1
jopt=8
socks="OFF"
howde="OFF"
modules="OFF"
native="ON"
conda=0
plugin="OFF"
plugin_only=0

while getopts "r:k:t:c:j:nmspqeh" option; do
  case "$option" in
    r)
      mode=${OPTARG}
      [[ ${mode} != ${rel} && ${mode} != ${deb} ]] && usage
      ;;
    k)
      ks=${OPTARG}
      ;;
    t)
      tests=${OPTARG}
      [[ ${tests} == 0 ]] || [[ ${tests} == 1 ]] || [[ ${tests} == 2 ]] || usage
      [[ ${tests} == 0 ]] && tests_run=0
      [[ ${tests} == 0 ]] && tests_str="OFF"
      [[ ${tests} == 1 ]] && tests_run=0
      ;;
    c)
      count=${OPTARG}
      [[ ${count} == 1 ]] || [[ ${count} == 2 ]] || [[ ${count} == 4 ]] || usage
      ;;
    j)
      jopt=${OPTARG}
      ;;
    m)
      modules="ON"
      ;;
    s)
      static="ON"
      ;;
    n)
      native="OFF"
      ;;
    e)
      conda=1
      ;;
    p)
      plugin="ON"
      ;;
    q)
      plugin_only=1
      tests_str="OFF"
      tests_run=0
      ;;
    h)
      usage
      ;;
    *)
      usage
      ;;
  esac
done

count=$((2**(${count}*8)-1))

if [[ ${conda} -eq 1 ]]; then
  check_conda
  conda_install_path=$(conda info | grep -i 'base environment')
  conda_install_path=$(echo ${conda_install_path} | cut -d' ' -f4)
  source ${conda_install_path}/etc/profile.d/conda.sh
  kmtricks_conda_build ${mode} ${modules} ${socks} ${howde} ${tests_str} ${static} "${ks}" ${count} ${jopt} ${tests_run} ${native} ${plugin} ${plugin_only}
else
  kmtricks_build ${mode} ${modules} ${socks} ${howde} ${tests_str} ${static} "${ks}" ${count} ${jopt} ${tests_run} ${native} ${plugin} ${plugin_only}
fi
