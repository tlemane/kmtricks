#!/bin/bash

SDSL_INC='-DSDSLINC=../../build/sdsl/src/SDSL-build/include'
SDSL_LIB='-DSDSLLIB=../../build/sdsl/src/SDSL-build/lib'
BIN_OUT='-DOUTPUT=../../../howdesbt'

cd ./HowDeSBT
mkdir build
cd build
cmake .. ${SDSL_INC} ${SDSL_LIB} ${BIN_OUT}
make -j4
cd -
