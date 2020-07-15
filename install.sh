#!/bin/bash

if [ -z "$1" ]; then
    mkdir build ; cd build ; cmake .. ; make -j8 ; cd -
    echo "Done. Executables at ./bin"
else
    mkdir build ; cd build ; cmake .. -DCMAKE_BUILD_TYPE=Debug ; make -j8 ; cd -
    echo "Done. Debug Mode"
fi
