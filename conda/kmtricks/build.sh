#!/bin/bash

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   kmdiff
#   Authors: T. Lemane
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

mkdir build-conda
cd build-conda
cmake .. -DNATIVE=OFF -DCONDA_BUILD=ON -DWITH_MODULES=ON -DWITH_HOWDE=ON -DKMER_LIST="32 64 96 128 160 192 224 256" -DWITH_SOCKS=ON
make -j8
cd ..

mkdir build-conda-debug
cd build-conda-debug
cmake .. -DNATIVE=OFF -DCMAKE_BUILD_TYPE=Debug -DCONDA_BUILD=ON -DWITH_MODULES=ON -DWITH_HOWDE=ON -DKMER_LIST="32 64 96 128 160 192 224 256" -DWITH_SOCKS=ON
make -j8
cd ..

mkdir -p $PREFIX/bin

cp -r ./bin/kmtricks $PREFIX/bin
cp -r ./bin/kmtricks-debug $PREFIX/bin
cp -r ./bin/kmtricks-socks $PREFIX/bin

mkdir build-conda-plugin
cd build-conda-plugin
cmake .. -DNATIVE=OFF -DWITH_PLUGIN=ON -DCONDA_BUILD=ON -DWITH_MODULES=ON -DKMER_LIST="32 64 96 128 160 192 224 256"
make -j8
cd ..

cp -r ./bin/kmtricks $PREFIX/bin/kmtricksp
