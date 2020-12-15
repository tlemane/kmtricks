#!/bin/bash

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   kmtricks
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
cmake .. -DKMER_NB_BIT=ALL -DCOUNT_NB_BIT=ALL -DHOWDE=1
make -j4
cd .. 

mkdir -p $PREFIX/bin
mkdir -p $PREFIX/lib
mkdir -p $PREFIX/include/kmtricks

cp -r ./bin/* $PREFIX/bin
cp ./bin/lib/* $PREFIX/lib
cp ./libs/kmtricks/*.hpp $PREFIX/include/kmtricks
cp ./kmtricks.py $PREFIX/bin
chmod +x $PREFIX/bin/kmtricks.py
