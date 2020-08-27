# kmtricks header-only

These headers allow to interact with some kmtricks files. Example are provided [here](../snippets) and in the [tests](../../tests/libs). Full documentation will be available soon.

## Available

* sequences.hpp : super-k-mers, k-mers and minimizers representation
* skreader.hpp : read super-k-mers partitions 
* repartition.hpp : request minimizer repartition
* merger.hpp : k-way merger
* bitmatrix.hpp : bit matrix representation with tranposition support
* kmlib.hpp : shortcut to include everything

## Usage

### As header-only

You only need to copy headers into your project and include "kmlib.hpp".

### As static library

If you plan to use it in multiple compilation units, please use the static library (after build at `bin/lib/libkmtricks.a`) with `_KM_LIB_INCLUDE_` to prevent multiple definitions and link it with `-lkmtricks`.

```cpp
#define _KM_LIB_INCLUDE_
#include "kmtricks/kmlib.hpp"
```

To build on the kmtricks library use :  
```bash
mkdir build ; cmake .. ; cd build ; cmake --build . --target kmtricks
```

