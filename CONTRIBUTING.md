# Contributing

Below are some checks to run before sending in Pull Requests (PR) and that we run 
when making new releases.

## Code formatting

```sh
pip install pre-commit
pre-commit install
```

Now hooks enforce, when you `git commit`, that:
    * the python code is formatted with [Black](https://github.com/psf/black>) 
    * the C/C++ code is formatted with [clang-format](https://clang.llvm.org/docs/ClangFormat.html)

## Code testing

### Front-end
```sh
pytest
```

from the root source dir.

### Back-end

```sh
mkdir cmake-build && cd cmake-build
cmake -DCMAKE_BUILD_TYPE=REL_WITH_ASSERTS ..
make -j 4 test_main
ctest -V
```
