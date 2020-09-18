# Kinetic Monte Carlo Code

1. The compile instruction is in the DOC directory.
2. The input example is test.config


## Build requirements

As of 18/09/20, the code is built and tested using :

- `gcc` (GCC) v9.3.1
- `cmake` version v3.17.4
- `boost` v1.69.0 and `boost-devel` v1.69.0 dnf packages 
- Fedora 31 OS

## Quick start

```
sudo dnf install gcc cmake boost boost-devel
git clone git@github.com:arn-all/Kinetic_Monte_Carlo_Code.git KMC
cd KMC
mkdir bin
cd bin
CC=gcc
CXX=g++
cmake ../src
make -j 4
```
