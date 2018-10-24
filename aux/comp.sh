#!/bin/bash

/usr/local/cuda-10.0/bin/nvcc -I/usr/local/fmt-5.2.1/include -I/home/pbarletta/chemfiles08/build/install/include -I../include -G -g -O0 -std=c++11 -gencode arch=compute_61,code=sm_61  -M -o utils.d ../src/utils.cpp
/usr/local/cuda-10.0/bin/nvcc -I/usr/local/fmt-5.2.1/include -I/home/pbarletta/chemfiles08/build/install/include -I../include -G -g -O0 -std=c++11 --compile  -x c++ -o  utils.o ../src/utils.cpp

/usr/local/cuda-10.0/bin/nvcc -I/usr/local/fmt-5.2.1/include -I/home/pbarletta/chemfiles08/build/install/include -I../include -G -g -O0 -std=c++11 -gencode arch=compute_61,code=sm_61 -M -o grid.d ../src/grid.cpp
/usr/local/cuda-10.0/bin/nvcc -I/usr/local/fmt-5.2.1/include -I/home/pbarletta/chemfiles08/build/install/include -I../include -G -g -O0 -std=c++11 --compile  -x c++ -o  grid.o ../src/grid.cpp

/usr/local/cuda-10.0/bin/nvcc -I/usr/local/fmt-5.2.1/include -I/home/pbarletta/chemfiles08/build/install/include -I../include -G -g -O0 -std=c++11 -gencode arch=compute_61,code=sm_61 -M -o continuous.d ../src/continuous.cpp
/usr/local/cuda-10.0/bin/nvcc -I/usr/local/fmt-5.2.1/include -I/home/pbarletta/chemfiles08/build/install/include -I../include -G -g -O0 -std=c++11 --compile  -x c++ -o continuous.o ../src/continuous.cpp

# main.cu
/usr/local/cuda-10.0/bin/nvcc -I/usr/local/fmt-5.2.1/include -I../include -I/home/pbarletta/chemfiles08/build/install/include -G -g -O0 -std=c++14 -gencode arch=compute_61,code=sm_61 -M -o main.d ../src/main.cu
/usr/local/cuda-10.0/bin/nvcc -I/usr/local/fmt-5.2.1/include -I../include -I/home/pbarletta/chemfiles08/build/install/include -G -g -O0 -std=c++14 --compile --relocatable-device-code=false -gencode arch=compute_61,code=compute_61 -gencode arch=compute_61,code=sm_61  -x cu -o main.o ../src/main.cu

#################
# link
################
/usr/local/cuda-10.0/bin/nvcc --cudart static -L/usr/local/fmt-5.2.1/build -L/home/pbarletta/chemfiles08/build/install/lib --relocatable-device-code=false -gencode arch=compute_61,code=sm_61 -link -o GANA main.o utils.o grid.o continuous.o -lfmt -lchemfiles -lnetcdf -lCGAL -lmpfr -lgmp -lm

exit 0
