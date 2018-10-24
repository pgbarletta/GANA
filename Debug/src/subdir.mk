################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../src/main.cu 

CPP_SRCS += \
../src/continuous.cpp \
../src/grid.cpp \
../src/utils.cpp 

OBJS += \
./src/continuous.o \
./src/grid.o \
./src/main.o \
./src/utils.o 

CU_DEPS += \
./src/main.d 

CPP_DEPS += \
./src/continuous.d \
./src/grid.d \
./src/utils.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-10.0/bin/nvcc -I/home/pbarletta/labo/18/cuGANA/include -I/home/pbarletta/labo/18/cuGANA/chemfiles/build/install/include -I/home/pbarletta/labo/18/cuGANA/fmt/include -G -g -O0 -gencode arch=compute_61,code=sm_61  -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-10.0/bin/nvcc -I/home/pbarletta/labo/18/cuGANA/include -I/home/pbarletta/labo/18/cuGANA/chemfiles/build/install/include -I/home/pbarletta/labo/18/cuGANA/fmt/include -G -g -O0 --compile  -x c++ -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/%.o: ../src/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-10.0/bin/nvcc -I/home/pbarletta/labo/18/cuGANA/include -I/home/pbarletta/labo/18/cuGANA/chemfiles/build/install/include -I/home/pbarletta/labo/18/cuGANA/fmt/include -G -g -O0 -gencode arch=compute_61,code=sm_61  -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-10.0/bin/nvcc -I/home/pbarletta/labo/18/cuGANA/include -I/home/pbarletta/labo/18/cuGANA/chemfiles/build/install/include -I/home/pbarletta/labo/18/cuGANA/fmt/include -G -g -O0 --compile --relocatable-device-code=false -gencode arch=compute_61,code=compute_61 -gencode arch=compute_61,code=sm_61  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


