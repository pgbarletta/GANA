################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../include/GANA/kernels.cu 

OBJS += \
./include/GANA/kernels.o 

CU_DEPS += \
./include/GANA/kernels.d 


# Each subdirectory must supply rules for building sources it contributes
include/GANA/%.o: ../include/GANA/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-10.0/bin/nvcc -I/home/pbarletta/labo/18/GANA/include -I/home/pbarletta/chemfiles08/build/install/include -I/usr/local/fmt-5.2.1/include -G -g -O0 -gencode arch=compute_61,code=sm_61 -m64 -odir "include/GANA" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-10.0/bin/nvcc -I/home/pbarletta/labo/18/GANA/include -I/home/pbarletta/chemfiles08/build/install/include -I/usr/local/fmt-5.2.1/include -G -g -O0 --compile --relocatable-device-code=false -gencode arch=compute_61,code=compute_61 -gencode arch=compute_61,code=sm_61 -m64  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


