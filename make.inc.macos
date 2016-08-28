#----------------------------------------------------------|
# Copyright (C) 2016 Yuan-Yen Tai, Hongchul Choi,          |
#                    Jian-Xin Zhu                          |
#                                                          |
# Los Alamos National Laboratory, T-4 Group                |
#                                                          |
#----------------------------------------------------------|

CC = g++

BOOST_DIR	= /usr/local/include/boost
CUDA_DIR	= /usr/local/cuda
MAGMA_DIR	= /Users/anrris/opt/magma

# MacOS library preference 
OSFLAG		= -std=c++11 -m64 -O3 -fPIC -DADD_ -DHAVE_CUBLAS -DMIN_CUDA_ARCH=300 -w\
			-lblas_fix -framework Accelerate -lcublas -lcudart -lstdc++ -lm -fopenmp

