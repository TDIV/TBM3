
OSFLAG		= -std=c++11 -m64 -O3 -fPIC -DADD_ -DHAVE_CUBLAS -DMIN_CUDA_ARCH=300 -w\
			-lblas_fix -framework Accelerate -lcublas -lcudart -lstdc++ -lm -fopenmp

CC = g++
BOOST_DIR	= /usr/local/include/boost
CUDA_DIR	= /usr/local/cuda
MAGMA_DIR	= /Users/anrris/opt/magma
