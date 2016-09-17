#----------------------------------------------------------|
# Copyright (C) 2016 Yuan-Yen Tai, Hongchul Choi,          |
#                    Jian-Xin Zhu                          |
#                                                          |
# Los Alamos National Laboratory, T-4 Group                |
#                                                          |
#----------------------------------------------------------|

include make.inc

IDIR = ./include

CFLAGS = 	-I$(IDIR)\
			-isystem$(BOOST_DIR)\
			-I$(MAGMA_DIR)/include\
			-L$(MAGMA_DIR)/lib -lmagma\
			-I$(CUDA_DIR)/include\
			-L$(CUDA_DIR)/lib\
			-L$(CUDA_DIR)/lib64\
			$(OSFLAG)

tbm-run: src/main-tbm-run.cpp
	$(CC) src/main-tbm-run.cpp -o bin/tbm-run $(CFLAGS)

tbm-wannier: src/main-tbm-wannier.cpp
	$(CC) src/main-tbm-wannier.cpp -o bin/tbm-wannier $(CFLAGS)

clean:
	rm bin/tbm-run
