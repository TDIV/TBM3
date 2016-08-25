
include make.inc

IDIR = ./include

CFLAGS = 	-I$(IDIR)\
			-isystem$(BOOST_INCLUDE)\
			-I$(MAGMA_DIR)/include\
			-L$(MAGMA_DIR)/lib\
			-I$(CUDA_DIR)/include\
			-L$(CUDA_DIR)/lib\
			$(OSFLAG)

tbm-run: src/main.cpp
	$(CC) src/main.cpp -o bin/tbm-run $(CFLAGS)

clean:
	rm bin/tbm-run
