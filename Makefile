
include make.inc

IDIR = ./include

CFLAGS = 	-I$(IDIR)\
			-isystem$(BOOST_DIR)\
			-I$(MAGMA_DIR)/include\
			-L$(MAGMA_DIR)/lib -lmagma\
			-I$(CUDA_DIR)/include\
			-L$(CUDA_DIR)/lib\
			$(OSFLAG)

tbm-run: src/main.cpp
	$(CC) src/main.cpp -o bin/tbm-run $(CFLAGS)

clean:
	rm bin/tbm-run
