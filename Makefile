LIB_WITH_CHOLMOD =
ifeq (,$(findstring -DNCHOLMOD, $(UMFPACK_CONFIG)))
    LIB_WITH_CHOLMOD = $(LIB_WITH_PARTITION) $(CUBLAS_LIB) $(CUDART_LIB)
endif
LIBS = -lm -lrt -Wl,-rpath=../SuiteSparse-5.10.1/lib -L../SuiteSparse-5.10.1/lib -L../SuiteSparse-5.10.1/UMFPACK/Lib -lumfpack -lamd -lsuitesparseconfig $(LIB_WITH_CHOLMOD) $(LAPACK)

PLATFORM = X11
HDR = graphics.h easygl_constants.h
ifeq ($(PLATFORM),X11)
   GRAPHICS_LIBS = -lX11
endif

FLAGS = -g -Wall -D$(PLATFORM) -I../SuiteSparse-5.10.1/UMFPACK/Include -I../SuiteSparse-5.10.1/include 

Part1: Part1.cpp $(HDR)
	g++ -std=c++11 Part1.cpp -c $(FLAGS) 
	g++ -c -g -Wall -D$(PLATFORM) graphics.cpp
	g++ -g -Wall -D$(PLATFORM) graphics.o Part1.o $(GRAPHICS_LIBS) -o Part1_exe.o $(LIBS) 

Part2: Part2.cpp $(HDR)
	g++ -std=c++11 Part2.cpp -c $(FLAGS) 
	g++ -c -g -Wall -D$(PLATFORM) graphics.cpp
	g++ -g -Wall -D$(PLATFORM) graphics.o Part2.o $(GRAPHICS_LIBS) -o Part2_exe.o $(LIBS) 

