#-------------------------------------------------------------------------------
# compile the UMFPACK demos
#-------------------------------------------------------------------------------

# UMFPACK Copyright (c) Timothy A. Davis.
# All Rights Reserved.  See ../Doc/License.txt for License.

default: libs run

all: libs run hb

all32: libs run hb fortran

all64: libs run hb fortran64

include ../SuiteSparse-5.10.1/SuiteSparse_config/SuiteSparse_config.mk

#-------------------------------------------------------------------------------
# UMFPACK optionally uses the CHOLMOD Partition module
LIB_WITH_CHOLMOD =
ifeq (,$(findstring -DNCHOLMOD, $(UMFPACK_CONFIG)))
    LIB_WITH_CHOLMOD = $(LIB_WITH_PARTITION) $(CUBLAS_LIB) $(CUDART_LIB)
endif

#-------------------------------------------------------------------------------
LIBS = -lm -lrt -Wl,-rpath=../SuiteSparse-5.10.1/lib -L../SuiteSparse-5.10.1/lib -L../SuiteSparse-5.10.1/UMFPACK/Lib -lumfpack -lamd -lsuitesparseconfig $(LIB_WITH_CHOLMOD) $(LAPACK)

libs: metis
	( cd ../SuiteSparse-5.10.1/SuiteSparse_config ; $(MAKE) )
	( cd ../SuiteSparse-5.10.1/AMD ; $(MAKE) library )
	( cd ../SuiteSparse-5.10.1/UMFPACK/Lib ; $(MAKE) )
	- ( cd ../SuiteSparse-5.10.1/CHOLMOD && $(MAKE) library )
	- ( cd ../SuiteSparse-5.10.1/COLAMD && $(MAKE) library )
	- ( cd ../SuiteSparse-5.10.1/CCOLAMD ; $(MAKE) library )
	- ( cd ../SuiteSparse-5.10.1/CAMD ; $(MAKE) library )

metis: ../SuiteSparse-5.10.1/include/metis.h

../SuiteSparse-5.10.1/include/metis.h:
	- ( cd  && $(MAKE) metis )

#-------------------------------------------------------------------------------
# COMPILATION
#-------------------------------------------------------------------------------

Part1: Part1.cpp
	g++ -std=c++11 -I../SuiteSparse-5.10.1/UMFPACK/Include -I../SuiteSparse-5.10.1/include  Part1.cpp $(LIBS)

