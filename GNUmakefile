SUBDIRS := func

include SoftRelTools/arch_spec_root.mk

LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk

LIBLINK := \
-L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) \
-L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR)  \
-l$(PACKAGE)Func \
-lGeometry_service \
-lRecoBase

override CPPFLAGS += -I$(HOME)/Healpix_3.31/src/cxx/optimized_gcc/include \
                     -Wno-unknown-pragmas -fopenmp

override LIBLIBS += $(LOADLIBES) \
 -L$(ART_LIB) \
 -lart_Framework_Services_Optional_TFileService_service \
 -L$(HOME)/lib \
 -lhealpix_cxx -lcxxsupport -lsharp -lfftpack -lc_utils -lcfitsio -lcurl -fopenmp
