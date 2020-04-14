SUBDIRS := func

include SoftRelTools/arch_spec_root.mk

LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_nutools.mk
include SoftRelTools/arch_spec_art.mk

LIBLINK := \
-L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) \
-L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR)  \
-l$(PACKAGE)Func \
-lGeometry_service \
-lRecoBase

# -pipe is *only* there to work-around people filling up /tmp,
# although I know no reason why it would be harmful.  I think it's
# just a memory/disk tradeoff.
override CPPFLAGS += -I$(HOME)/Healpix_3.31/src/cxx/basic_gcc/include \
                     -Wno-unknown-pragmas -fopenmp \
                     -pipe -O3

override LIBLIBS += $(LOADLIBES) \
 -L$(ART_LIB) \
 -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) \
 -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) \
 -lart_Framework_Services_Optional_TFileService_service \
 -lCelestialLocator_service \
 -L$(HOME)/lib \
 -lhealpix_cxx -lcxxsupport -lsharp -lfftpack -lc_utils -lcfitsio -lcurl -fopenmp
