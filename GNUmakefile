SUBDIRS := func

include SoftRelTools/arch_spec_root.mk

LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

override CPPFLAGS := -I$(NUTOOLS_INC) -I$(NOVADAQ_INC) $(CPPFLAGS)

include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk

LIBLINK := \
-L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) \
-L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR)  \
-l$(PACKAGE)Func \

override LIBLIBS += $(LOADLIBES)
