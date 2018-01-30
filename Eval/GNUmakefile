include SoftRelTools/arch_spec_root.mk

LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

override CPPFLAGS := -I$(NUTOOLS_INC) $(CPPFLAGS)

LIBLINK    := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) 


########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk
include SoftRelTools/arch_spec_novadaq.mk

override LIBLIBS+= $(LOADLIBES)
