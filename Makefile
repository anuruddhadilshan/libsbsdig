#------------------------------------------------------------------------------
include Makefile.arch

#------------------------------------------------------------------------------

NAME    := libsbsdig

SOLINCLUDE := -I$(shell pwd)/src -I$(shell pwd)/libsbsgem

#------------------------------------------------------------------------------
# Hall A Analyzer

# Analyzer default location,
ANALYZER ?= $(HOME)/ANALYZER
# Possible Analyzer header locations, will be used in the order found
ANAINCDIRS  := $(wildcard  $(addprefix $(ANALYZER)/, include src hana_decode hana_scaler evio))
ifeq ($(strip $(ANAINCDIRS)),)
  $(error No Analyzer header files found. Check $$ANALYZER)
endif
#SBSINCDIRS = $(SBS_ANALYSIS)
#ifeq ($(strip $(SBSINCDIRS)),)
#  $(error No Analyzer header files found. Check $$SBS_ANALYSIS)
#endif


#------------------------------------------------------------------------------
# EVIO => comment evio stuff
# 
# PLATFORM = $(shell uname -s)-$(shell uname -i)
# EVIO default locations as a last resort if $EVIO isn't set in build env
# EVIO ?= $(CODA)
# EVIO ?= ../libevio
# ifdef EVIO_INCDIR
# 	EVIOINC= $(EVIO_INCDIR)
# else
# 	# Possible EVIO header directories, will be used in the order found
# 	EVIOINC := $(wildcard $(addprefix $(EVIO)/, include src/libsrc src/libsrc++))
# endif
# ifdef EVIO_LIBDIR
# 	EVIOLIB= $(EVIO_LIBDIR)
# else
# 	# Possible EVIO library locations, the first one found will be used
# 	EVIOLIB := $(firstword $(wildcard $(addprefix $(EVIO)/, $(PLATFORM)/lib lib)))
# endif
# ifeq ($(strip $(EVIOINC)),)
#   $(error No EVIO header files found. Check $$EVIO)
# endif
# ifeq ($(strip $(EVIOLIB)),)
#   $(error No EVIO library directory found. Check $$EVIO)
# endif
# ifeq (debug,$(findstring debug,$(ROOTBUILD)))
#   DBGSFX = -dbg
#   CXXFLAGS += -O0
# else
#   CXXFLAGS += -O2 -g #-DNDEBUG
# endif
# SOLINCLUDE += $(addprefix -I, $(EVIOINC) )
# LDFLAGS  += -L$(EVIOLIB) -levioxx$(DBGSFX) -levio$(DBGSFX) -lz -lexpat

# Some of the analyzer include dirs conflict with headers in
# EVIO
SOLINCLUDE += $(addprefix -I, $(ANAINCDIRS) )
#SOLINCLUDE += $(addprefix -I, $(SBSINCDIRS) )

#------------------------------------------------------------------------------

CXXFLAGS += $(SOLINCLUDE)

DICT	= $(NAME)_dict

GEMSRC = libsbsgem/TGEMSBSBox.cxx \
         libsbsgem/TGEMSBSDBManager.cxx \
         libsbsgem/TGEMSBSGEMChamber.cxx \
         libsbsgem/TGEMSBSGEMHit.cxx \
         libsbsgem/TGEMSBSGEMPlane.cxx \
         libsbsgem/TGEMSBSGEMSimHitData.cxx \
         libsbsgem/TGEMSBSSimAuxi.cxx \
         libsbsgem/TGEMSBSSimDigitization.cxx \
         libsbsgem/TGEMSBSSpec.cxx
         #libsbsgem/TGEMSBSSimFile.cxx \
         #libsbsgem/TGEMSBSSimEvent.cxx \
         #libsbsgem/TGEMSBSGeant4File.cxx
         #libsbsgem/TGEMSBSSimDecoder.cxx

SRC   = src/g4sbs_tree.cxx \
        src/g4sbs_data.cxx \
        src/TSBSDBManager.cxx \
        src/TSBSGeant4File.cxx \
        src/TSBSSimAuxi.cxx \
        src/TSBSSimDecoder.cxx \
        src/TSBSSimEvent.cxx \
        src/TSBSSimFile.cxx \
        src/TSBSSimDigitizer.cxx \
        src/TSBSSimDetector.cxx \
        src/TSBSSimData.cxx \
	src/TSBSSimCher.cxx \
        src/TSBSSimECal.cxx \
        src/TSBSSimHCal.cxx \
        src/TSBSSimADC.cxx \
        src/TSBSSimTDC.cxx \
        src/TSBSSimMPD.cxx \
	src/TSBSSimScint.cxx \
        src/TSBSSpec.cxx \
        src/TSBSSimDataEncoder.cxx \
        src/TSBSSimGEM.cxx
SRC += $(GEMSRC)


OBJS	= $(SRC:.cxx=.$(ObjSuf)) $(DICT).o
HDR	= $(SRC:.cxx=.h) src/Linkdef.h
ROHDR	= $(SRC:.cxx=.h) src/Linkdef.h

LIBSBSDIG	= libsbsdig.so

PROGRAMS	= $(LIBSBSDIG)

all:	$(PROGRAMS)

$(LIBSBSDIG):	$(OBJS)
	$(LD) $(SOFLAGS) $^ -o $@ $(LDFLAGS) 

clean:
	@rm -f $(OBJS) $(PROGRAMS) *dict.*

$(DICT).cxx: $(ROHDR) 
	$(ROOTCINT) -f $@ -c $(SOLINCLUDE) $^ 

%.$(ObjSuf): %.$(SrcSuf)
	$(CXX)  $(CXXFLAGS) -c -o $@ $<
