TOP=../..

include $(TOP)/configure/CONFIG
#----------------------------------------
#  ADD MACRO DEFINITIONS AFTER THIS LINE

#=============================
# Build a support module called <APPNAME>

LIBRARY_IOC += tptlib

# <APPNAME>.dbd will be created and installed
DBD += tptlib.dbd

# tptlib.dbd will be made up from these files:
tptlib_DBD += base.dbd
# List local .dbd files here
# tptlib_DBD += xxx.dbd

# Include dbd files for external support modules:
# tptlib_DBD += genSubRecord.dbd

# specify all source files to be compiled
# including sequencer (.st) source files
#tptlib_SRCS += xxx.c
#tptlib_SRCS += xxx.st

tptlib_SRCS += addtrm.c
tptlib_SRCS += length.c
tptlib_SRCS += minit.c
tptlib_SRCS += ptappl.c
tptlib_SRCS += pterm.c
tptlib_SRCS += ptermg.c
tptlib_SRCS += pterml.c
tptlib_SRCS += pterms.c
tptlib_SRCS += ptpre.c
tptlib_SRCS += pup.c
tptlib_SRCS += scomp.c 
tptlib_SRCS += term.c
tptlib_SRCS += trmgen.c
tptlib_SRCS += trmloc.c
tptlib_SRCS += trmstd.c

INC += tptlib.h
#=============================

include $(TOP)/configure/RULES
#----------------------------------------
#  ADD RULES AFTER THIS LINE

