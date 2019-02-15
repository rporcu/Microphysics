# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME ?= ../../amrex

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ..

# compilation options
COMP    = gnu
USE_MPI = TRUE
USE_OMP = FALSE

USE_HYPRE = FALSE

DEBUG     = FALSE

PROFILE       = FALSE
TRACE_PROFILE = FALSE
COMM_PROFILE  = FALSE
TINY_PROFILE  = FALSE

# physics
DIM      = 3

# Checking for Make.local in current directory. Check the Make.local.template
# for more info.
ifneq ("$(wildcard Make.local)", "")
    $(info Loading Make.local (for ALGOIM and BLITZ home) ...)
    include ./Make.local
endif

include ./Make.mfix
