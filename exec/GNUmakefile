# define the location of the Microphysics top directory
MICROPHYSICS_HOME  := ../
AMREX_HOME  := ../subprojects/amrex

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP       = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE
USE_CUDA   = FALSE
CUDA_VERBOSE = FALSE
USE_REACT = TRUE

EBASE = main

BL_NO_FORT = TRUE

INTEGRATOR_DIR =  VODE

ifeq ($(USE_CUDA), TRUE)
  INTEGRATOR_DIR := VODE
endif

EXTERN_SEARCH += . ..

Bpack   := ./Make.package
Blocs   := .

include Make.Microphysics
