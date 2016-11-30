!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MFIX                                                    !
!  Author: M. Syamlal                                 Date: 29-JAN-92  !
!                                                                      !
!  Purpose: The main module in the MFIX program                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!
!> \mainpage Multiphase Flow with Interphase eXchanges
!!
!! MFIX is a general-purpose computer code developed at the National
!! Energy Technology Laboratory, NETL, for describing the hydrodynamics,
!! heat transfer, and chemical reactions in fluid-solid systems.
!!
!! It has been used for describing bubbling and circulating fluidized
!! beds and spouted beds. MFiX calculations give transient data on the
!! three-dimensional distribution of pressure, velocity, temperature,
!! and species mass fractions. MFiX code is based on a generally
!! accepted set of multiphase flow equations. The code is used as a
!! "test-stand" for testing and developing multiphase flow constitutive
!!  equations.
!!
!! \section Notice
!! Neither the United States Government nor any agency thereof, nor any
!! of their employees, makes any warranty, expressed or implied, or
!! assumes any legal liability or responsibility for the accuracy,
!! completeness, or usefulness of any information, apparatus, product,
!! or process disclosed or represents that its use would not infringe
!! privately owned rights.
!!
!! * MFIX is provided without any user support for applications in the
!!   user's immediate organization. It should not be redistributed in
!!   whole or in part.
!!
!! * The use of MFIX is to be acknowledged in any published paper based
!!   on computations using this software by citing the MFIX theory
!!   manual. Some of the submodels are being developed by researchers
!!   outside of NETL. The use of such submodels is to be acknowledged
!!   by citing the appropriate papers of the developers of the submodels.
!!
!! * The authors would appreciate receiving any reports of bugs or other
!!   difficulties with the software, enhancements to the software, and
!!   accounts of practical applications of this software.
!!
!! \section Disclaimer
!! This report was prepared as an account of work sponsored by an agency
!! of the United States Government. Neither the United States Government
!! nor any agency thereof, nor any of their employees, makes any
!! warranty, express or implied, or assumes any legal liability or
!! responsibility for the accuracy, completeness, or usefulness of any
!! information, apparatus, product, or process disclosed, or represents
!! that its use would not infringe privately owned rights. Reference
!! herein to any specific commercial product, process, or service by
!! trade name, trademark, manufacturer, or otherwise does not
!! necessarily constitute or imply its endorsement, recommendation, or
!! favoring by the United States Government or any agency thereof. The
!! views and opinions of authors expressed herein do not necessarily
!! state or reflect those of the United States Government or any
!! agency thereof.

      PROGRAM MFIX

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE discretelement
      USE fldvar
      USE fld_const, only: ro_g0
      USE functions
      USE funits
      USE machine
      USE output
      USE param
      USE param1
      USE run
      USE time_cpu

      use error_manager

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Final value of CPU time.
      DOUBLE PRECISION :: CPU1
! time used for computations.
      DOUBLE PRECISION :: CPUTIME_USED, WALLTIME_USED
! Save TIME in input file for RESTART_2
      DOUBLE PRECISION :: TIME_SAVE
! Temporary storage for DT
      DOUBLE PRECISION :: DT_tmp

!-----------------------------------------------

! Invoke MPI initialization routines and get rank info.
      ! CALL PARALLEL_INIT

! AEOLUS: stop trigger mechanism to terminate MFIX normally before batch
! queue terminates. timestep at the beginning of execution
      CALL CPU_TIME (CPU00)
      WALL0 = WALL_TIME()

! Read input data, check data, do computations for IC and BC locations
! and flows, and set geometry parameters such as X, X_E, DToDX, etc.
      CALL GET_DATA

! Write the initial part of the standard output file
      CALL WRITE_OUT0
!     CALL WRITE_FLAGS

! Write the initial part of the special output file(s)
      CALL WRITE_USR0

      CALL INIT_ERR_MSG('MFIX')

      DT_TMP = DT
      SELECT CASE (TRIM(RUN_TYPE))

      CASE ('NEW')
! Write the initial part of the restart files
         CALL WRITE_RES0

       CASE ('RESTART_1')
! Read the time-dependent part of the restart file
         CALL READ_RES1
         WRITE(ERR_MSG, 1010) TIME, NSTEP
         CALL FLUSH_ERR_MSG()

      CASE ('RESTART_2')
         TIME_SAVE = TIME

         CALL READ_RES1
         TIME = TIME_SAVE

         WRITE(ERR_MSG, 1010) TIME, NSTEP
         CALL FLUSH_ERR_MSG()

         CALL WRITE_RES0

! Writing the RES1 can only be done here when re-indexing is turned off
! This will be done after the cell re-indexing is done later in this file.
! This allows restarting independently of the re-indexing setting between
! the previous and current run.
         CALL WRITE_RES1

      CASE DEFAULT
         IF(DMP_LOG)WRITE (UNIT_LOG, *) &
            ' MFIX: Do not know how to process'
         IF(DMP_LOG)WRITE (UNIT_LOG, *) ' RUN_TYPE in data file'
         call mfix_exit(myPE)

      END SELECT

#ifdef MPI
      CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif

      IF (DT_TMP /= UNDEFINED) THEN
         DT = MAX(DT_MIN,MIN(DT_MAX,DT))
      ELSE
         DT = DT_TMP
      ENDIF

! Set the flags for wall surfaces impermeable and identify flow
! boundaries using FLAG_E, FLAG_N, and FLAG_T
      CALL SET_FLAGS1

! Calculate cell volumes and face areas
      CALL SET_GEOMETRY1

! Find corner cells and set their face areas to zero
      CALL GET_CORNER_CELLS()

! Set constant physical properties
      CALL SET_CONSTPROP(ro_g,lambda_g,mu_g)

! Set initial conditions
      CALL SET_IC(ep_g, p_g, u_g, v_g, w_g)

! Set point sources.
      CALL SET_PS

! Set boundary conditions
      CALL ZERO_NORM_VEL(u_g,v_g,w_g)
      CALL SET_BC0(p_g,ep_g,u_g,v_g,w_g,ro_g0)

! Set the pressure field for a fluidized bed
      IF (RUN_TYPE == 'NEW') CALL SET_FLUIDBED_P(p_g, ep_g)

! Initialize densities.
      IF (RUN_TYPE == 'NEW') CALL SET_RO_G(ro_g,rop_g,p_g,ep_g)

! Initialize time dependent boundary conditions
      CALL SET_BC1(p_g, ep_g, ro_g, rop_g, u_g, v_g, w_g)

! Check the field variable data and report errors.
      CALL CHECK_DATA_20(ep_g,p_g,ro_g,rop_g,u_g,v_g,w_g)

      IF(DISCRETE_ELEMENT) CALL MAKE_ARRAYS_DES

! Set the inflow/outflow BCs for DEM solids
      IF(DEM_SOLIDS) CALL SET_BC_DEM

! Initializations for CPU time calculations in iterate
      CPUOS = 0.
      CALL CPU_TIME (CPU1)
      CPU_NLOG = CPU1
      TIME_NLOG = TIME - DT




! Find the solution of the equations from TIME to TSTOP at
! intervals of DT
      CALL TIME_MARCH

! Call user-defined subroutine after time-loop.
      IF (CALL_USR) CALL USR3

! Compute the CPU time and write it out in the .OUT file.
      CPUTIME_USED = 0.0d0
      WALLTIME_USED = WALL_TIME() - WALL0
      CALL WRITE_OUT3 (CPUTIME_USED, WALLTIME_USED, CPU_IO)

! Finalize and terminate MPI
      ! call parallel_fin

      CALL FINL_ERR_MSG

      STOP

 1010 FORMAT('Message 1010: Read in data from .RES file for TIME = ',&
         G12.5,/'Time step number (NSTEP) =',I7)

      END PROGRAM MFIX
