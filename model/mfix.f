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
      USE cutcell
      USE dashboard
      USE discretelement
      USE fldvar
      USE functions
      USE funits
      USE machine
      USE mpi_utility
      USE output
      USE param
      USE param1
      USE quadric
      USE run
      USE time_cpu

      USE vtk, only : WRITE_VTK_FILES

      use error_manager

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Final value of CPU time.
      DOUBLE PRECISION :: CPU1
! time used for computations.
      DOUBLE PRECISION :: CPUTIME_USED, WALLTIME_USED
! CPU time unit.
      CHARACTER(LEN=4) :: TUNIT
! Save TIME in input file for RESTART_2
      DOUBLE PRECISION :: TIME_SAVE
! Temporary storage for DT
      DOUBLE PRECISION :: DT_tmp
! loop counter
      INTEGER :: L
! DISTIO variable for specifying the mfix version
      CHARACTER(LEN=512) :: version
! environment variable
!$      CHARACTER(LEN=512) :: omp_num_threads
!$      INTEGER :: length
!$      INTEGER :: status

!$      INTEGER num_threads, threads_specified, omp_id
!$      INTEGER omp_get_num_threads
!$      INTEGER omp_get_thread_num

! C Function
      INTERFACE
         SUBROUTINE INIT_CMD_SOCKET(port) BIND ( C )
           use, INTRINSIC :: iso_c_binding
           CHARACTER(KIND=C_CHAR), INTENT(IN) :: port(*)
         END SUBROUTINE INIT_CMD_SOCKET
         SUBROUTINE INIT_LOG_SOCKET(port) BIND ( C )
           use, INTRINSIC :: iso_c_binding
           CHARACTER(KIND=C_CHAR), INTENT(IN) :: port(*)
         END SUBROUTINE INIT_LOG_SOCKET
      END INTERFACE

!-----------------------------------------------

! DISTIO
! If you change the value below in this subroutine, you must also
! change it in write_res0.f and the value should also be consistent
! with the check in read_res0
      version = 'RES = 01.6'

! Invoke MPI initialization routines and get rank info.
      CALL PARALLEL_INIT
      CALL GEN_LOG_BASENAME

! we want only PE_IO to write out common error messages
      DMP_LOG = (myPE == PE_IO)

! set the version.release of the software
      ID_VERSION = '2015-2'


! specify the number of processors to be used
!$        call get_environment_variable("OMP_NUM_THREADS",omp_num_threads,length,status, .true.)
!$      if (status.eq.0 .and. length.ne.0) then
!$        read(omp_num_threads,*) threads_specified
!$      else
!$        WRITE(*,'(A,$)') 'Enter the number of threads to be used for SMP: '
!$        READ(*,*) threads_specified
!$      endif

!$      call omp_set_num_threads(threads_specified)

! Find the number of processors used
!$omp  parallel
!$      num_threads = omp_get_num_threads()
!$      omp_id = omp_get_thread_num()
!$      if(omp_id.eq.0) Write(*,*)' Number of threads used for SMP = ',  num_threads
!$omp  end parallel


! Set machine dependent constants
      CALL MACHINE_CONS

! Get the date and time. They give the unique run_id in binary output
! files
      CALL GET_RUN_ID

! AEOLUS: stop trigger mechanism to terminate MFIX normally before batch
! queue terminates. timestep at the beginning of execution
      CALL CPU_TIME (CPU00)
      WALL0 = WALL_TIME()

! Read input data, check data, do computations for IC and BC locations
! and flows, and set geometry parameters such as X, X_E, DToDX, etc.
      CALL GET_DATA

! Write the initial part of the standard output file
      CALL WRITE_OUT0
      IF(.NOT.CARTESIAN_GRID)  CALL WRITE_FLAGS

! Write the initial part of the special output file(s)
      CALL WRITE_USR0

!$    CALL START_LOG
!$    IF(DMP_LOG)WRITE (UNIT_LOG, *) ' '
!$    IF(DMP_LOG)WRITE (UNIT_LOG, *) ' Number of processors used = ', threads_specified
!$    IF(DMP_LOG)WRITE (UNIT_LOG, *) ' '
!$    CALL END_LOG

!  setup for PC quickwin application
      CALL PC_QUICKWIN


      CALL INIT_ERR_MSG('MFIX')


  101 CONTINUE


      DT_TMP = DT
      SELECT CASE (TRIM(RUN_TYPE))

      CASE ('NEW')
! Write the initial part of the restart files
         CALL WRITE_RES0
         DO L = 1, N_SPX
            CALL WRITE_SPX0 (L, 0)
         ENDDO

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

! Writing the RES1 and SPX1 can only be done here when re-indexing is turned off
! This will be done after the cell re-indexing is done later in this file.
! This allows restarting independently of the re-indexing setting between
! the previous and current run.
         IF(.NOT.RE_INDEXING) THEN
            CALL WRITE_RES1
            DO L = 1, N_SPX
               CALL WRITE_SPX0 (L, 0)
               CALL WRITE_SPX1 (L, 0)
            END DO
         ENDIF

      CASE DEFAULT
         CALL START_LOG
         IF(DMP_LOG)WRITE (UNIT_LOG, *) &
            ' MFIX: Do not know how to process'
         IF(DMP_LOG)WRITE (UNIT_LOG, *) ' RUN_TYPE in data file'
         CALL END_LOG
         call mfix_exit(myPE)

      END SELECT

      call MPI_Barrier(MPI_COMM_WORLD,mpierr)

      IF (DT_TMP /= UNDEFINED) THEN
         DT = MAX(DT_MIN,MIN(DT_MAX,DT))
      ELSE
         DT = DT_TMP
      ENDIF

! Set arrays for computing indices. A secondary call is made
! after cut cell-preprocessing to update array indices.
      IF(CARTESIAN_GRID) THEN
         CALL SET_INCREMENTS
      ENDIF


! Set the flags for wall surfaces impermeable and identify flow
! boundaries using FLAG_E, FLAG_N, and FLAG_T
      CALL SET_FLAGS1

!  Update flags for Cartesian_GRID.
      IF(CARTESIAN_GRID) CALL CHECK_BC_FLAGS

! Calculate cell volumes and face areas
      IF(.NOT.CARTESIAN_GRID) CALL SET_GEOMETRY1

! Find corner cells and set their face areas to zero
      IF(.NOT.CARTESIAN_GRID)  THEN
         CALL GET_CORNER_CELLS()
      ELSE
         IF (SET_CORNER_CELLS)  CALL GET_CORNER_CELLS ()
      ENDIF

! Set constant physical properties
      CALL SET_CONSTPROP

! Set initial conditions
      CALL SET_IC

! Set point sources.
      CALL SET_PS

! Set boundary conditions
      CALL ZERO_NORM_VEL
      CALL SET_BC0

! JFD: cartesian grid implementation
      IF(CARTESIAN_GRID) CALL CG_SET_BC0

! Set the pressure field for a fluidized bed
      IF (RUN_TYPE == 'NEW') CALL SET_FLUIDBED_P

! Initialize densities.
      IF (RUN_TYPE == 'NEW') CALL SET_RO_G

! Initialize time dependent boundary conditions
      CALL SET_BC1

! Check the field variable data and report errors.
      IF(.NOT.CARTESIAN_GRID)  CALL CHECK_DATA_20

! Setup VTK data for regular (no cut cells) grid
      IF(.NOT.CARTESIAN_GRID.AND.WRITE_VTK_FILES) CALL SETUP_VTK_NO_CUTCELL

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

! JFD: cartesian grid implementation
      IF(WRITE_DASHBOARD) THEN
         IF(DT>=DT_MIN) THEN
            RUN_STATUS = 'Complete.'
         ELSE
            RUN_STATUS = 'DT < DT_MIN.  Recovery not possible!'
         ENDIF
         CALL GET_TUNIT(CPUTIME_USED,TUNIT)
         CALL UPDATE_DASHBOARD(0,CPUTIME_USED,TUNIT)
      ENDIF
      IF(CARTESIAN_GRID)  CALL CLOSE_CUT_CELL_FILES

! Finalize and terminate MPI
      call parallel_fin

      CALL FINL_ERR_MSG

      STOP

 1000 FORMAT(/1X,'MFIX ',A,' Simulation:'/)

 1010 FORMAT('Message 1010: Read in data from .RES file for TIME = ',&
         G12.5,/'Time step number (NSTEP) =',I7)

      END PROGRAM MFIX

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GEN_LOG_BASENAME                                        !
!  Author: Aytekin Gel                                Date: 19-SEP-03  !
!                                                                      !
!  Purpose: Generate the file base for DMP logs.                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GEN_LOG_BASENAME

      use compar, only: myPE
      use compar, only: fbname

      implicit none

! Variables for generating file basename with processor id
      INTEGER :: i1, i10, i100, i1000, i10000

! PAR_I/O Generate file basename for LOG files
      i10000 = int(myPE/10000)
      i1000  = int((myPE-i10000*10000)/1000)
      i100   = int((myPE-i10000*10000-i1000*1000)/100)
      i10    = int((myPE-i10000*10000-i1000*1000-i100*100)/10)
      i1     = int((myPE-i10000*10000-i1000*1000-i100*100-i10*10)/1)

      i10000 = i10000 + 48
      i1000  = i1000  + 48
      i100   = i100   + 48
      i10    = i10    + 48
      i1     = i1     + 48

      fbname=char(i10000)//char(i1000)//char(i100)//char(i10)//char(i1)

      RETURN
      END SUBROUTINE GEN_LOG_BASENAME
