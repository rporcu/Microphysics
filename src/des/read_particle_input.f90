MODULE READ_PAR_INPUT_MODULE

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! Subroutine: READ_PAR_INPUT                                           !
!                                                                      !
! Purpose: Read the particle input and broadcasts the particle data to !
! respective processors.                                               !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE READ_PAR_INPUT(particle_state, &
         des_radius, ro_sol, &
         des_pos_new, des_vel_new)

      use compar, only: myPE, pe_io
      use discretelement, only: normal_particle, particles
      use error_manager, only: err_msg, ival, init_err_msg, flush_err_msg, finl_err_msg
      use exit_mod, only: mfix_exit
      use run, only: bdist_io

      implicit none

      real(c_real), DIMENSION(:), INTENT(OUT) :: des_radius, ro_sol
      real(c_real), DIMENSION(:,:), INTENT(OUT) :: des_vel_new, des_pos_new
      integer, DIMENSION(:), INTENT(OUT) :: particle_state

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      integer :: k
! index of particle
      integer :: lcurpar
! local unit
      integer, PARAMETER :: lunit=10
! local filename
      character(255) lfilename
! IO Status:
      integer :: IOS
! Flag to indicate if file exists.
      logical :: lEXISTS
! Read dimension: 2D vs 3D data
      integer :: RDMN
!-----------------------------------------------


      CALL INIT_ERR_MSG("READ_PAR_INPUT")


      IOS = 0
      RDMN = 3

! Setup the file name based on distributed or serial IO.
      lFILENAME= "particle_input.dat"

! Check the the file exists and open it.
      INQUIRE(FILE=lFILENAME, EXIST=lEXISTS)
      IF(.NOT.LEXISTS) THEN
         WRITE(ERR_MSG, 1100)
         CALL FLUSH_ERR_MSG
         IOS = 1
      ELSE
         OPEN(UNIT=lUNIT, FILE=lFILENAME, FORM="FORMATTED")
      ENDIF

! Collect the error message and quit.
      ! CALL GLOBAL_ALL_SUM(IOS)
      IF(IOS /= 0) CALL MFIX_EXIT(myPE)

 1100 FORMAT('Error 1100: FATAL - DEM particle input file not found!')

! Read the file
!----------------------------------------------------------------->>>
! In distributed IO the first line of the file will be number of
! particles in that processor
      DO lcurpar = 1,particles
         particle_state(lcurpar) = normal_particle
         read (lunit,*) (des_pos_new(lcurpar,k),k=1,RDMN),&
            des_radius(lcurpar), ro_sol(lcurpar),&
            (des_vel_new(lcurpar,k),k=1,RDMN)
      ENDDO

      IF(bDIST_IO .OR. myPE == PE_IO) CLOSE(lUNIT)

      CALL FINL_ERR_MSG()

      RETURN

      END SUBROUTINE READ_PAR_INPUT
END MODULE READ_PAR_INPUT_MODULE
