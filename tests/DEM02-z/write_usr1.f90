!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_USR1 (L)                                         C
!  Purpose: Write user-defined output                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine write_usr1(l, time, dt, max_pip, des_pos_new, des_vel_new, omega_new)

      use amrex_fort_module, only : c_real => amrex_real

      IMPLICIT NONE

      integer,      intent(in   ) :: l, max_pip
      real(c_real), intent(in   ) :: time, dt
      real(c_real), intent(in   ) :: des_pos_new(max_pip,3)
      real(c_real), intent(in   ) :: des_vel_new(max_pip,3)
      real(c_real), intent(in   ) :: omega_new(max_pip,3)

      SELECT CASE(L)
      CASE(1); CALL WRITE_DES_OUT()
      END SELECT

      RETURN
      END SUBROUTINE WRITE_USR1


!......................................................................!
!  Subroutine: WRITE_DES_Out                                           !
!                                                                      !
!  Purpose: Calculate the position and velocity (1D Y-axis) of a free  !
!  falling particle. Compare the results to the MFIX-DEM solultion.    !
!                                                                      !
!  Author: J.Musser                                   Date:  Jan-13    !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!......................................................................!
      SUBROUTINE WRITE_DES_Out

      use amrex_fort_module, only : c_real => amrex_real

      use compar, only: mype, pe_io
      Use usr, only: bounce_count, max_height_hs, max_height

      IMPLICIT NONE


! Local variables
!---------------------------------------------------------------------//
! file unit for heat transfer data
      INTEGER, PARAMETER :: uHeight = 2030
! Absolute relative percent differnece with hard sphere model.
      real(c_real) :: apDiff
! Max height of hard-sphere bounce
      real(c_real) :: maxHS
! Loop counter
      integer :: k
! Save the last location of the loop counter
      integer, save :: last_k = 0

      IF(myPE /= PE_IO) RETURN

! Open the files.
      OPEN(UNIT=uHeight,FILE='POST_HEIGHT.dat', &
         POSITION="APPEND",STATUS='OLD')


! Loop through the bounces and output the results.
      do k=last_k, BOUNCE_COUNT

! Calculate the hard sphere max bounce height.
         maxHS = MAX_HEIGHT_HS(k)

! Calculate the absolute relative percent difference.
         apDiff = (abs(MAX_HEIGHT(k) - maxHS)/abs(maxHS))*1.0d2

! Write the results to a file.
         WRITE(uHeight,"(3x,I2,5X,F15.8,2(3x,F15.8))") k, &
            maxHS, MAX_HEIGHT(k), apDiff

      enddo
      last_k = k


      CLOSE(uHeight)

      RETURN
      END SUBROUTINE WRITE_DES_Out
