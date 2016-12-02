!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS                                                    !
!  Author: J.Musser                                   Date: 17-Nov-14  !
!                                                                      !
!  REF: Di Renzo, A. and Di Maio F.P. "Comparison of contact-force     !
!       models for the simulation of collisions in DEM-based granular  !
!       flow codes," Chemical Engineering Science, 59(3), pg 525-541.  !
!                                                                      !
!  REF: Kharaz, A.H., Gorham, D.A., and Salman, A.D. "An experimental  !
!       study of the elastic rebound of spheres," Powder Technology,   !
!       120(3), pg 281-291.                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE usr

      Use param
      Use param1

! a dummy variable listed in usrnlst.inc
      DOUBLE PRECISION DUMMY_DP

! Initial collision angles
      DOUBLE PRECISION :: INIT_ANGLE(100)
! Initial Tangential velocity
      DOUBLE PRECISION :: INIT_Vel_T(100)

      contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE COLLECT_DEM05_DATA(pDATA, pTMP)

      use discretelement, only: iGlobal_ID
      use discretelement, only: MAX_PIP

      use discretelement, only: NORMAL_PARTICLE, PARTICLE_STATE

      double precision, intent(in) :: pDATA(:)
      double precision, intent(out) :: pTMP(62)

      integer :: np
      double precision :: lTMP(62)

! Map local proc data to global array
      lTMP = 0.0d0
      do np=1,max_pip
         if(normal_particle==particle_state(np) .and. iGlobal_ID(np) <=62) &
            lTMP(np) = pDATA(np)
      enddo

! Collect data on root proc
      pTMP = lTMP
      ! call global_sum(lTMP, pTMP)
      return

      END SUBROUTINE COLLECT_DEM05_DATA


      END MODULE usr
