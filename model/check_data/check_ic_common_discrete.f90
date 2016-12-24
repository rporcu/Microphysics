!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_IC_COMMON_DISCRETE                                !
!  Author:   R.Garg                                   Date: 11-Mar-14  !
!                                                                      !
!  Purpose: check the initial conditions input section for DEM.        !
!     - ensure the first IC is defined over the entire domain with     !
!        ep_g = 1 when more than one IC has solids                     !
!     - ensure the ICs are non-overlapping                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_IC_COMMON_DISCRETE

! Simulation dimension (2D/3D)
      USE discretelement, only: DIMN
! Flag indicating that the IC region is defined.
      USE ic, only: IC_DEFINED
! IC Region gas volume fraction.
      USE ic, only: IC_EP_G

      USE ic, only: IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b, IC_Z_t

      USE param1, only: UNDEFINED_I, ZERO, ONE

! direction wise spans of the domain and grid spacing in each direction
      Use geometry, only: xlength, ylength, zlength


! Maximum number of IC regions
      USE param, only: DIMENSION_IC

      ! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar, ival, err_msg

      implicit none

      INTEGER :: ICV, ICV2, IDIM
      INTEGER :: COUNT_IC, COUNT_IC_WITH_SOLS
      INTEGER :: FIRST_DEF_IC

      DOUBLE PRECISION :: IC_ORIG(3), IC_END(3), IC2_ORIG(3) , IC2_END(3)
      DOUBLE PRECISION :: IC_MIN, IC_MAX, IC2_MIN, IC2_MAX , TOL_IC_REG

      LOGICAL :: SEP_AXIS, first_ic_ok

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_IC_COMMON_DISCRETE")

! First check if multiple IC regions are defined for non-zero solids volume
! fraction, then check if the first IC is specified over the whole domain with IC_EP_g = 1

      !total count of defined ICs
      COUNT_IC           = 0
      !total count of defined IC's with solids
      COUNT_IC_WITH_SOLS = 0
      FIRST_DEF_IC = UNDEFINED_I
      DO ICV = 1, DIMENSION_IC

         IF (IC_DEFINED(ICV)) THEN
            COUNT_IC = COUNT_IC + 1
            FIRST_DEF_IC = MIN(FIRST_DEF_IC, ICV)

            IF(IC_EP_G(ICV).LT.ONE) COUNT_IC_WITH_SOLS &
            = COUNT_IC_WITH_SOLS  + 1

         ENDIF ! if(ic_defined(icv))
      end DO

      IF(COUNT_IC_WITH_SOLS >= 1 .AND. &
         COUNT_IC > COUNT_IC_WITH_SOLS+1) THEN

! If the number of IC's with solids is greater than one, make sure the
! first IC spans the entire domain with voidage of one. This ensures
! that the entire domain has valid ICs defined.
         ICV = FIRST_DEF_IC
         FIRST_IC_OK = .FALSE.
         IF(IC_EP_G(ICV).EQ.ONE &
           .AND.IC_X_W(ICV).LE.ZERO.AND.IC_X_E(ICV).GE.XLENGTH         &
           .AND.IC_Y_S(ICV).LE.ZERO.AND.IC_Y_N(ICV).GE.YLENGTH)        &
            FIRST_IC_OK = .TRUE.

         IF (FIRST_IC_OK .AND. IC_Z_B(ICV) <= ZERO .AND. &
            IC_Z_T(ICV) >= ZLENGTH) FIRST_IC_OK = .TRUE.

         IF(.NOT.FIRST_IC_OK) THEN
            WRITE(ERR_MSG, 1003)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1003 FORMAT(' Error 1003: Particle seeding with more than one IC ',   &
         'region requires',/'that IC 1 span the entire domain and ',   &
         'have IC_EP_g(1) = 1.0.',/'Please correct the mfix.dat file.')

      ENDIF

! Check if the ICs are non-overlapping.
      TOL_IC_REG  = 1E-04
      ICVLOOP : DO ICV = 1, DIMENSION_IC

         IF(.NOT.IC_DEFINED(ICV)) CYCLE ICVLOOP
         IF(IC_EP_G(ICV) == 1.d0) CYCLE ICVLOOP
         IC_ORIG(1) = IC_X_W(ICV)
         IC_ORIG(2) = IC_Y_S(ICV)
         IC_ORIG(3) = IC_Z_B(ICV)
         IC_END(1)  = IC_X_E(ICV)
         IC_END(2)  = IC_Y_N(ICV)
         IC_END(3)  = IC_Z_T(ICV)
         ICVTWOLOOP : DO ICV2 = ICV+1, DIMENSION_IC

            IF(.NOT.IC_DEFINED(ICV2)) CYCLE ICVTWOLOOP
            IF(IC_EP_G(ICV2) == 1.0d0) CYCLE ICVTWOLOOP

            IC2_ORIG(1) = IC_X_W(ICV2)
            IC2_ORIG(2) = IC_Y_S(ICV2)
            IC2_ORIG(3) = IC_Z_B(ICV2)
            IC2_END(1)  = IC_X_E(ICV2)
            IC2_END(2)  = IC_Y_N(ICV2)
            IC2_END(3)  = IC_Z_T(ICV2)

            sep_axis  = .false.
            DO idim = 1, dimn

               ic_min = IC_ORIG(idim)
               ic_max = IC_END(idim)
               ic2_min = IC2_ORIG(idim)
               ic2_max = ic2_END(idim)

! Check for separating axis. If the separating axis exists, then the IC
! regions can't overlap generally equality implies lack of sep_axis,
! and thus, overlapping. However, doing so will flag all IC's as
! overlapping since IC's have to share common edges. So here the
! equality is considered as existence of a separating axis, and hence,
! no overlap equality is also considered as separating axis which is
               if ((ic_min .ge. ic2_max)  .or. (ic_max .le. ic2_min) ) then
                  sep_axis = .true.
                  exit
               endif
            end DO

! Implies the IC regions could not find a separating axis and are
! therefore overlapping.
            IF(.NOT.sep_axis) THEN
               WRITE(ERR_MSG, 1004) ICV, ICV2
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF

 1004 FORMAT('Error 1004: Overlapping IC regions with nonzero solids ',&
         'volume',/'fraction detected. This is not supported for ',    &
         'discrete solids.',2/'Overlapping ICs: ',2(2x,I4),2/,         &
         'Please correct the mfix.dat file.')

         end DO ICVTWOLOOP
      end DO ICVLOOP



      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_IC_COMMON_DISCRETE
