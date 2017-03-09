MODULE CALC_FORCE_DEM_MODULE

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_FORCE_DEM                                          !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: Calculate contact force and torque on particle from        !
!           particle-particle and particle-wall collisions. Treats     !
!           wall interaction also as a two-particle interaction but    !
!           accounting for the wall properties                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_FORCE_DEM(particle_phase, des_radius, des_pos_new, des_vel_new, omega_new, pairs, pair_count, fc, tow)

         USE cfrelvel_module, only: cfrelvel
         USE discretelement, only: des_coll_model_enum, dtsolid
         USE discretelement, only: des_etan, des_etat, hert_kt, hert_kn
         USE discretelement, only: s_time, des_crossprdct
         USE discretelement, only: kn, kt, mew, hertzian

         USE drag_gs_des1_module, only: drag_gs_des
         USE error_manager, only: init_err_msg, flush_err_msg, err_msg, ival

         IMPLICIT NONE

      integer, intent(in) :: particle_phase(:)

      real(c_real), intent(in) :: des_radius(:)
      real(c_real), intent(in) :: des_pos_new(:,:)
      real(c_real), intent(in) :: des_vel_new(:,:)
      real(c_real), intent(in) :: omega_new(:,:)
      integer, intent(in) :: pairs(:,:)
      integer, intent(in) :: pair_count
      real(c_real), intent(inout) :: fc(:,:), tow(:,:)

! Local variables
!---------------------------------------------------------------------//
! percent of particle radius when excess overlap will be flagged
      real(c_real), PARAMETER :: flag_overlap = 0.20d0
! particle no. indices
      INTEGER :: I, LL, cc
! the overlap occuring between particle-particle or particle-wall
! collision in the normal direction
      real(c_real) :: OVERLAP_N, OVERLAP_T(3)
! square root of the overlap
      real(c_real) :: SQRT_OVERLAP
! distance vector between two particle centers or between a particle
! center and wall when the two surfaces are just at contact (i.e. no
! overlap)
      real(c_real) :: R_LM,DIST_CI,DIST_CL
! the normal and tangential components of the translational relative
! velocity
      real(c_real) :: V_REL_TRANS_NORM, rad
! distance vector between two particle centers or between a particle
! center and wall at current and previous time steps
      real(c_real) :: DIST(3), NORMAL(3), DIST_MAG, POS(3)
! tangent to the plane of contact at current time step
      real(c_real) :: VREL_T(3)
! normal and tangential forces
      real(c_real) :: FN(3), FT(3)
! temporary storage of force
      real(c_real) :: FC_TMP(3)
! temporary storage of force for torque
      real(c_real) :: TOW_FORCE(3)
! temporary storage of torque
      real(c_real) :: TOW_TMP(3,2)

! store solids phase index of particle (i.e. particle_phase(np))
      INTEGER :: PHASEI, PHASELL
! local values used spring constants and damping coefficients
      real(c_real) :: ETAN_DES, ETAT_DES
      real(c_real) :: KN_DES, KT_DES

      LOGICAL, PARAMETER :: report_excess_overlap = .FALSE.

      real(c_real) :: FNMD, MAG_OVERLAP_T, TANGENT(3)

!-----------------------------------------------

! Check particle LL neighbor contacts
!---------------------------------------------------------------------//

      DO cc = 1, pair_count

         ll = pairs(cc, 1)
         i = pairs(cc, 2)

         pos = DES_POS_NEW(LL,:)
         rad = DES_RADIUS(LL)

            R_LM = rad + DES_RADIUS(I)
            DIST(:) = DES_POS_NEW(I,:) - POS(:)
            DIST_MAG = dot_product(DIST,DIST)

            IF(abs(DIST_MAG) < epsilon(dist_mag)) THEN
               WRITE(*,8550) LL, I
               STOP "division by zero"
 8550 FORMAT('distance between particles is zero:',2(2x,I10))
            ENDIF

            DIST_MAG = SQRT(DIST_MAG)
            NORMAL(:)= DIST(:)/DIST_MAG

! Calcuate the normal overlap
            OVERLAP_N = R_LM-DIST_MAG
            IF(REPORT_EXCESS_OVERLAP) CALL PRINT_EXCESS_OVERLAP

! Calculate the components of translational relative velocity for a
! contacting particle pair and the tangent to the plane of contact
            CALL CFRELVEL(LL, I, V_REL_TRANS_NORM, VREL_T,            &
               NORMAL(:), DIST_MAG, DES_VEL_NEW, DES_RADIUS, OMEGA_NEW)

            phaseLL = particle_phase(LL)
            phaseI = particle_phase(I)

! Hertz spring-dashpot contact model
            IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN) THEN
               sqrt_overlap = SQRT(OVERLAP_N)
               KN_DES = hert_kn(phaseLL,phaseI)*sqrt_overlap
               KT_DES = hert_kt(phaseLL,phaseI)*sqrt_overlap
               sqrt_overlap = SQRT(sqrt_overlap)
               ETAN_DES = DES_ETAN(phaseLL,phaseI)*sqrt_overlap
               ETAT_DES = DES_ETAT(phaseLL,phaseI)*sqrt_overlap

! Linear spring-dashpot contact model
            ELSE
               KN_DES = KN
               KT_DES = KT
               ETAN_DES = DES_ETAN(phaseLL,phaseI)
               ETAT_DES = DES_ETAT(phaseLL,phaseI)
            ENDIF

! Calculate the normal contact force
            FN(:) =  -(KN_DES * OVERLAP_N * NORMAL(:) + &
               ETAN_DES * V_REL_TRANS_NORM * NORMAL(:))

! Calcuate the tangential overlap
            OVERLAP_T(:) = DTSOLID*VREL_T(:)
            MAG_OVERLAP_T = sqrt(dot_product(OVERLAP_T,OVERLAP_T))

! Calculate the tangential contact force.
            IF(MAG_OVERLAP_T > 0.0) THEN
! Max force before the on set of frictional slip.
               FNMD = MEW*sqrt(dot_product(FN,FN))
! Direction of tangential force.
               TANGENT = OVERLAP_T/MAG_OVERLAP_T
! Frictional slip
               FT = -FNMD * TANGENT
            ELSE
               FT = 0.0
            ENDIF

! calculate the distance from the particles' centers to the contact point,
! which is taken as the radical line
! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2
            DIST_CL = DIST_MAG/2.d0 + (DES_RADIUS(LL)**2 - &
               DES_RADIUS(I)**2)/(2.d0*DIST_MAG)

            DIST_CI = DIST_MAG - DIST_CL

            TOW_force(:) = DES_CROSSPRDCT(NORMAL(:), FT(:))
            TOW_TMP(:,1) = DIST_CL*TOW_force(:)
            TOW_TMP(:,2) = DIST_CI*TOW_force(:)

! Calculate the total force FC of a collision pair
! total contact force
            FC_TMP(:) = FN(:) + FT(:)

            FC(LL,:) = FC(LL,:) + FC_TMP(:)

            FC(I,1) = FC(I,1) - FC_TMP(1)
            FC(I,2) = FC(I,2) - FC_TMP(2)
            FC(I,3) = FC(I,3) - FC_TMP(3)

! for each particle the signs of norm and ft both flip, so add the same torque
            TOW(LL,:) = TOW(LL,:) + TOW_TMP(:,1)

            TOW(I,1)  = TOW(I,1)  + TOW_TMP(1,2)
            TOW(I,2)  = TOW(I,2)  + TOW_TMP(2,2)
            TOW(I,3)  = TOW(I,3)  + TOW_TMP(3,2)

      ENDDO

      RETURN

      contains

        include 'functions.inc'

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: print_excess_overlap                                    !
!                                                                      !
!  Purpose: Print overlap warning messages.                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE PRINT_EXCESS_OVERLAP

      IF(OVERLAP_N > flag_overlap*DES_RADIUS(LL) .OR.                  &
         OVERLAP_N > flag_overlap*DES_RADIUS(I)) THEN

         WRITE(ERR_MSG,1000) trim(iVAL(LL)), trim(iVAL(I)), S_TIME,    &
            DES_RADIUS(LL), DES_RADIUS(I), OVERLAP_N

         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

 1000 FORMAT('WARNING: Excessive overplay detected between ',          &
         'particles ',A,' and ',/A,' at time ',g11.4,'.',/             &
         'RADII:  ',g11.4,' and ',g11.4,4x,'OVERLAP: ',g11.4)

      END SUBROUTINE PRINT_EXCESS_OVERLAP

    END SUBROUTINE CALC_FORCE_DEM

END MODULE CALC_FORCE_DEM_MODULE
