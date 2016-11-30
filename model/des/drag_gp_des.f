!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DES_DRAG_GP                                             C
!  Purpose: Calculate the gas-particle drag coefficient using          C
!           the gas velocity interpolated to the particle position     C
!           and the particle velocity.                                 C
!           Invoked from des_drag_gs and calc_des_drag_gs              C
!                                                                      C
!  Comments: The BVK drag model and all drag models with the           C
!            polydisperse correction factor (i.e., suffix _PCF)        C
!            require an average particle diameter. This has been       C
!            loosely defined for discrete particles based on their     C
!            solids phase                                              C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DES_DRAG_GP(NP, PARTICLE_VEL, FLUID_VEL, EPg)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE constant
      USE discretelement
      USE drag
      USE fldvar
      USE functions
      USE geometry
      USE param
      USE param1
      USE physprop
      USE run
      USE funits

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! particle number id.
      INTEGER , INTENT(IN) :: NP
! particle velocity
      DOUBLE PRECISION, INTENT(IN) :: PARTICLE_VEL(3)
! fluid velocity interpolated to particle position
      DOUBLE PRECISION, INTENT(IN) :: FLUID_VEL(3)
! Gas phase volume fraction.
      DOUBLE PREcISION, INTENT(IN) :: EPg

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices, associated with current particle
      INTEGER :: IJK, I, J, K
! solids phase index, associated with current particle
      INTEGER :: M
! Slip velocity and its magnitude
      DOUBLE PRECISION :: VSLP(3), VREL
! gas laminar viscosity redefined here to set viscosity at pressure
! boundaries
      DOUBLE PRECISION :: Mu
! drag coefficient
      DOUBLE PRECISION :: DgA
! indices of solids phases (continuous, discrete)
      INTEGER :: lM
! correction factors for implementing polydisperse drag model
! proposed by van der Hoef et al. (2005)
      DOUBLE PRECISION :: F_cor, tSUM
! average particle diameter in polydisperse systems
      DOUBLE PRECISION :: DPA
! diameter ratio in polydisperse systems
      DOUBLE PRECISION :: Y_i
! total solids volume fraction
      DOUBLE PRECISION :: PHIS
! aliases for void fraction, gas density, gas bulk density,
! solids volume fraction, particle diameter, particle density
      DOUBLE PRECISION :: ROg, ROPg, DPM
!-----------------------------------------------


! values based on current particle
      IJK = PIJK(NP,4)
      I = PIJK(NP,1)
      J = PIJK(NP,2)
      K = PIJK(NP,3)
! solids phase index of current particle
      M = PIJK(NP,5)
! Gas material and bulk densities
      ROg = RO_G(I,J,K)
      ROPg = RO_G(I,J,K) * EPg
! Laminar viscosity.
      Mu = MU_G(I,J,K)
! Slip velocity and its magnitude
      VSLP = FLUID_VEL - PARTICLE_VEL
      VREL = SQRT(dot_product(VSLP, VSLP))
! assign variables for short dummy arguments
      DPM = 2.0d0*DES_RADIUS(NP)
! Total solids volume fraction.
      PHIS = ONE - EPg

! determine the drag coefficient
      SELECT CASE(DRAG_TYPE_ENUM)

      CASE (SYAM_OBRIEN)
         CALL DRAG_SYAM_OBRIEN(DgA,EPG,Mu,ROg,VREL,DPM)

      CASE (GIDASPOW)
         CALL DRAG_GIDASPOW(DgA,EPg,Mu,ROg,ROPg,VREL,DPM)

      CASE (GIDASPOW_BLEND)
         CALL DRAG_GIDASPOW_BLEND(DgA,EPg,Mu,ROg,ROPg,VREL,DPM)

      CASE (WEN_YU)
         CALL DRAG_WEN_YU(DgA,EPg,Mu,ROPg,VREL,DPM)

      CASE (KOCH_HILL)
         CALL DRAG_KOCH_HILL(DgA,EPg,Mu,ROPg,VREL,DPM,DPM,phis)

      CASE (USER_DRAG)
         CALL DRAG_USR(IJK,NP,DgA,EPg,Mu,ROg,VREL,DPM,RO_S0(M), &
            FLUID_VEL(1), FLUID_VEL(2), FLUID_VEL(3))

      CASE DEFAULT

! calculate the average particle diameter and particle ratio
         tSUM = ZERO
         DO lM = 1,MMAX
            IF(PHIS > ZERO) THEN
               tSUM = tSUM + DES_ROP_S(IJK,lM) / &
                  (PHIS*RO_S0(lM)*D_p0(lM))
             ELSE
               tSUM = tSUM + ONE/D_p0(lM)
             ENDIF
         ENDDO

         DPA = ONE / tSUM
         Y_i = DPM * tSUM

         SELECT CASE(DRAG_TYPE_ENUM)


         CASE (GIDASPOW_PCF)
            CALL DRAG_GIDASPOW(DgA,EPg,Mu,ROg,ROPg,VREL,DPA)
         CASE (GIDASPOW_BLEND_PCF)
            CALL DRAG_GIDASPOW_BLEND(DgA,EPg,Mu,ROg,ROPg,VREL,DPA)
         CASE (WEN_YU_PCF)
            CALL DRAG_WEN_YU(DgA,EPg,Mu,ROPg,VREL,DPA)
         CASE (KOCH_HILL_PCF)
            CALL DRAG_KOCH_HILL(DgA,EPg,Mu,ROPg,VREL,DPM,DPA,phis)
         CASE (BVK)
            CALL DRAG_BVK(DgA,EPg,Mu,ROPg,VREL,DPM,DPA,phis)

         CASE DEFAULT
            IF(DMP_LOG) WRITE (*, '(A,A)') &
               'Unknown DRAG_TYPE: ', DRAG_TYPE
            WRITE (UNIT_LOG, '(A,A)') 'Unknown DRAG_TYPE: ', DRAG_TYPE
            CALL mfix_exit(myPE)
         END SELECT   ! end selection of drag_type


! Modify drag coefficient to account for possible corrections and for
! differences between Model B and Model A; see erratum Beetstra (2007)
         F_cor = Y_i
         DgA = ONE/(Y_i*Y_i) * DgA * F_cor

      END SELECT   ! end selection of drag_type

! Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
      F_gp(NP) = DgA * PVOL(NP)


      RETURN
      END SUBROUTINE DES_DRAG_GP
