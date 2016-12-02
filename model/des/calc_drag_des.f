module calc_drag_des_module

      use comp_mean_fields_module, only: comp_mean_fields
      use compar, only:  istart3, iend3, jstart3, jend3, kstart3, kend3
      use discretelement, only: DES_CONTINUUM_COUPLED
      use discretelement, only: DES_EXPLICITLY_COUPLED
      use discretelement, only: DRAG_FC, FC, MAX_PIP
      use drag_gs_des0_module, only: drag_gs_des0, drag_gs_gas0
      use drag_gs_des1_module, only: drag_gs_des1, drag_gs_gas1
      use functions, only: IS_NORMAL
      use particle_filter, only: DES_INTERP_GARG
      use particle_filter, only: DES_INTERP_SCHEME_ENUM

  contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_DRAG_DES                                           !
!                                                                      !
!  Purpose: This subroutine is called from DES routines. It calls      !
!  functions that calcultate the drag force acting on particles. No    !
!  field variables are updated.                                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_DRAG_DES(ep_g,u_g,v_g,w_g,ro_g,mu_g,gradPg)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN   ) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: gradPg&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)

      INTEGER :: II

! Apply the drag force calculated by the gas phase.
      IF(DES_EXPLICITLY_COUPLED) THEN

         IF(DES_CONTINUUM_COUPLED) THEN
            DO II = 1, MAX_PIP
               IF(IS_NORMAL(II)) &
                  FC(II,:) = FC(II,:) + DRAG_FC(II,:)
            ENDDO
         ENDIF

      ELSE

! Calculate gas-solids drag force on particle
         IF(DES_CONTINUUM_COUPLED) THEN
            SELECT CASE(DES_INTERP_SCHEME_ENUM)
            CASE(DES_INTERP_GARG)
               CALL DRAG_GS_DES0(ep_g,u_g,v_g,w_g,ro_g,mu_g,gradPg)
            CASE DEFAULT
               CALL DRAG_GS_DES1(ep_g,u_g,v_g,w_g,ro_g,mu_g,gradPg)
            END SELECT
         ENDIF

      ENDIF

      RETURN
      END SUBROUTINE CALC_DRAG_DES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DES_DRAG_GS                                             !
!                                                                      !
!  Purpose: This subroutine is only called from the CONTINUUM side. It !
!  calls the correct routine for calculating the gas drag force.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_DRAG_DES_2FLUID(ep_g,u_g,v_g,w_g,ro_g,mu_g)


      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN   ) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

! Calculate gas-solids drag force.
      IF(DES_CONTINUUM_COUPLED) THEN
         SELECT CASE(DES_INTERP_SCHEME_ENUM)
         CASE(DES_INTERP_GARG) ; CALL DRAG_GS_GAS0(ep_g,u_g,v_g,w_g,ro_g,mu_g)
         CASE DEFAULT          ; CALL DRAG_GS_GAS1(ep_g,u_g,v_g,w_g,ro_g,mu_g)
         END SELECT
      ENDIF

      RETURN
      END SUBROUTINE CALC_DRAG_DES_2FLUID


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_DRAG_DES_EXPLICIT                                  !
!                                                                      !
!  Purpose: This subroutine is only called from the CONTINUUM side.    !
!  Moreover, it is only called once per time step so that the total    !
!  drag force is calculated explicitly for the fluid phase and DEM     !
!  particles. This is done to ensure momentum conservation and reduce  !
!  computational overhead.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_DRAG_DES_EXPLICIT(ep_g,u_g,v_g,w_g,ro_g,rop_g,mu_g)


      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

! Bin particles to the fluid grid.
      CALL PARTICLES_IN_CELL
! Calculate mean fields (EPg).
      CALL COMP_MEAN_FIELDS(ep_g,ro_g,rop_g)

! Calculate gas-solids drag force on particle
      IF(DES_CONTINUUM_COUPLED) CALL DRAG_GS_GAS1(ep_g, u_g, v_g, w_g, ro_g, mu_g)

      END SUBROUTINE CALC_DRAG_DES_EXPLICIT

end module calc_drag_des_module
