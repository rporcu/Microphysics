module calc_pg_grad_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_pg_grad                                            !
!  Purpose: Calculate cell centered pressure force exerted on the      !
!           particles in the cell by the gas/fluid phase               !
!           Note that gradPg is evaluated as -dp/dx                    !
!                                                                      !
!  Notes: This pressure force only needs to be calculated once during  !
!         the DEM loop (at the beginning) since the gas/fluid phase    !
!         is essentially static at that point (i.e., gas field is not  !
!         updated during DEM loop                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine calc_pg_grad(slo, shi, lo, hi, max_pip, &
                              p_g, gradPg,  particle_state, des_pos_new,&
                              pvol, drag_fc, dx, dy, dz, domlo, domhi)

      use calc_grad_des_module, only: calc_grad_des
      use discretelement, only: entering_particle, entering_ghost
      use discretelement, only: exiting_particle, exiting_ghost
      use discretelement, only: nonexistent

      ! Flags for cyclic BC with pressure drop
      use geometry, only: CYCLIC_X_PD, CYCLIC_Y_PD, CYCLIC_Z_PD

      ! Specified pressure drop
      use bc, only: DELP_X, DELP_Y, DELP_Z

      ! Domain length
      use geometry, only: XLENGTH, YLENGTH, ZLENGTH

      use discretelement, only: DES_EXPLICITLY_COUPLED

! Global Parameters:
!---------------------------------------------------------------------//
! Double precision values.
      use param1, only: ZERO

      implicit none

      integer, intent(in   ) :: slo(3),shi(3)
      integer, intent(in   ) ::  lo(3), hi(3)
      integer, intent(in   ) :: domlo(3), domhi(3)
      integer, intent(in   ) :: max_pip

      real(c_real), intent(in   ) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(out  ) :: gradpg&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

      integer     , intent(in   ) :: particle_state(max_pip)
      real(c_real), intent(in   ) :: pvol(max_pip)
      real(c_real), intent(in   ) :: des_pos_new(max_pip,3)
      real(c_real), intent(inout) :: drag_fc(max_pip,3)

      real(c_real), intent(in   ) :: dx, dy, dz


! Loop counters: Particle, fluid cell, neighbor cells
      integer :: NP, I, J, K
! mean pressure gradient for the case of periodic boundaries
      real(c_real) :: cPG(3)
! One over cell volume
      real(c_real) :: odx, ody, odz
!......................................................................!

! Calculate the gas phase pressure gradient. (dP/dx)
      call calc_grad_des(slo, shi, lo, hi, P_G, gradPg, dx, dy, dz, domlo, domhi)

! Add in cyclic BC pressure drop.
      cPG(1) = merge(DELP_X/XLENGTH, ZERO, CYCLIC_X_PD)
      cPG(2) = merge(DELP_Y/YLENGTH, ZERO, CYCLIC_Y_PD)
      cPG(3) = merge(DELP_Z/ZLENGTH, ZERO, CYCLIC_Z_PD)

      DO K = lo(3),hi(3)
         DO J = lo(2),hi(2)
            DO I = lo(1),hi(1)
               gradPg(I,J,K,:) = cPG - gradPg(I,J,K,:)
            ENDDO
         ENDDO
      ENDDO

      IF(DES_EXPLICITLY_COUPLED) THEN

         odx = 1.0d0/dx
         ody = 1.0d0/dy
         odz = 1.0d0/dz

            ! Calculate the gas phase forces acting on each particle.
         DO NP=1,MAX_PIP

            if(nonexistent==particle_state(np) .or.        &
               entering_particle==particle_state(np) .or.  &
               entering_ghost==particle_state(np) .or.     &
               exiting_particle==particle_state(np)  .or.  &
               exiting_ghost==particle_state(np)) cycle

            ! Fluid cell containing the particle
            i = floor(des_pos_new(np,1)*odx) - 1
            j = floor(des_pos_new(np,2)*ody) - 1
            k = floor(des_pos_new(np,3)*odz) - 1

            ! Include gas pressure and gas-solids drag
            drag_fc(NP,:) = drag_fc(NP,:) + gradPg(i,j,k,:)*PVOL(NP)

         ENDDO
      ENDIF

      end subroutine calc_pg_grad

end module calc_pg_grad_module
