!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: calc_outflow                                            C
!  Purpose: Calculate mass and volumetric out flow rates from an       C
!  outflow boundary                                                    C
!                                                                      C
!  Author: M. Syamlal                                 Date: 28-OCT-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_OUTFLOW(L)

! Modules
!--------------------------------------------------------------------//
      use bc, only: bc_plane
      use bc, only: bc_k_b, bc_k_t
      use bc, only: bc_j_s, bc_j_n
      use bc, only: bc_i_w, bc_i_e
      use bc, only: bc_out_n
      use bc, only: bc_mout_g, bc_vout_g
      use geometry, only: dx, dy, dz
      use fldvar, only: u_g, v_g, w_g
      use fldvar, only: rop_g, ep_g
      use functions, only: fluid_at
      use functions, only: is_on_mype_plus2layers
      use functions, only: funijk, iplus, iminus, jplus, jminus, kplus, kminus
      use compar, only: dead_cell_at
      IMPLICIT NONE

! Dummy arguments
!--------------------------------------------------------------------//
! Boundary condition number
      INTEGER, INTENT(IN) :: L

! Local variables
!--------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IJK
! ijk index of fluid cell adjacent to boundary cell
      INTEGER :: IJK2
!--------------------------------------------------------------------//

      BC_OUT_N(L) = BC_OUT_N(L) + 1
      DO K = BC_K_B(L), BC_K_T(L)
         DO J = BC_J_S(L), BC_J_N(L)
            DO I = BC_I_W(L), BC_I_E(L)
! Check if current i,j,k resides on this PE
               IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I,J,K)
               SELECT CASE (TRIM(BC_PLANE(L)))
               CASE ('W')
                  IJK2 = FUNIJK(iminus(i,j,k),j,k)
                  BC_MOUT_G(L) = BC_MOUT_G(L) + DY*DZ*&
                     U_G(iminus(i,j,k),j,k)*ROP_G(iminus(i,j,k),j,k)
                  BC_VOUT_G(L) = BC_VOUT_G(L) + DY*DZ*&
                     U_G(iminus(i,j,k),j,k)*EP_G(IJK2)
               CASE ('E')
                  IJK2 = FUNIJK(iplus(i,j,k),j,k)
                  BC_MOUT_G(L) = BC_MOUT_G(L) + DY*DZ*&
                     U_G(I,J,K)*ROP_G(iplus(i,j,k),j,k)
                  BC_VOUT_G(L) = BC_VOUT_G(L) + DY*DZ*&
                     U_G(I,J,K)*EP_G(IJK2)
               CASE ('S')
                  IJK2 = FUNIJK(i,jminus(i,j,k),k)
                  BC_MOUT_G(L) = BC_MOUT_G(L) + DX*DZ*&
                     V_G(i,jminus(i,j,k),k)*ROP_G(i,jminus(i,j,k),k)
                  BC_VOUT_G(L) = BC_VOUT_G(L) + DX*DZ*&
                     V_G(i,jminus(i,j,k),k)*EP_G(IJK2)
               CASE ('N')
                  IJK2 = FUNIJK(i,jplus(i,j,k),k)
                  BC_MOUT_G(L) = BC_MOUT_G(L) + DX*DZ*&
                     V_G(I,J,K)*ROP_G(i,jplus(i,j,k),k)
                  BC_VOUT_G(L) = BC_VOUT_G(L) + DX*DZ*&
                     V_G(I,J,K)*EP_G(IJK2)
               CASE ('B')
                  IJK2 = FUNIJK(i,j,kminus(i,j,k))
                  BC_MOUT_G(L) = BC_MOUT_G(L) + DX*DY*&
                     W_G(i,j,kminus(i,j,k))*ROP_G(i,j,kminus(i,j,k))
                  BC_VOUT_G(L)=BC_VOUT_G(L)+DX*DY*&
                     W_G(i,j,kminus(i,j,k))*EP_G(IJK2)
               CASE ('T')
                  IJK2 = FUNIJK(i,j,kplus(i,j,k))
                  BC_MOUT_G(L)=BC_MOUT_G(L)+DX*DY*&
                     W_G(I,J,K)*ROP_G(i,j,kplus(i,j,k))
                  BC_VOUT_G(L) = BC_VOUT_G(L) + DX*DY*&
                     W_G(I,J,K)*EP_G(IJK2)
               END SELECT

            ENDDO   ! end do loop (i=bc_i_w(l), bc_i_e(l))
         ENDDO   ! end do loop (j=bc_j_s(l), bc_j_n(l))
      ENDDO   ! end do loop (k=bc_k_b(l), bc_k_t(l))

      RETURN
      END SUBROUTINE CALC_OUTFLOW
