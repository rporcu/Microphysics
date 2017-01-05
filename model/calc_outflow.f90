MODULE CALC_OUTFLOW_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: calc_outflow                                            C
!  Purpose: Calculate mass and volumetric out flow rates from an       C
!  outflow boundary                                                    C
!                                                                      C
!  Author: M. Syamlal                                 Date: 28-OCT-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_OUTFLOW(L,u_g,v_g,w_g,rop_g,ep_g,dx,dy,dz)

! Modules
!--------------------------------------------------------------------//
      use bc, only: bc_plane
      use bc, only: bc_k_b, bc_k_t
      use bc, only: bc_j_s, bc_j_n
      use bc, only: bc_i_w, bc_i_e
      use bc, only: bc_out_n
      use bc, only: bc_mout_g, bc_vout_g

      use compar  , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      use functions, only: iplus, iminus, jplus, jminus, kplus, kminus

      implicit none

!--------------------------------------------------------------------//
      ! Boundary condition number
      INTEGER, intent(in) :: L

      real(c_real), intent(in) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in) :: dx, dy, dz

! Local variables
!--------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K
!--------------------------------------------------------------------//

      BC_OUT_N(L) = BC_OUT_N(L) + 1
      DO K = BC_K_B(L), BC_K_T(L)
         DO J = BC_J_S(L), BC_J_N(L)
            DO I = BC_I_W(L), BC_I_E(L)
! Check if current i,j,k resides on this PE
               SELECT CASE (TRIM(BC_PLANE(L)))
               CASE ('W')
                  BC_MOUT_G(L) = BC_MOUT_G(L) + DY*DZ*&
                     U_G(iminus(i,j,k),j,k)*ROP_G(iminus(i,j,k),j,k)
                  BC_VOUT_G(L) = BC_VOUT_G(L) + DY*DZ*&
                     U_G(iminus(i,j,k),j,k)*EP_G(iminus(i,j,k),j,k)
               CASE ('E')
                  BC_MOUT_G(L) = BC_MOUT_G(L) + DY*DZ*&
                     U_G(I,J,K)*ROP_G(iplus(i,j,k),j,k)
                  BC_VOUT_G(L) = BC_VOUT_G(L) + DY*DZ*&
                     U_G(I,J,K)*EP_G(iplus(i,j,k),j,k)
               CASE ('S')
                  BC_MOUT_G(L) = BC_MOUT_G(L) + DX*DZ*&
                     V_G(i,jminus(i,j,k),k)*ROP_G(i,jminus(i,j,k),k)
                  BC_VOUT_G(L) = BC_VOUT_G(L) + DX*DZ*&
                     V_G(i,jminus(i,j,k),k)*EP_G(i,jminus(i,j,k),k)
               CASE ('N')
                  BC_MOUT_G(L) = BC_MOUT_G(L) + DX*DZ*&
                     V_G(I,J,K)*ROP_G(i,jplus(i,j,k),k)
                  BC_VOUT_G(L) = BC_VOUT_G(L) + DX*DZ*&
                     V_G(I,J,K)*EP_G(i,jplus(i,j,k),k)
               CASE ('B')
                  BC_MOUT_G(L) = BC_MOUT_G(L) + DX*DY*&
                     W_G(i,j,kminus(i,j,k))*ROP_G(i,j,kminus(i,j,k))
                  BC_VOUT_G(L)=BC_VOUT_G(L)+DX*DY*&
                     W_G(i,j,kminus(i,j,k))*EP_G(i,j,kminus(i,j,k))
               CASE ('T')
                  BC_MOUT_G(L)=BC_MOUT_G(L)+DX*DY*&
                     W_G(I,J,K)*ROP_G(i,j,kplus(i,j,k))
                  BC_VOUT_G(L) = BC_VOUT_G(L) + DX*DY*&
                     W_G(I,J,K)*EP_G(i,j,kplus(i,j,k))
               END SELECT

            ENDDO   ! end do loop (i=bc_i_w(l), bc_i_e(l))
         ENDDO   ! end do loop (j=bc_j_s(l), bc_j_n(l))
      ENDDO   ! end do loop (k=bc_k_b(l), bc_k_t(l))

      RETURN
      END SUBROUTINE CALC_OUTFLOW
END MODULE CALC_OUTFLOW_MODULE
