MODULE CALC_OUTFLOW_MODULE

   use amrex_fort_module, only : c_real => amrex_real
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

      SUBROUTINE CALC_OUTFLOW(L,slo,shi,ulo,uhi,vlo,vhi,wlo,whi,u_g,v_g,w_g,rop_g,ep_g,dx,dy,dz)

! Modules
!--------------------------------------------------------------------//
      use bc, only: bc_plane
      use bc, only: bc_k_b, bc_k_t
      use bc, only: bc_j_s, bc_j_n
      use bc, only: bc_i_w, bc_i_e
      use bc, only: bc_out_n
      use bc, only: bc_mout_g, bc_vout_g

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)


      ! Boundary condition number
      integer, intent(in) :: L

      real(c_real), intent(in) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in) :: dx, dy, dz

! Local variables
!--------------------------------------------------------------------//
! Indices
      integer :: I, J, K
!--------------------------------------------------------------------//

      BC_OUT_N(L) = BC_OUT_N(L) + 1
      DO K = BC_K_B(L), BC_K_T(L)
         DO J = BC_J_S(L), BC_J_N(L)
            DO I = BC_I_W(L), BC_I_E(L)

               ! Check if current i,j,k resides on this PE
               SELECT CASE (TRIM(BC_PLANE(L)))
               CASE ('W')
                  BC_MOUT_G(L) = BC_MOUT_G(L) + DY*DZ*&
                     U_G(i-1,j,k)*ROP_G(i-1,j,k)
                  BC_VOUT_G(L) = BC_VOUT_G(L) + DY*DZ*&
                     U_G(i-1,j,k)*EP_G(i-1,j,k)
               CASE ('E')
                  BC_MOUT_G(L) = BC_MOUT_G(L) + DY*DZ*&
                     U_G(I,J,K)*ROP_G(i+1,j,k)
                  BC_VOUT_G(L) = BC_VOUT_G(L) + DY*DZ*&
                     U_G(I,J,K)*EP_G(i+1,j,k)
               CASE ('S')
                  BC_MOUT_G(L) = BC_MOUT_G(L) + DX*DZ*&
                     V_G(i,j-1,k)*ROP_G(i,j-1,k)
                  BC_VOUT_G(L) = BC_VOUT_G(L) + DX*DZ*&
                     V_G(i,j-1,k)*EP_G(i,j-1,k)
               CASE ('N')
                  BC_MOUT_G(L) = BC_MOUT_G(L) + DX*DZ*&
                     V_G(I,J,K)*ROP_G(i,j+1,k)
                  BC_VOUT_G(L) = BC_VOUT_G(L) + DX*DZ*&
                     V_G(I,J,K)*EP_G(i,j+1,k)
               CASE ('B')
                  BC_MOUT_G(L) = BC_MOUT_G(L) + DX*DY*&
                     W_G(i,j,k-1)*ROP_G(i,j,k-1)
                  BC_VOUT_G(L)=BC_VOUT_G(L)+DX*DY*&
                     W_G(i,j,k-1)*EP_G(i,j,k-1)
               CASE ('T')
                  BC_MOUT_G(L)=BC_MOUT_G(L)+DX*DY*&
                     W_G(I,J,K)*ROP_G(i,j,k+1)
                  BC_VOUT_G(L) = BC_VOUT_G(L) + DX*DY*&
                     W_G(I,J,K)*EP_G(i,j,k+1)
               END SELECT

            ENDDO   ! end do loop (i=bc_i_w(l), bc_i_e(l))
         ENDDO   ! end do loop (j=bc_j_s(l), bc_j_n(l))
      ENDDO   ! end do loop (k=bc_k_b(l), bc_k_t(l))

      RETURN
      END SUBROUTINE CALC_OUTFLOW
END MODULE CALC_OUTFLOW_MODULE
