module source_pp_module

  use bl_fort_module, only : c_real
  use iso_c_binding , only: c_int

  contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_Pp_g                                             C
!  Purpose: Determine source terms for Pressure correction equation.   C
!                                                                      C
!  Notes: The off-diagonal coefficients are positive. The center       C
!         coefficient and the source vector are negative. See          C
!         conv_Pp_g                                                    C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

subroutine source_pp_g(A_M, B_M, B_MMAX, dt, u_g, v_g, w_g, p_g, ep_g,&
                       rop_g, rop_go, ro_g, d_e, d_n, d_t, flag)

      USE bc, ONLY: SMALL_NUMBER, ONE, ZERO, IJK_P_G
      USE compar, ONLY: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE eos, ONLY: DROODP_G
      USE fld_const, ONLY: RO_G0
      USE geometry, ONLY: VOL
      USE matrix, ONLY: E, W, N, S, T, B
      USE param1, ONLY: IS_DEFINED, IS_UNDEFINED
      USE run, ONLY: UNDEFINED_I
      USE ur_facs, ONLY: UR_FAC
      USE write_error_module, only: write_error

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      real(c_real), INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      real(c_real), INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
! maximum term in b_m expression
      real(c_real), INTENT(INOUT) :: B_mmax&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      double precision, intent(in   ) :: dt

      real(c_real), INTENT(IN   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), INTENT(IN   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), INTENT(IN   ) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), INTENT(IN   ) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), INTENT(IN   ) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), INTENT(IN   ) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), INTENT(IN   ) :: rop_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), INTENT(IN   ) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), INTENT(IN   ) :: d_e&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), INTENT(IN   ) :: d_n&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), INTENT(IN   ) :: d_t&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      INTEGER, INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      integer :: i,j,k
! under relaxation factor for pressure
      real(c_real) fac
! terms of bm expression
      real(c_real) bma, bme, bmw, bmn, bms, bmt, bmb
! error message
      CHARACTER(LEN=80) :: LINE(1)
      double precision :: oDT
!-----------------------------------------------

      odt = 1.0d0/dt

! Calculate convection-diffusion fluxes through each of the faces

        do k = kstart3, kend3
           do j = jstart3, jend3
             do i = istart3, iend3

            IF (1.eq.flag(i,j,k,1)) THEN

            bma = (ROP_G(I,J,K)-ROP_GO(I,J,K))*VOL*ODT
            bme = A_M(I,J,K,E)*U_G(I,J,K)
            bmw = A_M(I,J,K,W)*U_G(iminus(i,j,k),j,k)
            bmn = A_M(I,J,K,N)*V_G(I,J,K)
            bms = A_M(I,J,K,S)*V_G(i,jminus(i,j,k),k)
            bmt = A_M(I,J,K,T)*W_G(I,J,K)
            bmb = A_M(I,J,K,B)*W_G(i,j,kminus(i,j,k))
            B_M(I,J,K) = -((-(bma + bme - bmw + bmn - bms + bmt - bmb )) )
            B_MMAX(I,J,K) = max(abs(bma), abs(bme), abs(bmw), abs(bmn), &
               abs(bms), abs(bmt), abs(bmb))

            A_M(I,J,K,E) = A_M(I,J,K,E)*D_E(i,j,k)
            A_M(I,J,K,W) = A_M(I,J,K,W)*D_E(iminus(i,j,k),j,k)
            A_M(I,J,K,N) = A_M(I,J,K,N)*D_N(i,j,k)
            A_M(I,J,K,S) = A_M(I,J,K,S)*D_N(i,jminus(i,j,k),k)
            A_M(I,J,K,T) = A_M(I,J,K,T)*D_T(i,j,k)
            A_M(I,J,K,B) = A_M(I,J,K,B)*D_T(i,j,kminus(i,j,k))

            A_M(I,J,K,0) = -(A_M(I,J,K,E)+A_M(I,J,K,W)+A_M(I,J,K,N)+A_M(I,J,K,S) + &
               A_M(I,J,K,T)+A_M(I,J,K,B))

            IF (ABS(A_M(I,J,K,0)) < SMALL_NUMBER) THEN
               IF (ABS(B_M(I,J,K)) < SMALL_NUMBER) THEN
                  A_M(I,J,K,0) = -ONE
                  B_M(I,J,K) = ZERO
               ELSEIF (IS_DEFINED(RO_G0)) THEN !This is an error only in incompressible flow
                  WRITE (LINE, '(A,I3,1x,I3,1x,I3,A,I1,A,G12.5)') 'Error: At IJK = ', i,j,k, &
                     ' M = ', 0, ' A = 0 and b = ', B_M(I,J,K)
                  CALL WRITE_ERROR ('SOURCE_Pp_g', LINE, 1)
               ENDIF
            ENDIF

         ELSE   ! if/else branch .not.1.eq.flag(i,j,k)
! set the value (correction) in all wall and flow boundary cells to zero
! note the matrix coefficients and source vector should already be zero
! from the initialization of A and B but the following ensures the zero
! value.
            A_M(I,J,K,E) = ZERO
            A_M(I,J,K,W) = ZERO
            A_M(I,J,K,N) = ZERO
            A_M(I,J,K,S) = ZERO
            A_M(I,J,K,T) = ZERO
            A_M(I,J,K,B) = ZERO
            A_M(I,J,K,0) = -ONE
            B_M(I,J,K) = ZERO
         ENDIF   ! end if/else branch 1.eq.flag(i,j,k)
      ENDDO
      ENDDO
      ENDDO

! make correction for compressible flows
      IF (IS_UNDEFINED(RO_G0)) THEN
         fac = UR_FAC(1)  !since p_g = p_g* + ur_fac * pp_g

         do k = kstart3, kend3
            do j = jstart3, jend3
              do i = istart3, iend3

               if (1.eq.flag(i,j,k,1)) THEN
                  A_M(I,J,K,0) = A_M(I,J,K,0) - &
                     fac*DROODP_G(RO_G(I,J,K),P_G(I,J,K))*&
                     EP_G(I,J,K)*VOL*ODT
               end if

             end do
           end do
         end do
      ENDIF   ! end if (ro_g0 == undefined); i.e., compressible flow


! Remove the asymmetry in matrix caused by the pressure outlet or inlet
! boundaries.  Because the P' at such boundaries is zero we may set the
! coefficient in the neighboring fluid cell to zero without affecting
! the linear equation set.
      do k = kstart3, kend3
         do j = jstart3, jend3
           do i = istart3, iend3

         IF (1.eq.flag(i,j,k,1)) THEN

! Cut the neighbor link between fluid cell and adjacent pressure flow cell
            if(flag(iminus(i,j,k),j,k,1) == 10 .OR. &
               flag(iminus(i,j,k),j,k,1) == 11) A_m(I,J,K,W) = ZERO
            if(flag(iplus(i,j,k),j,k ,1) == 10 .OR. &
               flag(iplus(i,j,k),j,k ,1) == 11) A_m(I,J,K,E) = ZERO
            if(flag(i,jminus(i,j,k),k,1) == 10 .OR. &
               flag(i,jminus(i,j,k),k,1) == 11) A_m(I,J,K,S) = ZERO
            if(flag(i,jplus(i,j,k),k ,1) == 10 .OR. &
               flag(i,jplus(i,j,k),k ,1) == 11) A_m(I,J,K,N) = ZERO
            if(flag(i,j,kminus(i,j,k),1) == 10 .OR. &
               flag(i,j,kminus(i,j,k),1) == 11) A_m(I,J,K,B) = ZERO
            if(flag(i,j,kplus(i,j,k) ,1) == 10 .OR. &
               flag(i,j,kplus(i,j,k) ,1) == 11) A_m(I,J,K,T) = ZERO

         ENDIF
          end do
        end do
      end do

! Specify P' to zero for incompressible flows. Check set_bc0
! for details on selection of IJK_P_g.
      IF (IJK_P_G(1) /= UNDEFINED_I) THEN
         B_M(ijk_p_g(1),ijk_p_g(2),ijk_p_g(3)) = ZERO
         A_M(ijk_p_g(1),ijk_p_g(2),ijk_p_g(3),:) = ZERO
         A_M(ijk_p_g(1),ijk_p_g(2),ijk_p_g(3),0) = -ONE
      ENDIF

      RETURN

CONTAINS

      INCLUDE 'functions.inc'

end subroutine source_pp_g

end module source_pp_module
