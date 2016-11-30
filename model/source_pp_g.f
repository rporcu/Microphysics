module source_pp_module

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

subroutine source_pp_g(A_M, B_M, B_MMAX, u_g, v_g, w_g, p_g, ep_g,&
                       rop_g, rop_go, ro_g, d_e, d_n, d_t)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc, ONLY: SMALL_NUMBER, ONE, ZERO, UNDEFINED, IJK_P_G
      USE compar, ONLY: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE eos, ONLY: DROODP_G
      USE geometry, ONLY: VOL
      use matrix, ONLY: E, W, N, S, T, B
      USE fld_const, ONLY: RO_G0
      USE run, ONLY: ODT, UNDEFINED_I
      USE ur_facs, ONLY: UR_FAC
      USE xsi_array, ONLY: LOCK_XSI_ARRAY, UNLOCK_XSI_ARRAY

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
! maximum term in b_m expression
      DOUBLE PRECISION, INTENT(INOUT) :: B_mmax&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      DOUBLE PRECISION, INTENT(INOUT) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: d_e&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: d_n&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: d_t&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      integer :: i,j,k
! under relaxation factor for pressure
      DOUBLE PRECISION fac
! terms of bm expression
      DOUBLE PRECISION bma, bme, bmw, bmn, bms, bmt, bmb
! error message
      CHARACTER(LEN=80) :: LINE(1)
! temporary use of global arrays:
! xsi_array: convection weighting factors
!      DOUBLE PRECISION XSI_e(istart3:iend3, jstart3:jend3, kstart3:kend3)
!      DOUBLE PRECISION XSI_n(istart3:iend3, jstart3:jend3, kstart3:kend3)
!      DOUBLE PRECISION XSI_t(istart3:iend3, jstart3:jend3, kstart3:kend3)

!-----------------------------------------------
      call lock_xsi_array

! Calculate convection-diffusion fluxes through each of the faces

        do k = kstart3, kend3
           do j = jstart3, jend3
             do i = istart3, iend3

            IF (fluid_at(i,j,k)) THEN

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
               ELSEIF (RO_G0 .NE. UNDEFINED) THEN !This is an error only in incompressible flow
                  WRITE (LINE, '(A,I3,1x,I3,1x,I3,A,I1,A,G12.5)') 'Error: At IJK = ', i,j,k, &
                     ' M = ', 0, ' A = 0 and b = ', B_M(I,J,K)
                  CALL WRITE_ERROR ('SOURCE_Pp_g', LINE, 1)
               ENDIF
            ENDIF

         ELSE   ! if/else branch .not.fluid_at(i,j,k)
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
         ENDIF   ! end if/else branch fluid_at(i,j,k)
      ENDDO
      ENDDO
      ENDDO

! make correction for compressible flows
      IF (RO_G0 == UNDEFINED) THEN
         fac = UR_FAC(1)  !since p_g = p_g* + ur_fac * pp_g

         do k = kstart3, kend3
            do j = jstart3, jend3
              do i = istart3, iend3

               if (fluid_at(i,j,k)) THEN
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

         IF (fluid_at(i,j,k)) THEN

! Cutting the neighbor link between fluid cell and adjacent p_flow_at cell
            if(p_flow_at(iminus(i,j,k),j,k)) A_m(I,J,K,W) = ZERO
            if(p_flow_at(iplus(i,j,k),j,k))  A_m(I,J,K,E) = ZERO
            if(p_flow_at(i,jminus(i,j,k),k)) A_m(I,J,K,S) = ZERO
            if(p_flow_at(i,jplus(i,j,k),k))  A_m(I,J,K,N) = ZERO
            if(p_flow_at(i,j,kminus(i,j,k))) A_m(I,J,K,B) = ZERO
            if(p_flow_at(i,j,kplus(i,j,k)))  A_m(I,J,K,T) = ZERO
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

      call unlock_xsi_array

      RETURN

CONTAINS

      INCLUDE 'functions.inc'

end subroutine source_pp_g

end module source_pp_module
