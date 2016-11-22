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

SUBROUTINE SOURCE_PP_G(A_M, B_M, B_MMAX)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc, ONLY: SMALL_NUMBER, ONE, ZERO, UNDEFINED, IJK_P_G, DIMENSION_3
      USE compar, ONLY: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE eos, ONLY: DROODP_G
      USE fldvar, ONLY: U_G, V_G, W_G,ROP_G, ROP_GO, RO_G, P_G, EP_G
      USE geometry, ONLY: VOL
      use matrix, ONLY: E, W, N, S, T, B
      USE fldvar, ONLY: D_E, D_N, D_T
      USE fldvar, ONLY: RO_G0
      USE run, ONLY: ODT, UNDEFINED_I
      USE ur_facs, ONLY: UR_FAC
      USE xsi_array, ONLY: LOCK_XSI_ARRAY, UNLOCK_XSI_ARRAY

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3)
! maximum term in b_m expression
      DOUBLE PRECISION, INTENT(INOUT) :: B_mmax(DIMENSION_3)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      integer :: i,j,k
      INTEGER :: IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
! under relaxation factor for pressure
      DOUBLE PRECISION fac
! terms of bm expression
      DOUBLE PRECISION bma, bme, bmw, bmn, bms, bmt, bmb
! error message
      CHARACTER(LEN=80) :: LINE(1)
! temporary use of global arrays:
! xsi_array: convection weighting factors
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),&
!                       XSI_t(DIMENSION_3)
!-----------------------------------------------
      call lock_xsi_array

! Calculate convection-diffusion fluxes through each of the faces

        do k = kstart3, kend3
           do j = jstart3, jend3
             do i = istart3, iend3

         ijk = funijk(i,j,k)

         IF (fluid_cell(i,j,k)) THEN
            IMJK = FUNIJK(iminus(i,j,k),j,k)
            IJMK = FUNIJK(i,jminus(i,j,k),k)
            IJKM = FUNIJK(i,j,kminus(i,j,k))

            bma = (ROP_G(IJK)-ROP_GO(I,J,K))*VOL*ODT
            bme = A_M(IJK,E)*U_G(IJK)
            bmw = A_M(IJK,W)*U_G(IMJK)
            bmn = A_M(IJK,N)*V_G(IJK)
            bms = A_M(IJK,S)*V_G(IJMK)
            bmt = A_M(IJK,T)*W_G(IJK)
            bmb = A_M(IJK,B)*W_G(IJKM)
            B_M(IJK) = -((-(bma + bme - bmw + bmn - bms + bmt - bmb )) )
            B_MMAX(IJK) = max(abs(bma), abs(bme), abs(bmw), abs(bmn), abs(bms), abs(bmt), abs(bmb))

            A_M(IJK,E) = A_M(IJK,E)*D_E(i,j,k)
            A_M(IJK,W) = A_M(IJK,W)*D_E(iminus(i,j,k),j,k)
            A_M(IJK,N) = A_M(IJK,N)*D_N(i,j,k)
            A_M(IJK,S) = A_M(IJK,S)*D_N(i,jminus(i,j,k),k)
            A_M(IJK,T) = A_M(IJK,T)*D_T(i,j,k)
            A_M(IJK,B) = A_M(IJK,B)*D_T(i,j,kminus(i,j,k))

            A_M(IJK,0) = -(A_M(IJK,E)+A_M(IJK,W)+A_M(IJK,N)+A_M(IJK,S) + &
               A_M(IJK,T)+A_M(IJK,B))

            IF (ABS(A_M(IJK,0)) < SMALL_NUMBER) THEN
               IF (ABS(B_M(IJK)) < SMALL_NUMBER) THEN
                  A_M(IJK,0) = -ONE
                  B_M(IJK) = ZERO
               ELSEIF (RO_G0 .NE. UNDEFINED) THEN !This is an error only in incompressible flow
                  WRITE (LINE, '(A,I6,A,I1,A,G12.5)') 'Error: At IJK = ', IJK, &
                     ' M = ', 0, ' A = 0 and b = ', B_M(IJK)
                  CALL WRITE_ERROR ('SOURCE_Pp_g', LINE, 1)
               ENDIF
            ENDIF

         ELSE   ! if/else branch .not.fluid_cell(i,j,k)
! set the value (correction) in all wall and flow boundary cells to zero
! note the matrix coefficients and source vector should already be zero
! from the initialization of A and B but the following ensures the zero
! value.
            A_M(IJK,E) = ZERO
            A_M(IJK,W) = ZERO
            A_M(IJK,N) = ZERO
            A_M(IJK,S) = ZERO
            A_M(IJK,T) = ZERO
            A_M(IJK,B) = ZERO
            A_M(IJK,0) = -ONE
            B_M(IJK) = ZERO
         ENDIF   ! end if/else branch fluid_cell(i,j,k)
      ENDDO
      ENDDO
      ENDDO

! make correction for compressible flows
      IF (RO_G0 == UNDEFINED) THEN
         fac = UR_FAC(1)  !since p_g = p_g* + ur_fac * pp_g

         do k = kstart3, kend3
            do j = jstart3, jend3
              do i = istart3, iend3

               ijk = funijk(i,j,k)

               if (fluid_cell(i,j,k)) THEN
                  A_M(IJK,0) = A_M(IJK,0) - &
                     fac*DROODP_G(RO_G(IJK),P_G(IJK))*&
                     EP_G(IJK)*VOL*ODT
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

           ijk = funijk(i,j,k)
         IF (fluid_cell(i,j,k)) THEN
            IMJK = FUNIJK(iminus(i,j,k),j,k)
            IPJK = FUNIJK(iplus(i,j,k),j,k)
            IJMK = FUNIJK(i,jminus(i,j,k),k)
            IJPK = FUNIJK(i,jplus(i,j,k),k)
            IJKM = FUNIJK(i,j,kminus(i,j,k))
            IJKP = FUNIJK(i,j,kplus(i,j,k))

! Cutting the neighbor link between fluid cell and adjacent p_flow_at cell
            if(p_flow_at(imjk)) A_m(IJK,W) = ZERO
            if(p_flow_at(ipjk)) A_m(IJK,E) = ZERO
            if(p_flow_at(ijmk)) A_m(IJK,S) = ZERO
            if(p_flow_at(ijpk)) A_m(IJK,N) = ZERO
            if(p_flow_at(ijkm)) A_m(IJK,B) = ZERO
            if(p_flow_at(ijkp)) A_m(IJK,T) = ZERO
         ENDIF
          end do
        end do
      end do

! Specify P' to zero for incompressible flows. Check set_bc0
! for details on selection of IJK_P_g.
      IF (IJK_P_G(1) /= UNDEFINED_I) THEN
         ijk = funijk(ijk_p_g(1),ijk_p_g(2),ijk_p_g(3))
         B_M(IJK) = ZERO
         A_M(IJK,:) = ZERO
         A_M(IJK,0) = -ONE
      ENDIF

      call unlock_xsi_array

      RETURN

CONTAINS

      INCLUDE 'functions.inc'

END SUBROUTINE SOURCE_PP_G
