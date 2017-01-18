MODULE CONV_ROP_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_ROP                                                C
!  Purpose: Calculate the face value of density used for calculating   C
!  convection fluxes. Master routine.                                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CONV_ROP(slo, shi, lo, hi, &
         u_g, ulo, uhi, v_g, vlo, vhi, w_g, wlo, whi, &
         rop_g, rop_ge, rop_gn, rop_gt, &
         flag, dt, dx, dy, dz) bind(C, name="conv_rop")

! Modules
!---------------------------------------------------------------------//
      USE run, only: discretize

      IMPLICIT NONE

      integer(c_int), intent(in   ) :: slo(3), shi(3), lo(3), hi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3), vlo(3), vhi(3), wlo(3), whi(3)

      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(  out) :: rop_ge&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(  out) :: rop_gn&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(  out) :: rop_gt&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(in   ) :: dt, dx, dy, dz
!---------------------------------------------------------------------//
      stop

      IF (DISCRETIZE(1) == 0) THEN       ! 0 & 1 => first order upwinding
         CALL CONV_ROP0 (slo, shi, lo, hi, ROP_g, &
                         u_g, ulo, uhi, v_g, vlo, vhi, w_g, wlo, whi, &
                         ROP_gE, ROP_gN, ROP_gT, flag)
      ELSE
         CALL CONV_ROP1 (DISCRETIZE(1), &
                         slo, shi, lo, hi, &
                         flag, rop_g, &
                         u_g, ulo, uhi, v_g, vlo, vhi, w_g, wlo, whi, &
                         rop_ge, rop_gn, rop_gt, dt, dx, dy, dz)
      ENDIF

      RETURN
   CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_ROP                                                C
!  Purpose: Calculate the face value of density used for calculating   C
!  convection fluxes. FOU routine.                                     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CONV_ROP0(slo, shi, lo, hi, ROP, &
                           u, ulo, uhi, v, vlo, vhi, w, wlo, whi, &
                           ROP_E, ROP_N, ROP_T, flag)

! Modules
!---------------------------------------------------------------------//
      USE functions, only: ieast, jnorth, ktop
      USE functions, only: iwest, jsouth, kbot
      USE functions, only: iminus, jminus, kminus
      USE param1, only: zero

      IMPLICIT NONE

      integer(c_int), intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3), vlo(3), vhi(3), wlo(3), whi(3)

! macroscopic density (rho_prime)
      real(c_real), intent(in   ) :: rop&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! Velocity components
      real(c_real), intent(in   ) :: u&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

! Face value of density (for calculating convective fluxes)
      real(c_real), intent(  out) :: rop_e&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(  out) :: rop_n&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(  out) :: rop_t&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      integer, intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K
!---------------------------------------------------------------------//

      DO K = ulo(3),uhi(3)
        DO J = ulo(2),uhi(2)
          DO I = ulo(1),uhi(1)

            IF (1.eq.flag(i,j,k,1)) THEN

            ! East face (i+1/2, j, k)
            IF (U(i,j,k) >= ZERO) THEN
               ROP_E(i,j,k) = ROP(i,j,k)
            ELSE
               ROP_E(i,j,k) = ROP(ieast(i,j,k),j,k)
            ENDIF

            ! West face (i-1/2, j, k)
            IF (.NOT.1.eq.flag(iminus(i,j,k),j,k,1)) THEN
               IF (U(iminus(i,j,k),j,k) >= ZERO) THEN
                  ROP_E(iminus(i,j,k),j,k) = ROP(iwest(i,j,k),j,k)
               ELSE
                  ROP_E(iminus(i,j,k),j,k) = ROP(i,j,k)
               ENDIF
            ENDIF

         ENDIF
      ENDDO
      ENDDO
      ENDDO

      DO K = vlo(3),vhi(3)
        DO J = vlo(2),vhi(2)
          DO I = vlo(1),vhi(1)

            IF (1.eq.flag(i,j,k,1)) THEN

            ! North face (i, j+1/2, k)
            IF (V(i,j,k) >= ZERO) THEN
               ROP_N(i,j,k) = ROP(i,j,k)
            ELSE
               ROP_N(i,j,k) = ROP(i,jnorth(i,j,k),k)
            ENDIF

            ! South face (i, j-1/2, k)
            IF (.NOT.1.eq.flag(i,jminus(i,j,k),k,1)) THEN
               IF (V(i,jminus(i,j,k),k) >= ZERO) THEN
                 ROP_N(i,jminus(i,j,k),k) = ROP(i,jsouth(i,j,k),k)
               ELSE
                 ROP_N(i,jminus(i,j,k),k) = ROP(i,j,k)
               ENDIF
            ENDIF

         ENDIF
      ENDDO
      ENDDO
      ENDDO

      DO K = wlo(3),whi(3)
        DO J = wlo(2),whi(2)
          DO I = wlo(1),whi(1)

           IF (1.eq.flag(i,j,k,1)) THEN

            ! Top face (i, j, k+1/2)
            IF (W(i,j,k) >= ZERO) THEN
               ROP_T(i,j,k) = ROP(i,j,k)
            ELSE
               ROP_T(i,j,k) = ROP(i,j,ktop(i,j,k))
            ENDIF

            ! Bottom face (i, j, k-1/2)
            IF (.NOT.1.eq.flag(i,j,kminus(i,j,k),1)) THEN
               IF (W(i,j,kminus(i,j,k)) >= ZERO) THEN
                  ROP_T(i,j,kminus(i,j,k)) = ROP(i,j,kbot(i,j,k))
               ELSE
                  ROP_T(i,j,kminus(i,j,k)) = ROP(i,j,k)
               ENDIF
            ENDIF

         ENDIF
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CONV_ROP0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_ROP1                                               C
!  Purpose: Calculate the face value of density used for calculating   C
!  convection fluxes. HR routine.  Here interpolate the face value of  C
!  density.                                                            C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CONV_ROP1(DISC, &
                           slo, shi, lo, hi, flag, rop, &
                           u, ulo, uhi, v, vlo, vhi, w, wlo, whi, &
                           rop_e, rop_n, rop_t, &
                           dt, dx, dy, dz)

! Modules
!---------------------------------------------------------------------//
      USE functions, only: ieast, jnorth, ktop
      USE functions, only: iwest, jsouth, kbot
      USE functions, only: iminus, jminus, kminus
      USE param1, only: one
      USE xsi, only: calc_xsi
      IMPLICIT NONE

      integer(c_int), intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3), vlo(3), vhi(3), wlo(3), whi(3)

! Discretization scheme
      INTEGER, INTENT(IN) :: DISC

! macroscopic density (rho_prime)
      real(c_real), INTENT(in) :: rop&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! Velocity components
      real(c_real), INTENT(IN   ) :: u&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), INTENT(IN   ) :: v&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), INTENT(IN   ) :: w&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

! Face value of density (for calculating convective fluxes)
      real(c_real), intent(  out) :: rop_e&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(  out) :: rop_n&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(  out) :: rop_t&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      integer, intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(in) :: dt, dx, dy, dz
!
! Local variables
!---------------------------------------------------------------------//
      INTEGER :: I,J,K
      Integer :: incr
      real(c_real), allocatable :: xsi_e(:,:,:), xsi_n(:,:,:), xsi_t(:,:,:)

!---------------------------------------------------------------------//

      allocate( xsi_e(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) )
      allocate( xsi_n(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) )
      allocate( xsi_t(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) )


! Calculate factors
      incr=0
      CALL CALC_XSI (DISC, slo, shi, lo, hi, &
                     ROP, U, V, W, XSI_E, XSI_N, XSI_T, incr, dt, dx, dy, dz)

      DO K = slo(3),shi(3)
        DO J = slo(2),shi(2)
          DO I = slo(1),shi(1)

         IF (1.eq.flag(i,j,k,1)) THEN

! East face (i+1/2, j, k)
            ROP_E(i,j,k) = ((ONE-XSI_E(i,j,k))*ROP(i,j,k) + &
               XSI_E(i,j,k) *ROP(ieast(i,j,k),j,k) )
! West face (i-1/2, j, k)
            IF (.NOT.1.eq.flag(iminus(i,j,k),j,k,1)) THEN
               ROP_E(iminus(i,j,k),j,k) = &
                  ((ONE - XSI_E(iminus(i,j,k),j,k))*ROP(iwest(i,j,k),j,k) + &
                  XSI_E(iminus(i,j,k),j,k) *ROP(i,j,k) )
            ENDIF


! North face (i, j+1/2, k)
            ROP_N(i,j,k) = ((ONE-XSI_N(i,j,k))*ROP(i,j,k)+&
               XSI_N(i,j,k) *ROP(i,jnorth(i,j,k),k))
! South face (i, j-1/2, k)
            IF (.NOT.1.eq.flag(i,jminus(i,j,k),k,1)) THEN
               ROP_N(i,jminus(i,j,k),k) = &
                  ((ONE - XSI_N(i,jminus(i,j,k),k))*ROP(i,jsouth(i,j,k),k) + &
                  XSI_N(i,jminus(i,j,k),k) *ROP(i,j,k) )
            ENDIF

! Top face (i, j, k+1/2)
            ROP_T(i,j,k) = ((ONE - XSI_T(i,j,k))*ROP(i,j,k) + &
               XSI_T(i,j,k) *ROP(i,j,ktop(i,j,k)) )
! Bottom face (i, j, k-1/2)
            IF (.NOT.1.eq.flag(i,j,kminus(i,j,k),1)) THEN
               ROP_T(i,j,kminus(i,j,k)) = &
                  ((ONE - XSI_T(i,j,kminus(i,j,k)))*ROP(i,j,kbot(i,j,k)) + &
                  XSI_T(i,j,kminus(i,j,k)) *ROP(i,j,k) )
            ENDIF

         ENDIF
      ENDDO
      ENDDO
      ENDDO

      deallocate( xsi_e, xsi_n, xsi_t)

      RETURN
      END SUBROUTINE CONV_ROP1

   END SUBROUTINE CONV_ROP

END MODULE CONV_ROP_MODULE
