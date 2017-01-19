module source_u_g_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int
   use param1        , only: zero, half, one, undefined, is_undefined, small_number

  contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_U_g                                              C
!  Purpose: Determine source terms for U_g momentum eq. The terms      C
!  appear in the center coefficient and RHS vector. The center         C
!  coefficient and source vector are negative.  The off-diagonal       C
!  coefficients are positive.                                          C
!  The drag terms are excluded from the source at this stage.          C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 14-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOURCE_U_G(slo, shi, lo, hi, &
         A_m, b_m, dt, p_g, ep_g, ro_g, rop_g, rop_go, &
         u_g, u_go, tau_u_g, flag, dx, dy, dz)

! Modules
!---------------------------------------------------------------------//
      USE constant, only: gravity
      USE bc      , only: delp_x

      USE functions, only: avg
      USE functions, only: iminus,iplus,jminus,jplus,kminus,kplus,ieast,iwest
      USE functions, only: zmax

      USE geometry, only: domlo, domhi, cyclic_x_pd

      use matrix, only: e, w, s, n, t, b

      USE scales, only: p_scale
      USE toleranc, only: dil_ep_s

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
      real(c_real), intent(in   ) :: dt, dx, dy, dz

      ! Septadiagonal matrix A_m
      real(c_real), INTENT(INOUT) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      ! Vector b_m
      real(c_real), INTENT(INOUT) :: b_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: rop_go&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: u_go&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: tau_u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer, intent(in   ) :: flag &
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

! Local Variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: i,j,k
! Pressure at east cell
      real(c_real) :: PgE
! Average volume fraction
      real(c_real) :: EPGA
! Average density
      real(c_real) :: ROPGA, ROGA
! Source terms (Surface)
      real(c_real) :: Sdp
! Source terms (Volumetric)
      real(c_real) :: V0, Vbf
! local stress tensor quantity
      real(c_real) :: ltau_u_g
      real(c_real) :: odt
      real(c_real) :: ayz
      real(c_real) :: vol
!---------------------------------------------------------------------//

      odt = 1.0d0/dt
      ayz = dy*dz
      vol = dx*dy*dz

      DO K = lo(3), hi(3)
        DO J = lo(2), hi(2)
          DO I = lo(1), hi(1)+1

         EPGA = AVG(EP_G(I,J,K),EP_G(ieast(i,j,k),j,k))

         ! Impermeable internal surface
         IF (flag(i,j,k,2) < 1000) THEN
            A_m(I,J,K,E) = ZERO
            A_m(I,J,K,W) = ZERO
            A_m(I,J,K,N) = ZERO
            A_m(I,J,K,S) = ZERO
            A_m(I,J,K,T) = ZERO
            A_m(I,J,K,B) = ZERO
            A_m(I,J,K,0) = -ONE
            b_m(I,J,K) = ZERO

         ! Dilute flow
         ELSEIF (EPGA <= DIL_EP_S) THEN
            A_m(I,J,K,E) = ZERO
            A_m(I,J,K,W) = ZERO
            A_m(I,J,K,N) = ZERO
            A_m(I,J,K,S) = ZERO
            A_m(I,J,K,T) = ZERO
            A_m(I,J,K,B) = ZERO
            A_m(I,J,K,0) = -ONE
            b_m(I,J,K) = ZERO

            ! set velocity equal to that of west or east cell if solids are present
            ! in those cells else set velocity equal to known value
            IF (EP_G(iwest(i,j,k),j,k) > DIL_EP_S) THEN
               A_m(I,J,K,W) = ONE
            ELSE IF (EP_G(ieast(i,j,k),j,k) > DIL_EP_S) THEN
               A_m(I,J,K,E) = ONE
            ELSE
               b_m(I,J,K) = -U_G(I,J,K)
            ENDIF

         ! Normal case
         ELSE

            ! Pressure term
            PGE = P_G(ieast(i,j,k),j,k)
            if ( CYCLIC_X_PD) then
              if ( (i .eq. domlo(1)-1) .or. (i .eq. domhi(1)) ) &
                PGE = P_G(ieast(i,j,k),j,k) - DELP_X
            end if

            SDP = -P_SCALE*EPGA*(PGE - P_G(I,J,K))*AYZ

            ! Volumetric forces
            ROGA  = HALF * (RO_G(I,J,K) + RO_G(ieast(i,j,k),j,k))
            ROPGA = HALF * (ROP_G(I,J,K) + ROP_G(ieast(i,j,k),j,k))

            ! Previous time step
            V0 = HALF * (ROP_GO(I,J,K) + ROP_GO(ieast(i,j,k),j,k))*ODT

            ! Body force
            VBF = ROGA*GRAVITY(1)

            ltau_u_g = tau_u_g(i,j,k)

! Collect the terms
            A_m(I,J,K,0) = -(A_m(I,J,K,E)+A_m(I,J,K,W)+&
               A_m(I,J,K,N)+A_m(I,J,K,S)+A_m(I,J,K,T)+A_m(I,J,K,B)+&
               V0*VOL)

            b_m(I,J,K) = b_m(I,J,K) -(SDP + lTAU_U_G + &
               ( (V0)*U_GO(I,J,K) + VBF)*VOL )

         ENDIF   ! end branching on cell type (ip/dilute/block/else branches)

          ENDDO   ! end do loop over ijk
        ENDDO   ! end do loop over ijk
     ENDDO   ! end do loop over ijk


     ! do k = slo(3),shi(3)
     !     do j = slo(2),shi(2)
     !        do i = slo(1),shi(1)
     !           write(5555,"(3(i3),8(es12.4))")i,j,k,a_m(i,j,k,:),b_m(i,j,k)
     !        enddo
     !     enddo
     !  enddo

      ! modifications for bc
      CALL SOURCE_U_G_BC (slo, shi, lo, hi, A_m, b_m, U_G, flag, dx, dy, dz)

     ! do k = slo(3),shi(3)
     !     do j = slo(2),shi(2)
     !        do i = slo(1),shi(1)
     !           write(6666,"(3(i3),8(es12.4))")i,j,k,a_m(i,j,k,:),b_m(i,j,k)
     !        enddo
     !     enddo
     !  enddo

     !  stop 99999

      RETURN
      END SUBROUTINE SOURCE_U_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_U_g_BC                                           C
!  Purpose: Determine source terms for U_g momentum eq. The terms      C
!     appear in the center coefficient and RHS vector. The center      C
!     coefficient and source vector are negative. The off-diagonal     C
!     coefficients are positive.                                       C
!     The drag terms are excluded from the source at this stage.       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 15-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_U_G_BC(slo,shi,lo,hi,A_m, b_m, U_G, flag, dx, dy, dz)

      use ic, only: NSW_, FSW_
      USE bc, only: bc_hw_g, bc_uw_g
      USE bc, only: bc_i_w, bc_i_e, bc_j_s, bc_j_n, bc_k_b, bc_k_t
      USE bc, only: dimension_bc, bc_type, bc_defined, bc_plane
      USE functions, only: ieast, iwest, jsouth, jnorth, kbot, ktop
      USE functions, only: iminus, iplus, im1
      USE matrix, only: e, w, s, n, t, b
      use geometry, only: domlo, domhi

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! Septadiagonal matrix A_m
      real(c_real), INTENT(INOUT) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      ! Vector b_m
      real(c_real), INTENT(INOUT) :: b_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Velocity u_g
      real(c_real), INTENT(IN   ) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      INTEGER, INTENT(IN   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(in   ) :: dx, dy, dz
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Boundary condition
      INTEGER :: L
! Indices
      INTEGER ::  I,  J, K, IM, I1, I2, J1, J2, K1, K2

      real(c_real) :: ody, odz
!-----------------------------------------------

     ody = 1.d0 / dy
     odz = 1.d0 / dz

     if(.false.) then

        do k=slo(3),shi(3)
           write(6,"(2/'Scalar flag k-plane: ',i2)") k
           do j=shi(2),slo(2),-1
              write(6,"(3x,'j=',i2,3x)",advance='no') j
              do i=slo(1),shi(1)-1
                 write(6,"(i6)",advance='no') flag(i,j,k,1)
              enddo
              write(6,"(i6)",advance='yes') flag(shi(1),j,k,1)
           enddo
        enddo

        do k=slo(3),shi(3)
           write(6,"(2/'Flag EAST flag k-plane: ',i2)") k
           do j=shi(2),slo(2),-1
              write(6,"(3x,'j=',i2,3x)",advance='no') j
              do i=slo(1),shi(1)-1
                 write(6,"(i6)",advance='no') flag(i,j,k,2)
              enddo
              write(6,"(i6)",advance='yes') flag(shi(1),j,k,2)
           enddo
        enddo

        do k=slo(3),shi(3)
           write(6,"(2/'Flag NORTH flag k-plane: ',i2)") k
           do j=shi(2),slo(2),-1
              write(6,"(3x,'j=',i2,3x)",advance='no') j
              do i=slo(1),shi(1)-1
                 write(6,"(i6)",advance='no') flag(i,j,k,3)
              enddo
              write(6,"(i6)",advance='yes') flag(shi(1),j,k,3)
           enddo
        enddo


        do k=slo(3),shi(3)
           write(6,"(2/'Flag TOP flag k-plane: ',i2)") k
           do j=shi(2),slo(2),-1
              write(6,"(3x,'j=',i2,3x)",advance='no') j
              do i=slo(1),shi(1)-1
                 write(6,"(i6)",advance='no') flag(i,j,k,4)
              enddo
              write(6,"(i6)",advance='yes') flag(shi(1),j,k,4)
           enddo
        enddo

        write(6,"(2/'  ')")
        write(6,"('  slo',3(i3),'      shi',3(i3))")   slo,   shi
        write(6,"('   lo',3(i3),'       hi',3(i3))")    lo,    hi
        write(6,"('domlo',3(i3),'    domhi',3(i3))") domlo, domhi

     endif




     do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

! No-slip walls
               if (flag(i,j-1,k,1) == NSW_) then
                  b_m(i,j-1,k) = zero
                  A_m(i,j-1,k,:) = zero
                  A_m(i,j-1,k,0) = -one
                  A_m(i,j-1,k,n) = -one

               else if (flag(i,j+1,k,1) == NSW_) then
                  b_m(i,j+1,k) = zero
                  A_m(i,j+1,k,:) = zero
                  A_m(i,j+1,k,0) = -one
                  A_m(i,j+1,k,s) = -one

               else if (flag(i,j,k-1,1) == NSW_) then
                  b_m(i,j,k-1) = zero
                  A_m(i,j,k-1,:) = zero
                  A_m(i,j,k-1,0) = -one
                  A_m(i,j,k-1,t) = -one

               else if (flag(i,j,k+1,1) == NSW_) then
                  b_m(i,j,k+1) = zero
                  A_m(i,j,k+1,:) = zero
                  A_m(i,j,k+1,0) = -one
                  A_m(i,j,k+1,b) = -one

! Free-slip walls
               else if (flag(i,j-1,k,1) == FSW_) then
                  b_m(i,j-1,k) = zero
                  A_m(i,j-1,k,:) = zero
                  A_m(i,j-1,k,0) = -one
                  A_m(i,j-1,k,n) =  one

               else if (flag(i,j+1,k,1) == FSW_) then
                  b_m(i,j+1,k) = zero
                  A_m(i,j+1,k,:) = zero
                  A_m(i,j+1,k,0) = -one
                  A_m(i,j+1,k,s) =  one

               else if (flag(i,j,k-1,1) == FSW_) then
                  b_m(i,j,k-1) = zero
                  A_m(i,j,k-1,:) = zero
                  A_m(i,j,k-1,0) = -one
                  A_m(i,j,k-1,t) =  one

               else if (flag(i,j,k+1,1) == FSW_) then
                  b_m(i,j,k+1) = zero
                  A_m(i,j,k+1,:) = zero
                  A_m(i,j,k+1,0) = -one
                  A_m(i,j,k+1,b) =  one

               endif
            enddo
         enddo
      enddo


! Setting user specified boundary conditions
      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN

! Setting wall boundary conditions
! ---------------------------------------------------------------->>>
            IF (BC_TYPE(L) == 'PAR_SLIP_WALL') THEN
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
               DO K = K1, K2
                  DO J = J1, J2
                     DO I = I1, I2
                        IF (flag(i,j,k,1)<100) CYCLE  ! skip redefined cells
                        A_m(I,J,K,E) = ZERO
                        A_m(I,J,K,W) = ZERO
                        A_m(I,J,K,N) = ZERO
                        A_m(I,J,K,S) = ZERO
                        A_m(I,J,K,T) = ZERO
                        A_m(I,J,K,B) = ZERO
                        A_m(I,J,K,0) = -ONE
                        b_m(I,J,K) = ZERO
                        if (flag(i,jnorth(i,j,k),k,1) == 1) THEN
                           IF (IS_UNDEFINED(BC_HW_G(L))) THEN
                              A_m(I,J,K,N) = -HALF
                              A_m(I,J,K,0) = -HALF
                              b_m(I,J,K) = -BC_UW_G(L)
                           ELSE
                              A_m(I,J,K,0) = -(HALF*BC_HW_G(L)+ODY)
                              A_m(I,J,K,N) = -(HALF*BC_HW_G(L)-ODY)
                              b_m(I,J,K) = -BC_HW_G(L)*BC_UW_G(L)
                           ENDIF




                        elseif (flag(i,jsouth(i,j,k),k,1)==1) THEN




                           IF (IS_UNDEFINED(BC_HW_G(L))) THEN

                              ! cn7 = A_m(i,jsouth(i,j,k),k,n)
                              ! cs8 = -half
                              ! c08 = -half
                              ! s8  = -bc_uw_g(l)


                              A_m(i,jsouth(i,j,k),k,0) = A_m(i,jsouth(i,j,k),k,0) - A_m(i,jsouth(i,j,k),k,n)
                              b_m(i,jsouth(i,j,k),k) = b_m(i,jsouth(i,j,k),k) - 2.0d0*a_m(i,jsouth(i,j,k),k,n)*bc_uw_g(l)

                              A_m(i,jsouth(i,j,k),k,n) = 0.0d0

                              A_m(I,J,K,:) = 0.0d0
                              A_m(I,J,K,0) = -1.0d0
                              b_m(I,J,K) = 0.0d0
                           ELSE
                              A_m(I,J,K,S) = -(HALF*BC_HW_G(L)-ODY)
                              A_m(I,J,K,0) = -(HALF*BC_HW_G(L)+ODY)
                              b_m(I,J,K) = -BC_HW_G(L)*BC_UW_G(L)
                           ENDIF





                        else if (1.eq.flag(i,j,ktop(i,j,k),1)) THEN
                           IF (IS_UNDEFINED(BC_HW_G(L))) THEN
                              A_m(I,J,K,T) = -HALF
                              A_m(I,J,K,0) = -HALF
                              b_m(I,J,K) = -BC_UW_G(L)
                           ELSE
                              A_m(I,J,K,0)=-(HALF*BC_HW_G(L)+ODZ)
                              A_m(I,J,K,T)=-(HALF*BC_HW_G(L)-ODZ)
                              b_m(I,J,K) = -BC_HW_G(L)*BC_UW_G(L)
                           ENDIF
                        else if (1.eq.flag(i,j,kbot(i,j,k),1)) THEN
                           IF (IS_UNDEFINED(BC_HW_G(L))) THEN
                              A_m(I,J,K,B) = -HALF
                              A_m(I,J,K,0) = -HALF
                              b_m(I,J,K) = -BC_UW_G(L)
                           ELSE
                              A_m(I,J,K,B) = -(HALF*BC_HW_G(L)-ODZ)
                              A_m(I,J,K,0) = -(HALF*BC_HW_G(L)+ODZ)
                              b_m(I,J,K) = -BC_HW_G(L)*BC_UW_G(L)
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO


! Setting p_inflow or p_outflow flow boundary conditions
! ---------------------------------------------------------------->>>
            ELSEIF (BC_TYPE(L)=='P_INFLOW' .OR. BC_TYPE(L)=='P_OUTFLOW') THEN
                IF (BC_PLANE(L) == 'W') THEN
! if the fluid cell is on the west side of the outflow/inflow boundary
! then set the velocity in the boundary cell equal to the velocity of
! the adjacent fluid cell
                  I1 = BC_I_W(L)
                  I2 = BC_I_E(L)
                  J1 = BC_J_S(L)
                  J2 = BC_J_N(L)
                  K1 = BC_K_B(L)
                  K2 = BC_K_T(L)
                  DO K = K1, K2
                     DO J = J1, J2
                        DO I = I1, I2
                           A_m(I,J,K,E) = ZERO
                           A_m(I,J,K,W) = ONE
                           A_m(I,J,K,N) = ZERO
                           A_m(I,J,K,S) = ZERO
                           A_m(I,J,K,T) = ZERO
                           A_m(I,J,K,B) = ZERO
                           A_m(I,J,K,0) = -ONE
                           b_m(I,J,K) = ZERO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
! end setting of p_inflow or p_otuflow flow boundary conditions
! ----------------------------------------------------------------<<<


! Setting bc that are defined but not nsw, fsw, psw, p_inflow,
! p_outflow, or outflow (at this time, this section addresses
! mass_inflow and mass_outflow type boundaries)
! ---------------------------------------------------------------->>>
            ELSE
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
               DO K = K1, K2
                  DO J = J1, J2
                     DO I = I1, I2
! setting the velocity in the boundary cell equal to what is known
                        A_m(I,J,K,E) = ZERO
                        A_m(I,J,K,W) = ZERO
                        A_m(I,J,K,N) = ZERO
                        A_m(I,J,K,S) = ZERO
                        A_m(I,J,K,T) = ZERO
                        A_m(I,J,K,B) = ZERO
                        A_m(I,J,K,0) = -ONE
                        b_m(I,J,K) = -U_G(I,J,K)
                        IF (BC_PLANE(L) == 'W') THEN
! if the fluid cell is on the west side of the outflow/inflow boundary
! then set the velocity in the adjacent fluid cell equal to what is
! known in that cell
                           A_m(iwest(i,j,k),j,k,E) = ZERO
                           A_m(iwest(i,j,k),j,k,W) = ZERO
                           A_m(iwest(i,j,k),j,k,N) = ZERO
                           A_m(iwest(i,j,k),j,k,S) = ZERO
                           A_m(iwest(i,j,k),j,k,T) = ZERO
                           A_m(iwest(i,j,k),j,k,B) = ZERO
                           A_m(iwest(i,j,k),j,k,0) = -ONE
                           b_m(iwest(i,j,k),j,k) = -U_G(iwest(i,j,k),j,k)
                           if (j.eq.lo(2)) print *,'5:SETTING B TO U_G ',iwest(i,j,k), b_m(iwest(i,j,k),j,k)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF   ! end if/else (bc_type)
                    ! ns, fs, psw; else
                    ! p_inflow, p_outflow, or outflow; else
! end setting of 'else' flow boundary conditions
! (mass_inflow/mass_outflow)
! ----------------------------------------------------------------<<<

         ENDIF   ! end if (bc_defined)
      ENDDO   ! end L do loop over dimension_bc

      END SUBROUTINE SOURCE_U_G_BC

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_U_G                                        C
!  Purpose: Adds point sources to the gas phase U-momentum equation.   C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE POINT_SOURCE_U_G(slo, shi, lo, hi, A_m, b_m, flag, dx, dy, dz)

      use ps, only: dimension_ps, ps_defined, ps_volume, ps_vel_mag_g, ps_massflow_g
      use ps, only: ps_u_g, ps_i_e, ps_i_w, ps_j_s, ps_j_n, ps_k_b, ps_k_t

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! Septadiagonal matrix A_m
      real(c_real), INTENT(IN   ) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      ! Vector b_m
      real(c_real), INTENT(INOUT) :: b_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer, intent(in   ) :: flag &
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), INTENT(IN   ) :: dx,dy,dz

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K
      INTEGER :: PSV, M
      INTEGER :: lIE, lIW
! terms of bm expression
      real(c_real) :: pSource
      real(c_real) :: vol
!-----------------------------------------------

      vol = dx*dy*dz

! Set reference phase to gas
      M = 0

! Calculate the mass going into each (i,j,k) cell. This is done for each
! call in case the point source is time dependent.
      PS_LP: do PSV = 1, DIMENSION_PS
         if(.NOT.PS_DEFINED(PSV)) cycle PS_LP
         if(abs(PS_U_g(PSV)) < small_number) cycle PS_LP

         if(PS_U_g(PSV) < 0.0d0) then
            lIW = PS_I_W(PSV) - 1
            lIE = PS_I_E(PSV) - 1
         else
            lIW = PS_I_W(PSV)
            lIE = PS_I_E(PSV)
         endif

         do k = PS_K_B(PSV), PS_K_T(PSV)
         do j = PS_J_S(PSV), PS_J_N(PSV)
         do i = lIW, lIE

            if(.NOT. 1.eq.flag(i,j,k,1)) cycle

            pSource =  PS_MASSFLOW_G(PSV) * (vol/PS_VOLUME(PSV))

            b_m(I,J,K) = b_m(I,J,K) - pSource *                        &
               PS_U_g(PSV) * PS_VEL_MAG_g(PSV)

         enddo
         enddo
         enddo
      enddo PS_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_U_G
end module source_u_g_module
