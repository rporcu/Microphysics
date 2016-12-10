module drag_gs_des0_module

  contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_DES_DRAG_GS                                        C
!  Purpose: This subroutine is only called from the DISCRETE side.     C
!     It performs the following functions:                             C
!     - If interpolated, then it calculates the fluid-particle         C
!       drag coefficient (F_GP) based on the particle velocity and     C
!       interpolated fluid velocity. This F_GP is used to calculate    C
!       the fluid-solids drag on the particle.                         C
!     - The total contact force on the particle is then updated to     C
!       include the gas-solids drag force and gas pressure force       C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DRAG_GS_DES0(ep_g,u_g,v_g,w_g,ro_g,mu_g,gradPg,des_radius,pvol,des_pos_new,des_vel_new,fc)

!-----------------------------------------------
! Modules
!-----------------------------------------------

      use compar        , only:  istart3, iend3, jstart3, jend3, kstart3, kend3
      use discretelement, only: xe, yn, zt, dimn, pic, pinc, interp_scheme
      use functions     , only: fluid_at,ip1,jp1,kp1
      use interpolation , only: set_interpolation_stencil, set_interpolation_scheme

      use des_drag_gp_module, only: des_drag_gp
      use discretelement, only: particle_state, entering_ghost, exiting_ghost, nonexistent, normal_ghost

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: des_radius, pvol

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

      DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: des_vel_new, des_pos_new
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: fc

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! general i, j, k indices
      INTEGER :: I, J, K
      INTEGER :: II, JJ, KK
! indices used for interpolation stencil (unclear why IE, JN, KTP are
! needed)
      INTEGER :: lli, llj, llk
      INTEGER :: IW, IE, JS, JN, KB, KTP
! i,j,k indices of the fluid cell the particle resides in minus 1
! (e.g., shifted 1 in west, south, bottom direction)
      INTEGER, DIMENSION(3) :: PCELL
! order of interpolation set in the call to set_interpolation_scheme
! unless it is re/set later through the call to set_interpolation_stencil
      INTEGER :: ONEW
! particle number index, used for looping
      INTEGER :: NP, nindx
! constant whose value depends on dimension of system
! avg_factor=0.125 (in 3D) or =0.25 (in 2D)
! avg_factor=0.250 (in 3D) or =0.50 (in 2D)
      DOUBLE PRECISION :: AVG_FACTOR

!Handan Liu added temporary variables on April 20 2012
          DOUBLE PRECISION, DIMENSION(2,2,2,3) :: gst_tmp,vst_tmp
          DOUBLE PRECISION, DIMENSION(2,2,2) :: weight_ft
          DOUBLE PRECISION :: velfp(3), desposnew(3)
          DOUBLE PRECISION :: D_FORCE(3)
          DOUBLE PRECISION :: VEL_NEW(3)
          DOUBLE PRECISION :: f_gp

! INTERPOLATED fluid-solids drag (the rest of this routine):
! Calculate the gas solids drag coefficient using the particle
! velocity and the fluid velocity interpolated to particle
! position.
!----------------------------------------------------------------->>>
! avg_factor=0.25 (in 3D) or =0.5 (in 2D)
      AVG_FACTOR = 0.25d0

! sets several quantities including interp_scheme, scheme, and
! order and allocates arrays necessary for interpolation
      call set_interpolation_scheme(2)

      DO llK = kstart3, kend3
      DO llJ = jstart3, jend3
      DO llI = istart3, iend3
         if(.not.fluid_at(lli,llj,llk) .or. pinc(lli,llj,llk).eq.0) cycle
         i = lli
         j = llj
         k = llk

! generally a particle may not exist in a ghost cell. however, if the
! particle is adjacent to the west, south or bottom boundary, then pcell
! may be assigned indices of a ghost cell which will be passed to
! set_interpolation_stencil
         pcell(1) = i-1
         pcell(2) = j-1
         pcell(3) = k-1

! setup the stencil based on the order of interpolation and factoring in
! whether the system has any periodic boundaries. sets onew to order.
         call set_interpolation_stencil(pcell,iw,ie,js,jn,kb,&
                                        ktp,interp_scheme,dimn,ordernew = onew)

! Compute velocity at grid nodes and set the geometric stencil
         DO k = 1, ONEW
            DO j = 1,onew
               DO i = 1,onew
                  ii = iw + i-1
                  jj = js + j-1
                  kk = kb + k-1

                  gst_tmp(i,j,k,1) = xe(ii)
                  gst_tmp(i,j,k,2) = yn(jj)
                  gst_tmp(i,j,k,3) = zt(kk)
                  vst_tmp(i,j,k,1) = avg_factor*&
                     (u_g(ii,jj,kk)+u_g(ii,jp1(jj),kk))
                  vst_tmp(i,j,k,2) = avg_factor*&
                     (v_g(ii,jj,kk)+v_g(ip1(ii),jj,kk))

                  vst_tmp(i,j,k,1) = vst_tmp(i,j,k,1) + avg_factor*&
                     (u_g(ii,jj,kp1(kk)) + u_g(ii,jp1(jj),kp1(kk)))
                  vst_tmp(i,j,k,2) = vst_tmp(i,j,k,2) + avg_factor*&
                     (v_g(ii,jj,kp1(kk)) + v_g(ip1(ii),jj,kp1(kk)))
                  vst_tmp(i,j,k,3) = avg_factor*(w_g(ii,jj,kk)+&
                     w_g(ii,jp1(jj),kk)+w_g(ip1(ii),jj,kk)+&
                     w_g(ip1(ii),jp1(jj),kk))
               ENDDO
            ENDDO
         ENDDO
! loop through particles in the cell
! interpolate the fluid velocity (VELFP) to the particle's position.
         DO nindx = 1,PINC(lli,llj,llk)
            NP = PIC(lli,llj,llk)%p(nindx)
! skipping indices that do not represent particles and ghost particles
            if(nonexistent==particle_state(np) .or. &
               normal_ghost==particle_state(np) .or. &
               entering_ghost==particle_state(np) .or. &
               exiting_ghost==particle_state(np)) cycle

            desposnew(:) = des_pos_new(np,:)
            call DRAG_INTERPOLATION(gst_tmp,vst_tmp,desposnew,velfp,weight_ft)

! Calculate the particle centered drag coefficient (F_GP) using the
! particle velocity and the interpolated gas velocity.  Note F_GP
! obtained from des_drag_gp subroutine is given as:
!    f_gp=beta*vol_p/eps where vol_p is the particle volume.
! The drag force on each particle is equal to:
!    beta(u_g-u_s)*vol_p/eps.
! Therefore, the drag force = f_gp*(u_g - u_s)
            VEL_NEW(:) = DES_VEL_NEW(NP,:)
            CALL DES_DRAG_GP(NP, VEL_NEW, VELFP, EP_G(llI,llJ,llK), ro_g, mu_g,f_gp, des_radius, pvol)

! Calculate the gas-solids drag force on the particle
            D_FORCE(1:3) = F_GP*(VELFP-VEL_NEW)

! Update the contact forces (FC) on the particle to include gas
! pressure and gas-solids drag
            FC(NP,:) = FC(NP,:) + D_FORCE(:3)

! gradPg is evaluated as -dp/dx
            FC(NP,:) = FC(NP,:) + gradPg(lli,llj,llk,:)*PVOL(NP)
         ENDDO       ! end do (nindx = 1,pinc)

      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE DRAG_GS_DES0

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DES_DRAG_GS                                             C
!  Purpose: This subroutine is only called from the CONTINUUM side.    C
!     It performs the following functions:                             C
!     - If interpolated then, it calculates the fluid-particle         C
!       drag coefficient (F_GP) based on the particle velocity and     C
!       interpolated fluid velocity. It then determines the            C
!       the contributions of fluid-particle drag to the center         C
!       coefficient of the A matrix and the b (source) vector in the   C
!       matrix equation (A*VELFP=b) equation for the fluid phase       C
!       x, y and z momentum balances using F_GP.                       C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DRAG_GS_GAS0(ep_g, u_g, v_g, w_g, ro_g, mu_g,&
         f_gds, drag_am, drag_bm, des_radius, pvol, des_pos_new, des_vel_new, particle_phase)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use compar        , only:  istart2, iend2, jstart2, jend2, kstart2, kend2
      use compar        , only:  istart3, iend3, jstart3, jend3, kstart3, kend3

      use discretelement, only: xe, yn, zt, dimn, pic, pinc, &
                                interp_scheme, des_vol_node
      use interpolation , only: set_interpolation_stencil, set_interpolation_scheme
      use param1  , only: zero, one
      use functions     , only: fluid_at,ip1,jp1,kp1
      use mpi_node_des, only: des_addnodevalues

      use des_drag_gp_module, only: des_drag_gp
      use discretelement, only: particle_state, nonexistent, &
         normal_ghost, entering_ghost, exiting_ghost
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: des_radius, pvol

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
      DOUBLE PRECISION, INTENT(OUT  ) :: f_gds&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(OUT  ) :: drag_am&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(OUT  ) :: drag_bm&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)

      DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: des_vel_new, des_pos_new
      INTEGER, DIMENSION(:), INTENT(IN) :: particle_phase

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! local variable used for debugging
      LOGICAL :: FOCUS
! general i, j, k indices
      INTEGER :: I, J, K
      INTEGER :: II, JJ, KK
! indices used for interpolation stencil (unclear why IE, JN, KTP are
! needed)
      INTEGER :: IW, IE, JS, JN, KB, KTP
! i,j,k indices of the fluid cell the particle resides in minus 1
! (e.g., shifted 1 in west, south, bottom direction)
      INTEGER, DIMENSION(3) :: PCELL
! order of interpolation set in the call to set_interpolation_scheme
! unless it is re/set later through the call to set_interpolation_stencil
      INTEGER :: ONEW
! index of solid phase that particle NP belongs to
      INTEGER :: M, lli, llj, llk
! particle number index, used for looping
      INTEGER :: NP, nindx
! one over the volume of fluid cell
      DOUBLE PRECISION :: OVOL
! volume of fluid cell particle resides in
      DOUBLE PRECISION :: VCELL
! constant whose value depends on dimension of system
! avg_factor=0.125 (in 3D) or =0.25 (in 2D)
! avg_factor=0.250 (in 3D) or =0.50 (in 2D)
      DOUBLE PRECISION :: AVG_FACTOR
!Handan Liu added temporary variables on April 20 2012
      DOUBLE PRECISION, DIMENSION(2,2,2,3) :: gst_tmp,vst_tmp
      DOUBLE PRECISION, DIMENSION(2,2,2) :: weight_ft
      DOUBLE PRECISION :: velfp(3), desposnew(3)
      DOUBLE PRECISION :: VEL_NEW(3)
      DOUBLE PRECISION :: f_gp

!-----------------------------------------------



! INTERPOLATED fluid-solids drag (the rest of this routine):
! Calculate the fluid solids drag coefficient using the particle
! velocity and the fluid velocity interpolated to particle
! position.
!----------------------------------------------------------------->>>
! initializations
      drag_am = ZERO
      drag_bm = ZERO

! avg_factor=0.25 (in 3D) or =0.5 (in 2D)
      AVG_FACTOR = 0.25d0

! sets several quantities including interp_scheme, scheme, and
! order and allocates arrays necessary for interpolation
      call set_interpolation_scheme(2)

      DO llK = kstart3, kend3
      DO llJ = jstart3, jend3
      DO llI = istart3, iend3
         IF(.NOT.fluid_at(lli,llj,llk) .OR. PINC(lli,llj,llk)==0) cycle
         i = lli
         j = llj
         k = llk

! generally a particle may not exist in a ghost cell. however, if the
! particle is adjacent to the west, south or bottom boundary, then pcell
! may be assigned indices of a ghost cell which will be passed to
! set_interpolation_stencil
         pcell(1) = i-1
         pcell(2) = j-1
         pcell(3) = k-1

! setup the stencil based on the order of interpolation and factoring in
! whether the system has any periodic boundaries. sets onew to order.
         call set_interpolation_stencil(pcell,iw,ie,js,jn,kb,&
              ktp,interp_scheme,dimn,ordernew = onew)

! Compute velocity at grid nodes and set the geometric stencil
         DO k = 1, ONEW
            DO j = 1,onew
               DO i = 1,onew
                  ii = iw + i-1
                  jj = js + j-1
                  kk = kb + k-1
                  GST_TMP(I,J,K,1) = XE(II)
                  GST_TMP(I,J,K,2) = YN(JJ)
                  GST_TMP(I,J,K,3) = ZT(KK)
                  VST_TMP(I,J,K,1) = AVG_FACTOR*&
                     (U_G(ii,jj,kk)+U_G(ii,jp1(jj),kk))
                  VST_TMP(I,J,K,2) = AVG_FACTOR*&
                     (V_G(ii,jj,kk)+V_G(ip1(ii),jj,kk))

                  VST_TMP(I,J,K,1) = VST_TMP(I,J,K,1) + AVG_FACTOR*&
                     (U_G(Ii,jj,kp1(kk)) + U_G(ii,jp1(jj),kp1(kk)))

                  VST_TMP(I,J,K,2) = VST_TMP(I,J,K,2) + AVG_FACTOR*&
                     (V_G(ii,jj,kp1(kk)) + V_G(ip1(ii),jj,kp1(kk)))

                  VST_TMP(I,J,K,3) = AVG_FACTOR*(W_G(ii,jj,kk)+&
                     W_G(ii,jp1(jj),kk)+W_G(IP1(ii),jj,kk)+&
                     W_G(IP1(ii),JP1(jj),kK))
               ENDDO
            ENDDO
         ENDDO
! loop through particles in the cell
! interpolate the fluid velocity (VELFP) to the particle's position.
         DO nindx = 1,PINC(lli,llj,llk)
            NP = PIC(lli,llj,llk)%p(nindx)
! skipping indices that do not represent particles and ghost particles
            if(nonexistent==particle_state(np)) cycle
            if(normal_ghost==particle_state(np) .or. &
               entering_ghost==particle_state(np) .or. &
               exiting_ghost==particle_state(np)) cycle
            desposnew(:) = des_pos_new(np,:)
            call DRAG_INTERPOLATION(gst_tmp,vst_tmp,desposnew,velfp,weight_ft)
!
! Calculate the particle centered drag coefficient (F_GP) using the
! particle velocity and the interpolated gas velocity.  Note F_GP
! obtained from des_drag_gp subroutine is given as:
!    f_gp=beta*vol_p/eps where vol_p is the particle volume.
! The drag force on each particle is equal to:
!    beta(u_g-u_s)*vol_p/eps.
! Therefore, the drag force = f_gp*(u_g - u_s)
            VEL_NEW(:) = DES_VEL_NEW(NP,:)
            CALL DES_DRAG_GP(NP, VEL_NEW, VELFP, EP_G(llI,llJ,llK), ro_g, mu_g, f_gp, des_radius, pvol)

!-----------------------------------------------------------------<<<
! Calculate the corresponding gas solids drag force that is used in
! the gas phase momentum balances.
!----------------------------------------------------------------->>>
            focus = .false.
            M = particle_phase(np)

            DO k = 1, ONEW
               DO j = 1, onew
                  DO i = 1, onew
! shift loop index to new variables for manipulation
                     ii = iw + i-1
                     jj = js + j-1
                     kk = kb + k-1
! The interpolation is done using node. so one should use consistent
! numbering system. in the current version imap_c is used instead of
! iplus or iminus

! Replacing the volume of cell to volume at the node
                     vcell = des_vol_node(ii, jj, kk)
                     ovol = one/vcell

                     drag_am(ii, jj, kk) = drag_am(ii, jj, kk) + &
                        f_gp*weight_ft(i,j,k)*ovol

                     drag_bm(ii, jj, kk,1:3) = &
                        drag_bm(ii, jj, kk,1:3) + &
                        f_gp * vel_new(1:3) * &
                        weight_ft(i,j,k)*ovol
                  ENDDO
               ENDDO
            ENDDO
         ENDDO   ! end do (nindx = 1,pinc)

      ENDDO
      ENDDO
      ENDDO


! At the interface drag_am and drag_bm have to be added
! send recv will be called and the node values will be added
! at the junction. drag_am are drag_bm are altered by the
! routine when periodic boundaries are invoked. so both
! quantities are needed at the time of this call.
      call des_addnodevalues(drag_am, drag_bm)
!-----------------------------------------------------------------<<<
! Calculate/update the cell centered drag coefficient F_GDS for use
! in the pressure correction equation
!----------------------------------------------------------------->>>
! avg_factor=0.125 (in 3D) or =0.25 (in 2D)
      AVG_FACTOR = 0.125D0

      DO K = kstart3, kend3
         DO J = jstart3, jend3
            DO I = istart3, iend3

               if(fluid_at(i,j,k)) then
                  if (i.lt.istart2 .or. i.gt.iend2) cycle
                  if (j.lt.jstart2 .or. j.gt.jend2) cycle
                  if (k.lt.kstart2 .or. k.gt.kend2) cycle

                  f_gds(i,j,k) = avg_factor*(&
                     drag_am(i,j,k)      + drag_am(i,j-1,k)   + &
                     drag_am(i-1,j-1,k)  + drag_am(i-1,j,k) + &
                     drag_am(i,j,k-1)    + drag_am(i,j-1,k-1) + &
                     drag_am(i-1,j-1,k-1)+ drag_am(i-1,j,k-1))
               endif
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE DRAG_GS_GAS0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!  Subroutine: DRAG_INTERPOLATION                                       C
!  Purpose: DES - Calculate the fluid velocity interpolated at the      C
!           particle's location and weights. Replace 'interpolator'     C
!                       interface for OpenMP implementation.            C
!                                                                       C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_INTERPOLATION(GSTEN,VSTEN,DESPOS,VELFP,WEIGHTFACTOR)

        IMPLICIT NONE

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
        DOUBLE PRECISION, DIMENSION(2,2,2,3), INTENT(IN):: GSTEN
        DOUBLE PRECISION, DIMENSION(2,2,2,3), INTENT(IN):: VSTEN
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: DESPOS
        DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: VELFP
        DOUBLE PRECISION, DIMENSION(2,2,2), INTENT(OUT) :: WEIGHTFACTOR
        INTEGER :: II, JJ, KK

        DOUBLE PRECISION, DIMENSION(2) :: XXVAL, YYVAL, ZZVAL
        DOUBLE PRECISION :: DXX, DYY, DZZ
        DOUBLE PRECISION, DIMENSION(3) :: ZETAA

        DXX = GSTEN(2,1,1,1) - GSTEN(1,1,1,1)
        DYY = GSTEN(1,2,1,2) - GSTEN(1,1,1,2)

        ZETAA(1:2) = DESPOS(1:2) - GSTEN(1,1,1,1:2)

        ZETAA(1) = ZETAA(1)/DXX
        ZETAA(2) = ZETAA(2)/DYY

        XXVAL(1)=1-ZETAA(1)
        YYVAL(1)=1-ZETAA(2)
        XXVAL(2)=ZETAA(1)
        YYVAL(2)=ZETAA(2)

        VELFP(:) = 0.D0

        DZZ = GSTEN(1,1,2,3) - GSTEN(1,1,1,3)
        ZETAA(3) = DESPOS(3) - GSTEN(1,1,1,3)
        ZETAA(3) = ZETAA(3)/DZZ
        ZZVAL(1)=1-ZETAA(3)
        ZZVAL(2)=ZETAA(3)
        DO KK=1,2
           DO JJ=1,2
              DO II=1,2
                 WEIGHTFACTOR(II,JJ,KK) = XXVAL(II)*YYVAL(JJ)*ZZVAL(KK)
                 VELFP(1:3) = VELFP(1:3) + VSTEN(II,JJ,KK,1:3)*WEIGHTFACTOR(II,JJ,KK)
              ENDDO
           ENDDO
        ENDDO

      END SUBROUTINE DRAG_INTERPOLATION

end module drag_gs_des0_module
