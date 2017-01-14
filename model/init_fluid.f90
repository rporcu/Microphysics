!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: init_fluid                                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine init_fluid(slo, shi, lo, hi, ep_g, ro_g, rop_g, p_g, u_g, v_g, w_g, &
                         mu_g, lambda_g, flag, dx, dy, dz) &
      bind(C, name="init_fluid")

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

      use calc_ro_g_module, only: calc_ro_g
      use calc_mu_g_module, only: calc_mu_g

      use param1, only: is_undefined, undefined
      use fld_const, only: ro_g0, mu_g0

      implicit none

! Dummy arguments .....................................................//
      integer(c_int), intent(in) :: slo(3), shi(3), lo(3), hi(3)

      real(c_real), intent(inout) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: lambda_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(in   ) :: dx, dy, dz

! Local variables .....................................................//
      real(c_real), parameter :: f2o3 = 2.d0/3.d0

!-----------------------------------------------------------------------!

      ! Set user specified initial conditions (IC)
      call set_ic(slo, shi, lo, hi, p_g, u_g, v_g, w_g, flag)

      ! Set the initial pressure field
      call set_p_g(slo, shi, lo, hi, p_g, ep_g, flag(:,:,:,1), dx, dy, dz)

      ! Set the initial fluid density
      if (is_undefined(ro_g0)) then
         call calc_ro_g(slo, shi, ro_g,rop_g,p_g,ep_g,flag)
      else
         where (flag(:,:,:,1) < 100) ro_g = ro_g0
         where (flag(:,:,:,1) < 100) rop_g = ro_g0*ep_g
      endif

      ! Remove undefined values at wall cells for scalars
      where(rop_g == undefined) rop_g = 0.0

      ! Set the initial viscosity
      if (is_undefined(ro_g0)) then
         call calc_mu_g(slo,shi,lambda_g,mu_g,flag)
      else
         where (flag(:,:,:,1) == 1) mu_g = mu_g0
         where (flag(:,:,:,1) == 1) lambda_g = -(2.0d0/3.0d0)*mu_g0
      endif

   end subroutine init_fluid

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_IC                                                  !
!  Author: M. Syamlal                                 Date: 21-JAN-92  !
!                                                                      !
!  Purpose: This module sets all the initial conditions.               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine set_ic(slo, shi, lo, hi, p_g, u_g, v_g, w_g, flag)

      use ic, only: dimension_ic, ic_type, ic_defined
      use ic, only: ic_i_w, ic_j_s, ic_k_b, ic_i_e, ic_j_n, ic_k_t
      use ic, only: ic_ep_g, ic_p_g, ic_u_g, ic_v_g, ic_w_g
      use scales, only: scale_pressure
      use param1, only: undefined, is_defined
      use run, only: dem_solids

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

      IMPLICIT NONE

      integer, intent(in) :: slo(3), shi(3), lo(3), hi(3)

      real(c_real), intent(inout) ::  p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) ::  u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) ::  v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) ::  w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer, intent(in) ::  flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      integer :: i, j, k
! local index for initial condition
      integer :: l
! Temporary variables for storing IC values
      real(c_real) :: pgx, ugx, vgx, wgx
!-----------------------------------------------

!  Set the initial conditions.
      do l = 1, dimension_ic
         if (ic_defined(l)) then

! Use the volume fraction already calculated from particle data
            pgx = ic_p_g(l)
            ugx = ic_u_g(l)
            vgx = ic_v_g(l)
            wgx = ic_w_g(l)

            do k = ic_k_b(l), ic_k_t(l)
               do j = ic_j_s(l), ic_j_n(l)
                  do i = ic_i_w(l), ic_i_e(l)

                     if (flag(i,j,k,1)<100) then

                        p_g(i,j,k) = merge(scale_pressure(pgx),&
                           undefined, is_defined(pgx))

                        if (is_defined(ugx)) u_g(i,j,k) = ugx
                        if (is_defined(vgx)) v_g(i,j,k) = vgx
                        if (is_defined(wgx)) w_g(i,j,k) = wgx

                     endif
                  enddo
               enddo
            enddo
         endif
      enddo

      return
   end subroutine set_ic

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: set_p_g                                                 !
!  Author: M. Syamlal                                 Date: 21-JAN-92  !
!                                                                      !
!  Purpose: Set the pressure field inside the bed assuming gravity     !
!           is acting in the negative y-direction.                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE set_p_g(slo, shi, lo, hi, p_g, ep_g, flag, dx, dy, dz)

      USE bc, only: delp_x, delp_y, delp_z
      USE bc, only: dimension_ic
      USE bc, only: dimension_bc, bc_type, bc_p_g, bc_defined
      USE compar, only: myPE
      USE constant , only: gravity
      USE eos, ONLY: EOSG
      USE fld_const, only: mw_avg, ro_g0
      USE geometry, only: domlo,domhi
      USE geometry, only: xlength, ylength, zlength
      USE ic       , only: ic_p_g, ic_defined
      USE param1   , only: is_defined, zero, undefined, is_undefined
      USE scales   , only: scale_pressure
      use exit_mod, only: mfix_exit
      use funits   , only: dmp_log, unit_log

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

      IMPLICIT NONE

      integer, intent(in) :: slo(3), shi(3), lo(3), hi(3)

      real(c_real), intent(inout) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: dx, dy, dz
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      INTEGER :: I, J, K
! Local loop counter
      INTEGER :: L
! Gas pressure at the axial location j
      real(c_real) :: PJ
! Bed weight per unit area
      real(c_real) :: BED_WEIGHT
! Total area of a x-z plane
      real(c_real) :: AREA
! x-z plane area of one cell
      real(c_real) :: dAREA
! Average pressure drop per unit length
      real(c_real) :: DPoDX, DPoDY, DPoDZ
!-----------------------------------------------

! If any initial pressures are unspecified skip next section
! calculations.
      do l = 1, dimension_ic
         if (ic_defined(l)) then
            if (is_undefined(ic_p_g(l))) goto 60
            pj = ic_p_g(l)
         endif
      enddo

! Here the pressure in each cell is determined from a specified pressure
! drop across the domain length. This section requires that the pressure
! is already defined in all initial condition regions (otherwise this
! section would be skipped)
! ---------------------------------------------------------------->>>
      if (is_defined(delp_x)) then
         dpodx = delp_x/xlength
         pj = pj - dpodx*dx
         do i = domhi(1), domlo(1), -1
            pj = pj + dpodx*dx
            do k = domlo(3), domhi(3)
               do j = domlo(2), domhi(2)
                  if (flag(i,j,k)==1) p_g(i,j,k) = scale_pressure(pj)
               enddo
            enddo
         enddo
      endif

      if (is_defined(delp_y)) then
         dpody = delp_y/ylength
         pj = pj - dpody*dy
         do j = domhi(2), domlo(2), -1
            pj = pj + dpody*dy
            do k = domlo(3), domhi(3)
               do i = domlo(1), domhi(1)
                  if (flag(i,j,k)==1) p_g(i,j,k) = scale_pressure(pj)
               enddo
            enddo
         enddo
      endif

      if (is_defined(delp_z)) then
         dpodz = delp_z/zlength
         pj = pj - dpodz*dz
         do k = domhi(3), domlo(3), -1
            pj = pj + dpodz*dz
            do j = domlo(2), domhi(2)
               do i = domlo(1), domhi(1)
                  if (flag(i,j,k)==1) p_g(i,j,k) = scale_pressure(pj)
               enddo
            enddo
         enddo
      endif
! ----------------------------------------------------------------<<<
      GOTO 100   ! pressure in all intial condition region cells was defined

   60 CONTINUE   ! pressure in an initial condition region cell was undefined


! ---------------------------------------------------------------->>>
! Search for an outflow boundary condition where pressure is specified
      pj = undefined
      do l = 1, dimension_bc
         if (bc_defined(l) .and. bc_type(l)=='P_OUTFLOW') pj = bc_p_g(l)
      enddo

      if (is_undefined(pj)) then
! either a PO was not specified and/or a PO was specified but not the
! pressure at the outlet
         if (is_defined(ro_g0)) then

            ! If incompressible flow set P_g to zero
            where (flag .eq. 1) p_g = zero
            goto 100

         else   ! compressible case

! Error condition -- no pressure outflow boundary condition is specified
! if a case is compressible and pressure in any of the initial
! conditions regions is unspecified, then a PO is effectively required
! (i.e., is specifies a bc_p_g).
            IF(DMP_LOG)WRITE (UNIT_LOG, 1000)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF


! Set an approximate pressure field assuming that the pressure drop
! balances the weight of the bed, if the initial pressure-field is not
! specified
      do j = domhi(2)+1, domlo(2), -1

! Find the average weight per unit area over an x-z slice
         bed_weight = 0.0
         area = 0.0
         do k = domlo(3), domhi(3)
            do i = domlo(1), domhi(1)
               if (flag(i,j,k)==1) then
                  darea = dx*dz
                  area = area + darea
                  if (is_undefined(ro_g0)) then
                     bed_weight = bed_weight - dy*gravity(2)*ep_g(i,j,k)*eosg(&
                        mw_avg,pj,295.15d0)*darea
                  else
                     bed_weight = bed_weight - dy*gravity(2)*ep_g(i,j,k)*ro_g0&
                        *darea
                  endif
               endif
            enddo
         enddo

! Global Sum
         ! call global_all_sum(bed_weight)
         ! call global_all_sum(area)
         IF (0.0 < ABS(AREA)) BED_WEIGHT = BED_WEIGHT/AREA

         PJ = PJ + BED_WEIGHT
         DO K = domlo(3),domhi(3)
            DO I = domlo(1),domhi(1)
               IF(flag(i,j,k) ==1 .AND. IS_UNDEFINED(P_G(I,J,K)))&
                  P_G(I,J,K)=SCALE_PRESSURE(PJ)
            ENDDO
         ENDDO
      ENDDO   ! end do (j=jmax2,jimn1, -1)
! end setting an undefined pressure in an initial condition region
! ----------------------------------------------------------------<<<

  100 CONTINUE

      ! call send_recv(P_G,2)

      RETURN

 1000 FORMAT(/1X,70('*')//' From: SET_FLUIDBED_P'/' Message: Outflow ',&
         'pressure boundary condition (P_OUTFLOW) not found.',/&
         'All the initial pressures (IC_P_g) or at least one P_OUTFLOW',/&
         'condition need to be specified',/1X,70('*')/)

   end subroutine set_p_g
