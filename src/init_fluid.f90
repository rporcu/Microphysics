module init_fluid_module
   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: init_fluid                                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine init_fluid(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, lo, hi, &
                         ep_g, ro_g, rop_g, p_g, u_g, v_g, w_g, &
                         mu_g, lambda_g, dx, dy, dz) &
      bind(C, name="init_fluid")

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      use calc_ro_g_module, only: calc_ro_g
      use calc_mu_g_module, only: calc_mu_g

      use param1, only: is_undefined, undefined
      use fld_const, only: ro_g0, mu_g0

      implicit none

! Dummy arguments .....................................................//
      integer(c_int), intent(in   ) :: slo(3), shi(3), lo(3), hi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)

      real(c_real), intent(inout) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(inout) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(inout) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(inout) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(inout) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: lambda_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: dx, dy, dz

      ! Set user specified initial conditions (IC)
      call set_ic(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, p_g, u_g, v_g, w_g)

      ! Set the initial pressure field
      call set_p_g(slo, shi, lo, hi, p_g, ep_g, dx, dy, dz)

      ! Set the initial fluid density
      if (is_undefined(ro_g0)) then
         call calc_ro_g(slo,shi,lo,hi,ro_g,rop_g,p_g,ep_g)
      else
         ro_g = ro_g0
         rop_g = ro_g0*ep_g
      endif

      ! Remove undefined values at wall cells for scalars
      where(rop_g == undefined) rop_g = 0.0

      ! Set the initial viscosity
      if (is_undefined(mu_g0)) then
         call calc_mu_g(slo,shi,lambda_g,mu_g)
      else
             mu_g(:,:,:) = mu_g0
         lambda_g(:,:,:) = -(2.0d0/3.0d0)*mu_g0
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
   subroutine set_ic(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, p_g, u_g, v_g, w_g)

      use geometry, only: domlo, domhi
      use ic, only: dimension_ic, ic_defined
      use ic, only: ic_i_w, ic_j_s, ic_k_b, ic_i_e, ic_j_n, ic_k_t
      use ic, only: ic_p_g, ic_u_g, ic_v_g, ic_w_g
      use scales, only: scale_pressure
      use param1, only: undefined, is_defined

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      IMPLICIT NONE

      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)

      real(c_real), intent(inout) ::  p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) ::  u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(inout) ::  v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(inout) ::  w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      integer :: i, j, k
      integer :: istart, iend
      integer :: jstart, jend
      integer :: kstart, kend
      ! local index for initial condition
      integer :: icv

      ! Temporary variables for storing IC values
      real(c_real) :: pgx, ugx, vgx, wgx

!  Set the initial conditions.
      do icv = 1, dimension_ic
         if (ic_defined(icv)) then

            ! Use the volume fraction already calculated from particle data
            pgx = ic_p_g(icv)
            ugx = ic_u_g(icv)
            vgx = ic_v_g(icv)
            wgx = ic_w_g(icv)

            if (is_defined(ugx)) then
               istart = max(ulo(1), ic_i_w(icv))
               jstart = max(ulo(2), ic_j_s(icv))
               kstart = max(ulo(3), ic_k_b(icv))
               iend   = min(uhi(1), ic_i_e(icv))
               jend   = min(uhi(2), ic_j_n(icv))
               kend   = min(uhi(3), ic_k_t(icv))
               u_g(istart:iend,jstart:jend,kstart:kend) = ugx
               if (ulo(1).lt.domlo(1)) &
                  u_g(ulo(1):istart-1,jstart:jend,kstart:kend) = ugx
               if (uhi(1).gt.domhi(1)) &
                  u_g(iend+1:uhi(1)  ,jstart:jend,kstart:kend) = ugx
            end if

            if (is_defined(ugx)) then
               istart = max(vlo(1), ic_i_w(icv))
               jstart = max(vlo(2), ic_j_s(icv))
               kstart = max(vlo(3), ic_k_b(icv))
               iend   = min(vhi(1), ic_i_e(icv))
               jend   = min(vhi(2), ic_j_n(icv))
               kend   = min(vhi(3), ic_k_t(icv))
               v_g(istart:iend,jstart:jend,kstart:kend) = vgx
               if (vlo(2).lt.domlo(2)) &
                  v_g(istart:iend,vlo(2):jstart-1,kstart:kend) = vgx
               if (vhi(2).gt.domhi(2)) &
                  v_g(istart:iend,jend+1:vhi(2)  ,kstart:kend) = vgx
            end if

            if (is_defined(ugx)) then
               istart = max(wlo(1), ic_i_w(icv))
               jstart = max(wlo(2), ic_j_s(icv))
               kstart = max(wlo(3), ic_k_b(icv))
               iend   = min(whi(1), ic_i_e(icv))
               jend   = min(whi(2), ic_j_n(icv))
               kend   = min(whi(3), ic_k_t(icv))
               w_g(istart:iend,jstart:jend,kstart:kend) = wgx
               if (wlo(3).lt.domlo(3)) &
                  w_g(istart:iend,jstart:jend,wlo(3):kstart-1) = wgx
               if (whi(3).gt.domhi(3)) &
                  w_g(istart:iend,jstart:jend,kend+1:whi(3)  ) = wgx
            end if

            istart = max(slo(1), ic_i_w(icv))
            jstart = max(slo(2), ic_j_s(icv))
            kstart = max(slo(3), ic_k_b(icv))
            iend   = min(shi(1), ic_i_e(icv))
            jend   = min(shi(2), ic_j_n(icv))
            kend   = min(shi(3), ic_k_t(icv))
            do k = kstart, kend
               do j = jstart, jend
                  do i = istart, iend
                        p_g(i,j,k) = merge(scale_pressure(pgx),&
                                      undefined, is_defined(pgx))
                  enddo
               enddo
            enddo

         endif
      enddo

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
      SUBROUTINE set_p_g(slo, shi, lo, hi, p_g, ep_g, dx, dy, dz)

      USE bc, only: delp_x, delp_y, delp_z
      USE bc, only: dimension_ic
      USE bc, only: dimension_bc, bc_type, bc_p_g, bc_defined, bc_plane
      USE compar, only: myPE
      USE constant , only: gravity
      USE eos, ONLY: EOSG
      USE fld_const, only: mw_avg, ro_g0
      USE geometry, only: xlength, ylength, zlength
      USE ic       , only: ic_p_g, ic_defined
      USE scales   , only: scale_pressure
      use exit_mod, only: mfix_exit
      use funits   , only: dmp_log, unit_log
      USE geometry, only: domhi

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int
      USE param1   , only: zero, undefined
      USE param1   , only: is_defined, is_undefined

      IMPLICIT NONE

      integer, intent(in) :: slo(3), shi(3), lo(3), hi(3)

      real(c_real), intent(inout) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: ep_g&
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
         pj = pj - dpodx*dx*(hi(1)-domhi(1)+2)
         do i = shi(1), slo(1), -1
            pj = pj + dpodx*dx
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  p_g(i,j,k) = scale_pressure(pj)
               enddo
            enddo
         enddo
      endif

      if (is_defined(delp_y)) then
         dpody = delp_y/ylength
         pj = pj - dpody*dy*(hi(2)-domhi(2)+1)
         do j = hi(2), lo(2), -1
            pj = pj + dpody*dy
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  p_g(i,j,k) = scale_pressure(pj)
               enddo
            enddo
         enddo
      endif

      if (is_defined(delp_z)) then
         dpodz = delp_z/zlength
         pj = pj - dpodz*dz*(hi(3)-domhi(3)+1)
         do k = hi(3), lo(3), -1
            pj = pj + dpodz*dz
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  p_g(i,j,k) = scale_pressure(pj)
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
            p_g = zero
            goto 100

         else   ! compressible case

! Error condition -- no pressure outflow boundary condition is specified
! if a case is compressible and pressure in any of the initial
! conditions regions is unspecified, then a PO is effectively required
! (i.e., is specifies a bc_p_g).
            if(dmp_log)write (unit_log, 1000)
            call mfix_exit(mype)
         endif
      endif

! Set an approximate pressure field assuming that the pressure drop
! balances the weight of the bed, if the initial pressure-field is not
! specified

      if(abs(gravity(1)) > epsilon(0.0d0)) then
         do i = hi(1), lo(1), -1

! Find the average weight per unit area over an x-z slice
            bed_weight = 0.0
            area = 0.0
            darea = dy*dz
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  area = area + darea
                  if (is_undefined(ro_g0)) then
                     bed_weight = bed_weight - dx*gravity(1)*&
                        ep_g(i,j,k)*eosg(mw_avg,pj,295.15d0)*darea
                  else
                     bed_weight = bed_weight - dx*gravity(1)*&
                        ep_g(i,j,k)*ro_g0*darea
                  endif
               enddo
            enddo

! Global Sum
            if (0.0 < abs(area)) bed_weight = bed_weight/area

            pj = pj + bed_weight
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  if(is_undefined(p_g(i,j,k))) p_g(i,j,k)=scale_pressure(pj)
               enddo
            enddo
         enddo

      else if(abs(gravity(3)) > epsilon(0.0d0)) then
         do k = hi(3), lo(3), -1

! Find the average weight per unit area over an x-z slice
            bed_weight = 0.0
            area = 0.0
            darea = dx*dy
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  area = area + darea
                  if (is_undefined(ro_g0)) then
                     bed_weight = bed_weight - dz*gravity(3)*&
                        ep_g(i,j,k)*eosg(mw_avg,pj,295.15d0)*darea
                  else
                     bed_weight = bed_weight - dz*gravity(3)*&
                        ep_g(i,j,k)*ro_g0*darea
                  endif
               enddo
            enddo

! Global Sum
            if (0.0 < abs(area)) bed_weight = bed_weight/area

            pj = pj + bed_weight
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  if(is_undefined(p_g(i,j,k))) p_g(i,j,k)=scale_pressure(pj)
               enddo
            enddo
         enddo

      else
         do j = hi(2), lo(2), -1

! Find the average weight per unit area over an x-z slice
            bed_weight = 0.0
            area = 0.0
            darea = dx*dz
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  area = area + darea
                  if (is_undefined(ro_g0)) then
                     bed_weight = bed_weight - dy*gravity(2)*ep_g(i,j,k)*&
                        eosg(mw_avg,pj,295.15d0)*darea
                  else
                     bed_weight = bed_weight - dy*gravity(2)*ep_g(i,j,k)*&
                        ro_g0*darea
                  endif
               enddo
            enddo

! Global Sum
         ! call global_all_sum(bed_weight)
         ! call global_all_sum(area)
            IF (0.0 < ABS(AREA)) BED_WEIGHT = BED_WEIGHT/AREA

            pj = pj + bed_weight
            do k = lo(3),hi(3)
               do i = lo(1),hi(1)
                  if (is_undefined(p_g(i,j,k)))&
                     p_g(i,j,k)=scale_pressure(pj)
               enddo
            enddo
         enddo
      endif
! end setting an undefined pressure in an initial condition region
! ----------------------------------------------------------------<<<

  100 CONTINUE

      RETURN

 1000 FORMAT(/1X,70('*')//' From: SET_FLUIDBED_P'/' Message: Outflow ',&
         'pressure boundary condition (P_OUTFLOW) not found.',/&
         'All the initial pressures (IC_P_g) or at least one P_OUTFLOW',/&
         'condition need to be specified',/1X,70('*')/)

   end subroutine set_p_g
end module init_fluid_module
