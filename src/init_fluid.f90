module init_fluid_module
   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: init_fluid                                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine init_fluid(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, lo, hi, &
                         domlo, domhi, ep_g, ro_g, rop_g, p_g, u_g, v_g, w_g, &
                         mu_g, lambda_g, dx, dy, dz, xlength, ylength, zlength) &
      bind(C, name="init_fluid")

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      use calc_ro_g_module, only: calc_ro_g
      use calc_mu_g_module, only: calc_mu_g

      use param, only: is_undefined, undefined

      implicit none

! Dummy arguments .....................................................//
      integer(c_int), intent(in   ) ::  lo(3),  hi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)

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
      real(c_real), intent(in   ) :: xlength, ylength, zlength

      ! Set user specified initial conditions (IC)
      call set_ic(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, &
                  domlo, domhi, dx, dy, dz, p_g, u_g, v_g, w_g)

      ! Set the initial pressure field
      call set_p_g(slo, shi, lo, hi, p_g, dx, dy, dz, &
                   xlength, ylength, zlength, domlo, domhi)

      ! Set the initial fluid density and viscosity
      call calc_ro_g(slo, shi, lo, hi, ro_g, rop_g, p_g, ep_g)

      call calc_mu_g(slo, shi, lo, hi, mu_g, lambda_g)

   end subroutine init_fluid

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: init_fluid_restart                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine init_fluid_restart(slo, shi, lo, hi, mu_g, lambda_g) &
      bind(C, name="init_fluid_restart")

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      use calc_ro_g_module, only: calc_ro_g
      use calc_mu_g_module, only: calc_mu_g

      implicit none

! Dummy arguments .....................................................//
      integer(c_int), intent(in   ) ::  lo(3),  hi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)

      real(c_real), intent(inout) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: lambda_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      call calc_mu_g(slo, shi, lo, hi, mu_g, lambda_g)

    end subroutine init_fluid_restart

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_IC                                                  !
!  Author: M. Syamlal                                 Date: 21-JAN-92  !
!                                                                      !
!  Purpose: This module sets all the initial conditions.               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine set_ic(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, &
                     domlo, domhi, dx, dy, dz, p_g, u_g, v_g, w_g)

      use ic, only: dim_ic, ic_defined
      use ic, only: ic_p_g, ic_u_g, ic_v_g, ic_w_g
      use ic, only: ic_x_e, ic_y_n, ic_z_t
      use ic, only: ic_x_w, ic_y_s, ic_z_b
      use scales, only: scale_pressure
      use param, only: undefined, is_defined

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      use calc_cell_module, only: calc_cell_ic

      implicit none

      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)
      real(c_real), intent(in   ) :: dx, dy, dz

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
      integer :: istart, iend
      integer :: jstart, jend
      integer :: kstart, kend
      ! local index for initial condition
      integer :: icv

      ! Temporary variables for storing IC values
      real(c_real) :: pgx, ugx, vgx, wgx, pval

      integer :: i_w, j_s, k_b
      integer :: i_e, j_n, k_t

!  Set the initial conditions.
      do icv = 1, dim_ic
         if (ic_defined(icv)) then

            call calc_cell_ic(dx, dy, dz, &
              ic_x_w(icv), ic_y_s(icv), ic_z_b(icv), &
              ic_x_e(icv), ic_y_n(icv), ic_z_t(icv), &
              i_w, i_e, j_s, j_n, k_b, k_t)

            ! Use the volume fraction already calculated from particle data
            pgx = ic_p_g(icv)
            ugx = ic_u_g(icv)
            vgx = ic_v_g(icv)
            wgx = ic_w_g(icv)

            if (is_defined(ugx)) then
               istart = max(ulo(1), i_w)
               jstart = max(ulo(2), j_s)
               kstart = max(ulo(3), k_b)
               iend   = min(uhi(1), i_e)
               jend   = min(uhi(2), j_n)
               kend   = min(uhi(3), k_t)
               u_g(istart:iend,jstart:jend,kstart:kend) = ugx
               if (ulo(1).lt.domlo(1)) &
                  u_g(ulo(1):istart-1,jstart:jend,kstart:kend) = ugx
               if (uhi(1).gt.domhi(1)) &
                  u_g(iend+1:uhi(1)  ,jstart:jend,kstart:kend) = ugx
            end if

            if (is_defined(vgx)) then
               istart = max(vlo(1), i_w)
               jstart = max(vlo(2), j_s)
               kstart = max(vlo(3), k_b)
               iend   = min(vhi(1), i_e)
               jend   = min(vhi(2), j_n)
               kend   = min(vhi(3), k_t)
               v_g(istart:iend,jstart:jend,kstart:kend) = vgx
               if (vlo(2).lt.domlo(2)) &
                  v_g(istart:iend,vlo(2):jstart-1,kstart:kend) = vgx
               if (vhi(2).gt.domhi(2)) &
                  v_g(istart:iend,jend+1:vhi(2)  ,kstart:kend) = vgx
            end if

            if (is_defined(wgx)) then
               istart = max(wlo(1), i_w)
               jstart = max(wlo(2), j_s)
               kstart = max(wlo(3), k_b)
               iend   = min(whi(1), i_e)
               jend   = min(whi(2), j_n)
               kend   = min(whi(3), k_t)
               w_g(istart:iend,jstart:jend,kstart:kend) = wgx
               if (wlo(3).lt.domlo(3)) &
                  w_g(istart:iend,jstart:jend,wlo(3):kstart-1) = wgx
               if (whi(3).gt.domhi(3)) &
                  w_g(istart:iend,jstart:jend,kend+1:whi(3)  ) = wgx
            end if

            istart = max(slo(1), i_w)
            jstart = max(slo(2), j_s)
            kstart = max(slo(3), k_b)
            iend   = min(shi(1), i_e)
            jend   = min(shi(2), j_n)
            kend   = min(shi(3), k_t)
            pval = merge(scale_pressure(pgx),undefined,is_defined(pgx))
            p_g(istart:iend,jstart:jend,kstart:kend) = pval

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
      subroutine set_p_g(slo, shi, lo, hi, p_g, dx, dy, dz, &
                         xlength, ylength, zlength, domlo, domhi)

      use bc, only: delp_x, delp_y, delp_z
      use bc, only: dim_bc, bc_type, bc_p_g, bc_defined
      use constant , only: gravity
      use eos, ONLY: EOSG
      use fld_const, only: mw_avg, ro_g0
      use ic       , only: ic_p_g, ic_defined
      use scales   , only: scale_pressure


      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int
      use param   , only: zero, undefined
      use param   , only: is_defined, is_undefined
      use param, only: dim_ic

      implicit none

      integer, intent(in) :: slo(3), shi(3), lo(3), hi(3)
      integer, intent(in) :: domlo(3), domhi(3)

      real(c_real), intent(inout) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: dx, dy, dz
      real(c_real), intent(in   ) :: xlength, ylength, zlength
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      integer :: I, J, K
! Local loop counter
      integer :: L
! Gas pressure at the axial location j
      real(c_real) :: PJ
! Average pressure drop per unit length
      real(c_real) :: DPoDX, DPoDY, DPoDZ
!-----------------------------------------------

! If any initial pressures are unspecified skip next section
! calculations.
      do l = 1, dim_ic
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
      if (abs(delp_x) > epsilon(zero)) then
         dpodx = delp_x/xlength
         pj = pj - dpodx*dx*(hi(1)-domhi(1)+1)
         do i = hi(1), lo(1), -1
            pj = pj + dpodx*dx
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  p_g(i,j,k) = scale_pressure(pj)
               enddo
            enddo
         enddo
      endif

      if (abs(delp_y) > epsilon(zero)) then
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

      if (abs(delp_z) > epsilon(zero)) then
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
      do l = 1, dim_bc
         if (bc_defined(l)) then
            if(bc_type(l)=='P_OUTFLOW' .or. bc_type(l)=='PO') &
               pj = bc_p_g(l)
         endif
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
            write (*, 1000)
            stop 20014
         endif
      endif

! Set an approximate pressure field assuming that the pressure drop
! balances the weight of the bed, if the initial pressure-field is not
! specified

      if (abs(gravity(1)) > epsilon(0.0d0)) then

         ! Find the average weight per unit area over an x-z slice
         if (is_undefined(ro_g0)) then
            dpodx = -gravity(1)*eosg(mw_avg,pj,295.15d0)
         else
            dpodx = -gravity(1)*ro_g0
         endif

         if (gravity(1) <= 0.0d0) then
            do i = domhi(1)+1, domlo(1), -1
               if (i <= hi(1)+1 .and. i >= lo(1)) &
                  p_g(i,:,:) = scale_pressure(pj)
               pj = pj + dpodx*dx
            enddo
         else
            do i = domlo(1), domhi(1)+1
               if (i <= hi(1)+1 .and. i >= lo(1)) &
                  p_g(i,:,:) = scale_pressure(pj)
               pj = pj - dpodx*dx
            enddo
         endif

      else if (abs(gravity(2)) > epsilon(0.0d0)) then

         if (is_undefined(ro_g0)) then
            dpody = -gravity(2)*eosg(mw_avg,pj,295.15d0)
         else
            dpody = -gravity(2)*ro_g0
         endif

         if (gravity(2) <= 0.0d0) then
            do j = domhi(2)+1, domlo(2), -1
               if (j <= hi(2)+1 .and. j>= lo(2)) &
                  p_g(:,j,:) = scale_pressure(pj)
               pj = pj + dpody*dy
            enddo
         else
            do j = domlo(2),domhi(2)+1
               if (j <= hi(2)+1 .and. j >= lo(2)) &
                  p_g(:,j,:) = scale_pressure(pj)
               pj = pj - dpody*dy
            enddo
         endif

      else if (abs(gravity(3)) > epsilon(0.0d0)) then

         if (is_undefined(ro_g0)) then
            dpodz = -gravity(3)*eosg(mw_avg,pj,295.15d0)
         else
            dpodz = -gravity(3)*ro_g0
         endif

         if(gravity(3) <= 0.0d0) then
            do k = domhi(3)+1, domlo(3), -1
               if (k <= hi(3)+1 .and. k >= lo(3)) &
                  p_g(:,:,k) = scale_pressure(pj)
               pj = pj + dpodz*dz
            enddo
         else
            do k = domlo(3),domhi(3)+1
               if (k <= hi(3)+1 .and. k >= lo(3)) &
                  p_g(:,:,k) = scale_pressure(pj)
               pj = pj - dpodz*dz
            enddo
         endif
      endif

  100 continue

      return

 1000 FORMAT(/1X,70('*')//' From: SET_FLUIDBED_P'/' Message: Outflow ',&
         'pressure boundary condition (P_OUTFLOW) not found.',/&
         'All the initial pressures (IC_P_g) or at least one P_OUTFLOW',/&
         'condition need to be specified',/1X,70('*')/)

   end subroutine set_p_g
end module init_fluid_module
