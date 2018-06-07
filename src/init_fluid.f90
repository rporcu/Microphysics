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

      use amrex_fort_module, only: rt => amrex_real
      use iso_c_binding    , only: c_int

      use calc_mu_g_module, only: calc_mu_g
      use fld_const       , only: ro_g0

      implicit none

! Dummy arguments .....................................................//
      integer(c_int), intent(in   ) ::  lo(3),  hi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)

      real(rt), intent(inout) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(inout) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(inout) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(inout) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(rt), intent(inout) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(rt), intent(inout) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(rt), intent(inout) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(rt), intent(inout) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(inout) :: lambda_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(rt), intent(in   ) :: dx, dy, dz
      real(rt), intent(in   ) :: xlength, ylength, zlength

      integer i

      ! Set user specified initial conditions (IC)
      call set_ic(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, &
                  domlo, domhi, dx, dy, dz, u_g, v_g, w_g)

      ro_g  = ro_g0
      rop_g = ro_g0 * ep_g

      call calc_mu_g(slo, shi, lo, hi, mu_g, lambda_g)

   end subroutine init_fluid

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: init_fluid_restart                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine init_fluid_restart(slo, shi, lo, hi, mu_g, lambda_g) &
      bind(C, name="init_fluid_restart")

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int
      use calc_mu_g_module, only: calc_mu_g

      implicit none

! Dummy arguments .....................................................//
      integer(c_int), intent(in   ) ::  lo(3),  hi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)

      real(rt), intent(inout) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(inout) :: lambda_g&
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
                     domlo, domhi, dx, dy, dz, u_g, v_g, w_g)

      use ic, only: dim_ic, ic_defined
      use ic, only: ic_p_g, ic_u_g, ic_v_g, ic_w_g
      use ic, only: ic_x_e, ic_y_n, ic_z_t
      use ic, only: ic_x_w, ic_y_s, ic_z_b
      use bc, only: delp_x, delp_y, delp_z
      use scales, only: scale_pressure
      use param, only: undefined, is_defined, zero

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

      use calc_cell_module, only: calc_cell_ic

      implicit none

      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)
      real(rt), intent(in   ) :: dx, dy, dz

      real(rt), intent(inout) ::  u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(rt), intent(inout) ::  v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(rt), intent(inout) ::  w_g&
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
      real(rt) :: pgx, ugx, vgx, wgx, pval

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
               if (ulo(1).lt.domlo(1) .and. domlo(1) == istart) &
                  u_g(ulo(1):istart-1,jstart:jend,kstart:kend) = ugx
               if (uhi(1).gt.domhi(1) .and. domhi(1) == iend  ) &
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
               if (vlo(2).lt.domlo(2) .and. domlo(2) == jstart) &
                  v_g(istart:iend,vlo(2):jstart-1,kstart:kend) = vgx
               if (vhi(2).gt.domhi(2) .and. domhi(2) == jend  ) &
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
               if (wlo(3).lt.domlo(3) .and. domlo(3) == kstart) &
                  w_g(istart:iend,jstart:jend,wlo(3):kstart-1) = wgx
               if (whi(3).gt.domhi(3) .and. domhi(3) == kend  ) &
                  w_g(istart:iend,jstart:jend,kend+1:whi(3)  ) = wgx
            end if

         endif
      enddo

   end subroutine set_ic

end module init_fluid_module
