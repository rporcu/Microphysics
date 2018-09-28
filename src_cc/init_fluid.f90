module init_fluid_module
   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: init_fluid                                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine init_fluid(slo, shi, lo, hi, &
                         domlo, domhi, ep_g, ro_g, rop_g, p_g, vel_g, &
                         mu_g, lambda_g, dx, dy, dz, xlength, ylength, zlength) &
      bind(C, name="init_fluid")

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

      use fld_const       , only: ro_g0
      use calc_mu_g_module, only: calc_mu_g

      implicit none

! Dummy arguments .....................................................//
      integer(c_int), intent(in   ) ::  lo(3),  hi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)

      real(rt), intent(inout) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(inout) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(inout) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(inout) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(rt), intent(inout) :: vel_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

      real(rt), intent(inout) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(inout) :: lambda_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(rt), intent(in   ) :: dx, dy, dz
      real(rt), intent(in   ) :: xlength, ylength, zlength

      ! Set user specified initial conditions (IC)
      call set_ic(slo, shi, domlo, domhi, dx, dy, dz, vel_g)

      ! call init_periodic_vortices ( lo, hi, vel_g, slo, shi, dx, dy, dz, domlo)

      ! Set the initial fluid density and viscosity
      ro_g  = ro_g0
      rop_g = ro_g0 * ep_g

      call calc_mu_g(slo, shi, lo, hi, mu_g, lambda_g)

   end subroutine init_fluid

   subroutine init_periodic_vortices ( lo, hi, vel, slo, shi, dx, dy, dz, domlo)

      use amrex_fort_module, only: ar => amrex_real
      use iso_c_binding ,    only: c_int
      use param,             only: zero, half, one
 
      implicit none

      ! Array bounds
      integer(c_int),   intent(in   ) :: slo(3), shi(3)

      ! Tile bounds
      integer(c_int),   intent(in   ) ::  lo(3),  hi(3)
      
      ! Grid and domain lower bound
      integer(c_int),   intent(in   ) :: domlo(3)
      real(ar),         intent(in   ) :: dx, dy, dz
      
      ! Arrays
      real(ar),         intent(inout) ::                   &
           & vel(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

      ! Local variables
      integer(c_int)                  :: i, j, k, plane
      real(ar)                        :: x, y, z
      real(ar)                        :: twopi = 8.0_ar * atan(one) 

      plane = 1

      select case ( plane )

      case (1)  ! x-y plane
         
         ! x-direction
         do j = lo(2), hi(2)
            y =  ( real(j,ar) + half ) * dy
            do i = lo(1), hi(1)
               x =  ( real(i,ar) + half ) * dx
               do k = lo(3), hi(3)
                  vel(i,j,k,1) = tanh ( 30.0_ar * (0.25_ar - abs ( y - 0.5_ar ) ) )
                  vel(i,j,k,2) = 0.05_ar * sin ( twopi * x )
                  vel(i,j,k,3) = zero
               end do
            end do
         end do
         
      case (2)  ! x-z plane
         
         ! x-direction
         do k = lo(3), hi(3)
            z =  ( real(k,ar) + half ) * dz
            do i = lo(1), hi(1)
               x =  ( real(i,ar) + half ) * dx
               do j = lo(2), hi(2)
                  vel(i,j,k,1) = tanh ( 30.0_ar * (0.25_ar - abs ( z - 0.5_ar ) ) )
                  vel(i,j,k,2) = zero
                  vel(i,j,k,3) = 0.05_ar * sin ( twopi * x )
               end do
            end do
         end do
         
      case (3)  ! y-z plane
         
         ! x-direction
         do k = lo(3), hi(3)
            z =  ( real(k,ar) + half ) * dz
            do j = lo(2), hi(2)
               y =  ( real(j,ar) + half ) * dy
               do i = lo(1), hi(1)
                  vel(i,j,k,1) = zero 
                  vel(i,j,k,2) = tanh ( 30.0_ar * (0.25_ar - abs ( z - 0.5_ar ) ) )
                  vel(i,j,k,3) = 0.05_ar * sin ( twopi * y )
               end do
            end do
         end do
         
      end select 
         
   end subroutine init_periodic_vortices

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
   subroutine set_ic(slo, shi, domlo, domhi, dx, dy, dz, vel_g)

      use ic, only: dim_ic, ic_defined
      use ic, only: ic_p_g, ic_u_g, ic_v_g, ic_w_g
      use ic, only: ic_x_e, ic_y_n, ic_z_t
      use ic, only: ic_x_w, ic_y_s, ic_z_b
      use scales, only: scale_pressure
      use param, only: undefined, is_defined

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

      use calc_cell_module, only: calc_cell_ic

      implicit none

      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)
      real(rt), intent(in   ) :: dx, dy, dz

      real(rt), intent(inout) ::  vel_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

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
      real(rt) :: pgx, ugx, vgx, wgx

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
               istart = max(slo(1), i_w)
               jstart = max(slo(2), j_s)
               kstart = max(slo(3), k_b)
               iend   = min(shi(1), i_e)
               jend   = min(shi(2), j_n)
               kend   = min(shi(3), k_t)
               vel_g(istart:iend,jstart:jend,kstart:kend,1) = ugx
               if (slo(1).lt.domlo(1) .and. domlo(1) == istart) &
                  vel_g(slo(1):istart-1,jstart:jend,kstart:kend,1) = ugx
               if (shi(1).gt.domhi(1) .and. domhi(1) == iend  ) &
                  vel_g(iend+1:shi(1)  ,jstart:jend,kstart:kend,1) = ugx
            end if

            if (is_defined(vgx)) then
               istart = max(slo(1), i_w)
               jstart = max(slo(2), j_s)
               kstart = max(slo(3), k_b)
               iend   = min(shi(1), i_e)
               jend   = min(shi(2), j_n)
               kend   = min(shi(3), k_t)
               vel_g(istart:iend,jstart:jend,kstart:kend,2) = vgx
               if (slo(2).lt.domlo(2) .and. domlo(2) == jstart) &
                  vel_g(istart:iend,slo(2):jstart-1,kstart:kend,2) = vgx
               if (shi(2).gt.domhi(2) .and. domhi(2) == jend  ) &
                  vel_g(istart:iend,jend+1:shi(2)  ,kstart:kend,2) = vgx
            end if

            if (is_defined(wgx)) then
               istart = max(slo(1), i_w)
               jstart = max(slo(2), j_s)
               kstart = max(slo(3), k_b)
               iend   = min(shi(1), i_e)
               jend   = min(shi(2), j_n)
               kend   = min(shi(3), k_t)
               vel_g(istart:iend,jstart:jend,kstart:kend,3) = wgx
               if (slo(3).lt.domlo(3) .and. domlo(3) == kstart) &
                  vel_g(istart:iend,jstart:jend,slo(3):kstart-1,3) = wgx
               if (shi(3).gt.domhi(3) .and. domhi(3) == kend  ) &
                  vel_g(istart:iend,jstart:jend,kend+1:shi(3)  ,3) = wgx
            end if

         endif
      enddo

   end subroutine set_ic

end module init_fluid_module
