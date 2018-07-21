!
!
!  This module contains the subroutines to perform some of the steps of the
!  projection method.
!
!
module projection_mod

   use amrex_fort_module, only: ar => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one

   implicit none
   private

   ! Define here the unit vectors
   ! This is used to shift index  based on how the variable is staggered
   ! Check e_x, e_y and e_z in mfix_level.H
   integer(c_int), parameter :: e_i(3,3) = reshape ( [1,0,0,0,1,0,0,0,1], [3,3] )

contains

   !
   ! Compute new dt by using the formula derived in
   ! "A Boundary Condition Capturing Method for Multiphase Incompressible Flow"
   ! by Kang et al. (JCP).
   !
   !  dt/2 * ( C+V + sqrt( (C+V)**2 + 4Fx/dx + 4Fy/dy + 4Fz/dz )
   !
   ! where
   !
   ! C = max(|U|)/dx + max(|V|)/dy + max(|W|)/dz    --> Convection
   !
   ! V = 2 * max(mu/ro) * (1/dx^2 + 1/dy^2 +1/dz^2) --> Diffusion
   !
   ! Fx, Fy, Fz = net acceleration due to external forces
   !
   ! WARNING: We use a slightly modified version of C in the implementation below
   !
   subroutine compute_new_dt ( umax, vmax, wmax, romin, mumax, gradp0max, &
                               dx, cfl, steady_state, time, stop_time, dt ) &
        & bind(C)

      use constant, only: gravity

      integer(c_int), intent(in   ) :: steady_state
      real(ar),       intent(in   ) :: umax, vmax, wmax
      real(ar),       intent(in   ) :: mumax, romin
      real(ar),       intent(in   ) :: dx(3), cfl
      real(ar),       intent(in   ) :: gradp0max(3)
      real(ar),       intent(in   ) :: time, stop_time
      real(ar),       intent(inout) :: dt
      real(ar)                      :: old_dt
      real(ar)                      :: c_cfl, v_cfl, f_cfl
      real(ar)                      :: odx, ody, odz
      real(ar)                      :: tmp
      real(ar),       parameter     :: two = 2.0_ar, four = two*two
      real(ar),       parameter     :: eps = epsilon (zero)

      odx    = one / dx(1)
      ody    = one / dx(2)
      odz    = one / dx(3)
      c_cfl  = zero
      v_cfl  = zero
      f_cfl  = zero
      old_dt = dt

      ! Convection
      c_cfl = max ( umax*odx, vmax*ody, wmax*odz )

      ! Viscous
      v_cfl = two * ( mumax / romin ) * ( odx**2 + ody**2 + odz**2 )

      ! Gravity and/or gradient of p0
      f_cfl = abs(gravity(1)-gradp0max(1)) * odx + &
              abs(gravity(2)-gradp0max(2)) * ody + &
              abs(gravity(3)-gradp0max(3)) * odz

      ! Put all together
      tmp = (c_cfl + v_cfl)  + sqrt ( (c_cfl + v_cfl)**2 + four * f_cfl )
      dt  = cfl * two / tmp

      ! Protect against tmp very small
      ! This may happen, for example, when the initial velocity field
      ! is zero for an inviscid flow with no external forcing
      if ( tmp <= eps ) then
         dt = .5 * old_dt
      end if

      ! Don't let the timestep grow by more than 1% per step.
      dt = min ( dt, 1.01*old_dt )

      ! Don't overshoot the final time if not running to steady state
      if (steady_state .eq. 0 .and. stop_time .ge. 0.) then
        if (time+dt .gt. stop_time) &
           dt = stop_time - time
      end if

   end subroutine compute_new_dt

   !
   ! Computes  vel = vel + c * (1/rho) grad(phi)
   !
   ! vel  = velocity            defined at cell centers
   ! ro_g = density field       defined at cell centers
   ! phi  = pressure correction defined at cell centers
   ! c    = real constant
   !
   subroutine add_grad_phicc ( lo, hi, vel, ulo, uhi, ro_g, slo, shi, &
        & phi, dx, c ) bind (C)

      ! Loop bounds
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array bounds
      integer(c_int),  intent(in   ) :: ulo(3), uhi(3)
      integer(c_int),  intent(in   ) :: slo(3), shi(3)

      ! Grid and time spacing
      real(ar),        intent(in   ) :: c, dx(3)

      ! Arrays
      real(ar),        intent(in   ) ::                       &
            ro_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),  &
             phi(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),        intent(inout) ::                       &
             vel(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)

      ! Local variables
      integer(c_int)                 :: i, j, k
      real(ar)                       :: odx, ody, odz
      real(ar)                       :: oro_x_lo, dp_x_lo
      real(ar)                       :: oro_x_hi, dp_x_hi
      real(ar)                       :: oro_y_lo, dp_y_lo
      real(ar)                       :: oro_y_hi, dp_y_hi
      real(ar)                       :: oro_z_lo, dp_z_lo
      real(ar)                       :: oro_z_hi, dp_z_hi

      odx = one / dx(1)
      ody = one / dx(2)
      odz = one / dx(3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               oro_x_lo  = half * ( one/ro_g(i,j,k) + one/ro_g(i-1,j,k) )
                dp_x_lo  =  phi(i,j,k) -  phi(i-1,j,k)

               oro_x_hi  = half * ( one/ro_g(i+1,j,k) + one/ro_g(i,j,k) )
                dp_x_hi  =  phi(i+1,j,k) -  phi(i,j,k)

               vel(i,j,k,1) = vel(i,j,k,1) + half * c * odx * (     &
                   oro_x_hi * dp_x_hi + oro_x_lo * dp_x_lo )

               oro_y_lo  = half * ( one/ro_g(i,j,k) + one/ro_g(i,j-1,k) )
                dp_y_lo  =  phi(i,j,k) -  phi(i,j-1,k)

               oro_y_hi  = half * ( one/ro_g(i,j+1,k) + one/ro_g(i,j,k) )
                dp_y_hi  =  phi(i,j+1,k) -  phi(i,j,k)

               vel(i,j,k,2) = vel(i,j,k,2) + half * c * ody * (     &
                   oro_y_hi * dp_y_hi + oro_y_lo * dp_y_lo )

               oro_z_lo  = half * ( one/ro_g(i,j,k) + one/ro_g(i,j,k-1) )
                dp_z_lo  =  phi(i,j,k) -  phi(i,j,k-1)

               oro_z_hi  = half * ( one/ro_g(i,j,k+1) + one/ro_g(i,j,k) )
                dp_z_hi  =  phi(i,j,k+1) -  phi(i,j,k)

               vel(i,j,k,3) = vel(i,j,k,3) + half * c * odz * (     &
                   oro_z_hi * dp_z_hi + oro_z_lo * dp_z_lo )

            end do
         end do
      end do

   end subroutine add_grad_phicc

   !
   ! Computes  vel = vel + c * (1/rho) grad(phi)
   !
   ! vel  = velocity            defined at cell centers
   ! ro_g = density field       defined at cell centers
   ! phi  = pressure correction defined at nodes
   ! c    = real constant
   !
   subroutine add_grad_phind ( lo, hi, vel, ulo, uhi, ro_g, slo, shi, &
                               phi, rlo, rhi, dx, c ) bind (C)

      ! Loop bounds
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array bounds
      integer(c_int),  intent(in   ) :: ulo(3), uhi(3)
      integer(c_int),  intent(in   ) :: slo(3), shi(3)
      integer(c_int),  intent(in   ) :: rlo(3), rhi(3)

      ! Grid and time spacing
      real(ar),        intent(in   ) :: c, dx(3)

      ! Arrays
      real(ar),        intent(in   ) ::                       &
            ro_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),  &
              phi(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))

      real(ar),        intent(inout) ::                       &
             vel(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)

      ! Local variables
      integer(c_int)                 :: i, j, k
      real(ar)                       :: odx, ody, odz, oro
      real(ar)                       :: phix, phiy, phiz

      odx = one / dx(1)
      ody = one / dx(2)
      odz = one / dx(3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               phix = 0.25d0 * ( &
                  phi(i+1,j,k) +  phi(i+1,j+1,k) +  phi(i+1,j,k+1) +  phi(i+1,j+1,k+1) &
               -  phi(i  ,j,k) -  phi(i  ,j+1,k) -  phi(i  ,j,k+1) -  phi(i  ,j+1,k+1) )

               phiy = 0.25d0 * ( &
                  phi(i,j+1,k) +  phi(i+1,j+1,k) +  phi(i,j+1,k+1) +  phi(i+1,j+1,k+1) &
               -  phi(i,j  ,k) -  phi(i+1,j  ,k) -  phi(i,j  ,k+1) -  phi(i+1,j  ,k+1) )

               phiz = 0.25d0 * ( &
                  phi(i,j,k+1) +  phi(i+1,j,k+1) +  phi(i,j+1,k+1) +  phi(i+1,j+1,k+1) &
               -  phi(i,j,k  ) -  phi(i+1,j,k  ) -  phi(i,j+1,k  ) -  phi(i+1,j+1,k  ) )

               oro = 1.d0/ro_g(i,j,k)

               vel(i,j,k,1) = vel(i,j,k,1) + c * odx * oro * phix
               vel(i,j,k,2) = vel(i,j,k,2) + c * ody * oro * phiy
               vel(i,j,k,3) = vel(i,j,k,3) + c * odz * oro * phiz

            end do
         end do
      end do

   end subroutine add_grad_phind

   !
   ! Add forcing (acceleration) terms to velocity
   ! These terms include the volumetric forces and the explicit part of the
   ! particle/fluid momentum exchange
   !
   subroutine add_forcing ( lo, hi, vel, ulo, uhi, drag, dlo, dhi, &
        & ro_g, slo, shi, rop_g, domlo, domhi, dx, dt )  bind(C)

      use constant, only: gravity

      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3)
      integer(c_int), intent(in   ) :: dlo(3), dhi(3)

      ! Time step width
      real(ar),       intent(in   ) :: dt

      ! Domain bounds
      integer(c_int), intent(in   ) :: domlo(3), domhi(3)

      ! Grid
      real(ar),       intent(in   ) :: dx(3)

      ! Arrays
      real(ar),       intent(in   ) :: &
           ro_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),  &
           rop_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           drag(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),3)

      real(ar),       intent(inout) :: &
           vel(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)

      ! Local variables
      integer(c_int)                :: i, j, k , n

      do n = 1, 3
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               vel(i,j,k,n) = vel(i,j,k,n) + dt * ( gravity(n) + drag(i,j,k,n) / rop_g(i,j,k) )

            end do
         end do
      end do
      end do

   end subroutine add_forcing

   !
   ! This part takes care of performing an implicit solve for the
   ! intermediate velocity.
   ! Currently, this is equivalent to dividing by a diagonal coefficient
   ! since the only implicit term is the fluid/particle momentum exchange.
   !
   ! Upon entry, vel contains the rhs of the system to be solved and on exit
   ! the solution of the system itself. This means that we are solving;
   !
   !    A*vel = rhs  with rhs = vel upon entry
   !
   ! So far the above system reduces to:
   !
   !    vel = vel / (A)_diagonal
   !
   subroutine compute_intermediate_velocity ( lo, hi, vel, ulo, uhi, &
        & f_gds, flo, fhi, rop, slo, shi, dt ) bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3)
      integer(c_int), intent(in   ) :: flo(3), fhi(3)

      ! Time step width
      real(ar),       intent(in   ) :: dt

            ! Arrays
      real(ar),       intent(in   ) :: &
             rop(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           f_gds(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

      real(ar),       intent(inout) :: &
           vel(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)

      ! Local variables
      integer(c_int)                :: i, j, k
      real(ar)                      :: orop

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               orop       = one / rop(i,j,k)

               vel(i,j,k,1) = vel(i,j,k,1) / (one + dt * f_gds(i,j,k) * orop)
               vel(i,j,k,2) = vel(i,j,k,2) / (one + dt * f_gds(i,j,k) * orop)
               vel(i,j,k,3) = vel(i,j,k,3) / (one + dt * f_gds(i,j,k) * orop)

            end do
         end do
      end do

   end subroutine compute_intermediate_velocity

   !
   ! Compute the coefficients of the PPE, i.e. ep_g / ro_g,
   ! at the faces of the pressure cells along the "dir"-axis.
   !
   subroutine compute_bcoeff_cc ( lo, hi, bcoeff, blo, bhi, &
        ro_g, slo, shi, ep_g, dir )  bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: blo(3),bhi(3)

      ! Direction
      integer(c_int), intent(in   ) :: dir

      ! Arrays
      real(ar),       intent(in   ) :: &
           ro_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),       intent(  out) :: &
           bcoeff(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3))

      integer      :: i, j, k, i0, j0, k0

      i0 = e_i(dir,1)
      j0 = e_i(dir,2)
      k0 = e_i(dir,3)

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               bcoeff(i,j,k) = half * (     ep_g(i,j,k) +     ep_g(i-i0,j-j0,k-k0) ) * &
                             & half * ( one/ro_g(i,j,k) + one/ro_g(i-i0,j-j0,k-k0) )
            end do
         end do
      end do

   end subroutine compute_bcoeff_cc

   subroutine compute_bcoeff_nd ( lo, hi, bcoeff, blo, bhi, &
        ro_g, slo, shi, ep_g, dir )  bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: blo(3),bhi(3)

      ! Direction
      integer(c_int), intent(in   ) :: dir

      ! Arrays
      real(ar),       intent(in   ) :: &
           ro_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),       intent(  out) :: &
           bcoeff(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3))

      integer      :: i, j, k

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               bcoeff(i,j,k) =  ep_g(i,j,k) / ro_g(i,j,k)
            end do
         end do
      end do

   end subroutine compute_bcoeff_nd

   !
   ! Set the boundary condition for Pressure Poisson Equation (PPE)
   !
   ! MLMG expects the BC type to be the uniform on each domain wall.
   ! Since mfix allows for BC patches on each wall, we first check that
   ! the user-provided BCs are uniform, and then return a single BC type for
   ! each domain wall.
   !
   subroutine set_ppe_bc ( bc_lo, bc_hi, domlo, domhi, ng, bct_ilo, bct_ihi, &
        & bct_jlo, bct_jhi, bct_klo, bct_khi, singular )  bind(C)

      use amrex_lo_bctypes_module
      use bc

      ! Array of global BC types
      integer(c_int), intent(  out) :: bc_lo(3), bc_hi(3)

      ! Domain bounds
      integer(c_int), intent(in   ) :: domlo(3), domhi(3), ng

      ! Whether the system is singular or not
      integer(c_int), intent(  out) :: singular

      ! Arrays of point-by-point BC types
      integer(c_int), intent(in   )  ::                                 &
           & bct_ilo(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bct_ihi(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bct_jlo(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bct_jhi(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bct_klo(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
           & bct_khi(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

      ! Local variables
      integer(c_int)                :: bc_face

      !
      ! By default, all the BCs are Neumann
      !
      singular = 1
      bc_lo    = amrex_lo_neumann
      bc_hi    = amrex_lo_neumann

      !
      ! BC -- X direction
      !
      if ( cyclic_x ) then
         bc_lo(1) = amrex_lo_periodic
         bc_hi(1) = amrex_lo_periodic
      else

         ! X at domlo(1)
         bc_face = get_bc_face(bct_ilo)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) ) then
            bc_lo(1) = amrex_lo_dirichlet
         end if

         ! X at domhi(1)
         bc_face = get_bc_face(bct_ihi)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) ) then
            bc_hi(1) = amrex_lo_dirichlet
         end if

      end if


      !
      ! BC -- Y direction
      !
      if ( cyclic_y ) then
         bc_lo(2) = amrex_lo_periodic
         bc_hi(2) = amrex_lo_periodic
      else

         ! Y at domlo(2)
         bc_face = get_bc_face(bct_jlo)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) ) then
            bc_lo(2) = amrex_lo_dirichlet
         end if

         ! Y at domhi(2)
         bc_face = get_bc_face(bct_jhi)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) ) then
            bc_hi(2) = amrex_lo_dirichlet
         end if

      end if

      !
      ! BC -- Z direction
      !
      if ( cyclic_z ) then
         bc_lo(3) = amrex_lo_periodic
         bc_hi(3) = amrex_lo_periodic
      else

         ! Z at domlo(3)
         bc_face = get_bc_face(bct_klo)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) ) then
            bc_lo(3) = amrex_lo_dirichlet
         end if

         ! Z at domhi(3)
         bc_face = get_bc_face(bct_khi)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) ) then
            bc_hi(3) = amrex_lo_dirichlet
         end if

      end if

      !
      ! Check whether the system is non-singular
      !
      if ( any ( bc_hi == amrex_lo_dirichlet ) .or. &
           any ( bc_lo == amrex_lo_dirichlet ) )   singular = 0

      contains

         !
         ! Test whether the BC type is the same everywhere on
         ! the face. If BC is uniform on face, it returns its value
         !
         function get_bc_face (bct_array) result (bc_face)
            integer(c_int), intent(in   ) :: bct_array(:,:,:)
            integer                       :: bc_face
            integer                       :: is, ie, js, je

            ! Do not consider the edges: they may cause problems
            is = 3
            ie = size (bct_array,1) - 2
            js = 3
            je = size (bct_array,2) - 2

            bc_face = bct_array(is,js,1)

            if ( .not. all (bct_array(is:ie,js:je,1) == bc_face) ) then
               stop "BC type must be uniform on each face of the domain"
            end if

         end function get_bc_face

   end subroutine set_ppe_bc

   !
   ! Compute the cell-centered divergence of ep_g * u_g
   !
   subroutine compute_diveucc ( lo, hi, diveu, slo, shi, ep_g, vel, ulo, uhi, dx) &
      bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Arrays bounds
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3)

      ! Grid
      real(ar),       intent(in   ) :: dx(3)

      ! Array
      real(ar),       intent(  out) :: &
           diveu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),       intent(in   ) :: &
           vel(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3), &
          ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Local variables
      integer  :: i, j, k
      real(ar) :: odx, ody, odz
      real(ar) :: eu_e, eu_w, ev_n, ev_s, ew_t, ew_b

      odx = one / dx(1)
      ody = one / dx(2)
      odz = one / dx(3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! Face values
               eu_e = half * ( ep_g(i+1,j,k  ) + ep_g(i,j,k  ) ) * &
                      half * (  vel(i+1,j,k,1) +  vel(i,j,k,1) )
               eu_w = half * ( ep_g(i-1,j,k  ) + ep_g(i,j,k  ) ) * &
                      half * (  vel(i-1,j,k,1) +  vel(i,j,k,1) )

               ev_n = half * ( ep_g(i,j+1,k  ) + ep_g(i,j,k  ) ) * &
                      half * (  vel(i,j+1,k,2) +  vel(i,j,k,2) )
               ev_s = half * ( ep_g(i,j-1,k  ) + ep_g(i,j,k  ) ) * &
                      half * (  vel(i,j-1,k,2) +  vel(i,j,k,2) )

               ew_t = half * ( ep_g(i,j,k+1  ) + ep_g(i,j,k  ) ) * &
                      half * (  vel(i,j,k+1,3) +  vel(i,j,k,3) )
               ew_b = half * ( ep_g(i,j,k-1  ) + ep_g(i,j,k  ) ) * &
                      half * (  vel(i,j,k-1,3) +  vel(i,j,k,3) )

               ! Divergence
               diveu(i,j,k) = (eu_e - eu_w) * odx + (ev_n - ev_s) * ody + &
                    &         (ew_t - ew_b) * odz

            end do
         end do
      end do

   end subroutine compute_diveucc

   subroutine compute_diveund ( lo, hi, diveu, slo, shi, vec, ulo, uhi, dx) &
      bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Arrays bounds
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3)

      ! Grid
      real(ar),       intent(in   ) :: dx(3)

      ! Array
      real(ar),       intent(  out) :: &
           diveu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),       intent(in   ) :: &
           vec(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)

      ! Local variables
      integer  :: i, j, k
      real(ar) :: odx, ody, odz
      real(ar) :: eu_x, eu_y, eu_z

      odx = one / dx(1)
      ody = one / dx(2)
      odz = one / dx(3)

      ! Note these are nodal indices
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! Divergence
               eu_x = ( vec(i  ,j  ,k  ,1) + vec(i  ,j-1,k  ,1) &
                       +vec(i  ,j  ,k-1,1) + vec(i  ,j-1,k-1,1) &
                       -vec(i-1,j  ,k  ,1) - vec(i-1,j-1,k  ,1) &
                       -vec(i-1,j  ,k-1,1) - vec(i-1,j-1,k-1,1) )

               eu_y = ( vec(i  ,j  ,k  ,2) + vec(i-1,j  ,k  ,2) &
                       +vec(i  ,j  ,k-1,2) + vec(i-1,j  ,k-1,2) &
                       -vec(i  ,j-1,k  ,2) - vec(i-1,j-1,k  ,2) &
                       -vec(i  ,j-1,k-1,2) - vec(i-1,j-1,k-1,2) )

               eu_z = ( vec(i  ,j  ,k  ,3) + vec(i-1,j  ,k  ,3) &
                       +vec(i  ,j-1,k  ,3) + vec(i-1,j-1,k  ,3) &
                       -vec(i  ,j  ,k-1,3) - vec(i-1,j  ,k-1,3) &
                       -vec(i  ,j-1,k-1,3) - vec(i-1,j-1,k-1,3) )

                diveu(i,j,k) = 0.25d0 * (eu_x*odx + eu_y*ody + eu_z*odz)

!               if (i.eq.0) print *,'DIVEU ',j,k,diveu(i,j,k), eu_x, eu_y, eu_z

            end do
         end do
      end do

   end subroutine compute_diveund

   subroutine compute_gradp0_max ( lo, hi, p0, slo, shi, gp0_max, dx, nodal_pressure) &
      bind (C)

      ! Loop bounds
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array bounds
      integer(c_int),  intent(in   ) :: slo(3), shi(3)

      ! Grid and time spacing
      real(ar),        intent(in   ) :: dx(3)
      real(ar),        intent(  out) :: gp0_max(3)

      ! Arrays
      real(ar),        intent(in   ) ::                       &
              p0(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer(c_int),  intent(in   ) :: nodal_pressure

      ! Local variables
      integer(c_int) :: i, j, k
      real(ar)       :: odx, ody, odz

      odx = one / dx(1)
      ody = one / dx(2)
      odz = one / dx(3)

      gp0_max(:) = zero

      if ( nodal_pressure == 1 ) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)

                  gp0_max(1) = max( gp0_max(1), abs( &
                       p0(i+1,j,k) + p0(i+1,j+1,k) + p0(i+1,j,k+1) + p0(i+1,j+1,k+1) &
                       - p0(i  ,j,k) - p0(i  ,j+1,k) - p0(i  ,j,k+1) - p0(i  ,j+1,k+1) ) )
                  
                  gp0_max(2) = max( gp0_max(2), abs( &
                       p0(i,j+1,k) + p0(i+1,j+1,k) + p0(i,j+1,k+1) + p0(i+1,j+1,k+1) &
                       - p0(i,j  ,k) - p0(i+1,j  ,k) - p0(i,j  ,k+1) - p0(i+1,j  ,k+1) ) )
                  
                  gp0_max(3) = max( gp0_max(3), abs( &
                       p0(i,j,k+1) + p0(i+1,j,k+1) + p0(i,j+1,k+1) + p0(i+1,j+1,k+1) &
                       - p0(i,j,k  ) - p0(i+1,j,k  ) - p0(i,j+1,k  ) - p0(i+1,j+1,k  ) ) )
                  

               end do
            end do
         end do
         
         gp0_max(1) = gp0_max(1) * odx * 0.25d0
         gp0_max(2) = gp0_max(2) * ody * 0.25d0
         gp0_max(3) = gp0_max(3) * odz * 0.25d0

      else

         ! The box is cell-centered
         do k = lo(3), hi(3)+1
            do j = lo(2), hi(2)+1
               do i = lo(1), hi(1)+1
                  
                  gp0_max(1) = max( gp0_max(1), abs(p0(i,j,k) - p0(i-1,j,k)))
                  gp0_max(2) = max( gp0_max(2), abs(p0(i,j,k) - p0(i,j-1,k)))
                  gp0_max(3) = max( gp0_max(3), abs(p0(i,j,k) - p0(i,j,k-1)))
                  
               end do
            end do
         end do
         
         gp0_max(1) = gp0_max(1) * odx
         gp0_max(2) = gp0_max(2) * ody
         gp0_max(3) = gp0_max(3) * odz
         
      end if

   end subroutine compute_gradp0_max

end module projection_mod
