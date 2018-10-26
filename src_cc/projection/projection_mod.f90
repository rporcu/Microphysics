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
   ! Add gravitational acceleration terms to velocity
   !
   subroutine add_gravity ( lo, hi, vel, ulo, uhi, dt )  bind(C) 

      use constant, only: gravity

      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: ulo(3), uhi(3)

      ! Time step width
      real(ar),       intent(in   ) :: dt

      real(ar),       intent(inout) :: &
           vel(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)

      ! Local variables
      integer(c_int)                :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               vel(i,j,k,1) = vel(i,j,k,1) + dt * gravity(1) 
               vel(i,j,k,2) = vel(i,j,k,2) + dt * gravity(2)
               vel(i,j,k,3) = vel(i,j,k,3) + dt * gravity(3)

            end do
         end do
      end do

   end subroutine add_gravity

   !
   ! This adds both components of the drag term
   ! Here f_gds = beta
   !      drag  = beta * particle_velocity
   ! 
   ! So the drag term we add is beta * (particle_velocity - fluid_velocity)
   !                          = drag - f_gds * fluid_velocity
   !
   subroutine add_drag_terms ( lo, hi, vel, ulo, uhi, &
                               f_gds  , flo, fhi, &
                               drag   , dlo, dhi, &
                               rop    , slo, shi, &
                               volfrac, vlo, vhi, dt ) bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: ulo(3), uhi(3)
      integer(c_int), intent(in   ) :: dlo(3), dhi(3)
      integer(c_int), intent(in   ) :: flo(3), fhi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: vlo(3), vhi(3)

      ! Time step width
      real(ar),       intent(in   ) :: dt

      ! Arrays
      real(ar),       intent(in   ) :: &
               rop(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
             f_gds(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3)), &
              drag(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),3), &
           volfrac(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

      real(ar),       intent(inout) :: &
           vel(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)

      ! Local variables
      integer(c_int)                :: i, j, k
      real(ar)                      :: orop

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               orop       = dt * volfrac(i,j,k) / rop(i,j,k)

               vel(i,j,k,1) = ( vel(i,j,k,1) + drag(i,j,k,1) * orop) / (one + f_gds(i,j,k) * orop)
               vel(i,j,k,2) = ( vel(i,j,k,2) + drag(i,j,k,2) * orop) / (one + f_gds(i,j,k) * orop)
               vel(i,j,k,3) = ( vel(i,j,k,3) + drag(i,j,k,3) * orop) / (one + f_gds(i,j,k) * orop)

            end do
         end do
      end do

   end subroutine add_drag_terms

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
        & bct_jlo, bct_jhi, bct_klo, bct_khi)  bind(C)

      use amrex_lo_bctypes_module
      use bc

      ! Array of global BC types
      integer(c_int), intent(  out) :: bc_lo(3), bc_hi(3)

      ! Domain bounds
      integer(c_int), intent(in   ) :: domlo(3), domhi(3), ng

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
         bc_face = get_bc_face(bct_ilo,ng)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) ) then
            bc_lo(1) = amrex_lo_dirichlet
         end if

         ! X at domhi(1)
         bc_face = get_bc_face(bct_ihi,ng)
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
         bc_face = get_bc_face(bct_jlo,ng)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) ) then
            bc_lo(2) = amrex_lo_dirichlet
         end if

         ! Y at domhi(2)
         bc_face = get_bc_face(bct_jhi,ng)
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
         bc_face = get_bc_face(bct_klo,ng)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) ) then
            bc_lo(3) = amrex_lo_dirichlet
         end if

         ! Z at domhi(3)
         bc_face = get_bc_face(bct_khi,ng)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) ) then
            bc_hi(3) = amrex_lo_dirichlet
         end if

      end if

   contains

      !
      ! Test whether the BC type is the same everywhere on
      ! the face. If BC is uniform on face, it returns its value
      !
      function get_bc_face (bct_array,nghost) result (bc_face)
         integer(c_int), intent(in   ) :: bct_array(:,:,:)
         integer(c_int), intent(in   ) :: nghost
         integer                       :: bc_face
         integer                       :: is, ie, js, je

         ! Do not consider the edges: they may cause problems
         is = nghost+1
         ie = size(bct_array,1) - nghost
         js = nghost+1
         je = size(bct_array,2) - nghost

         bc_face = bct_array(is,js,1)

         if ( .not. all (bct_array(is:ie,js:je,1) == bc_face) ) then
            stop "BC type must be uniform on each face of the domain"
         end if

      end function get_bc_face

   end subroutine set_ppe_bc

end module projection_mod
