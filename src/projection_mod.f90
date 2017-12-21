! 
!              
!  This module contains the subroutines to perform some of the steps of the
!  projection method.
!  
!  Author: Michele Rosso
! 
!  Date: October 13, 2017
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
   subroutine compute_new_dt ( umax, vmax, wmax, romin, mumax, dx, cfl, dt ) &
        & bind(C)
      
      ! subroutine compute_new_dt ( umax, vmax, wmax, fgdsumax, fgdsvmax fgdswmax, &
      !      dragumax, dragvmax, dragwmax, mumax, romin, dx, dt ) &
      

      use constant, only: gravity 

      real(ar),       intent(in   ) :: umax, vmax, wmax
      ! real(ar),       intent(in   ) :: fgdsumax, fgdsvmax fgdswmax
      ! real(ar),       intent(in   ) :: dragumax, dragvmax, dragwmax
      real(ar),       intent(in   ) :: mumax, romin
      real(ar),       intent(in   ) :: dx(3), cfl
      real(ar),       intent(inout) :: dt
      real(ar)                      :: old_dt
      real(ar)                      :: c_cfl, v_cfl, f_cfl
      real(ar)                      :: odx, ody, odz
      real(ar)                      :: fp_x, fp_y, fp_z, tmp
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
      ! c_cfl = umax*odx + vmax*ody + wmax*odz  <-- Too restricting
      c_cfl = max ( umax*odx, vmax*ody, wmax*odz )
      
      ! Viscous 
      v_cfl = two * ( mumax / romin ) * ( odx**2 + ody**2 + odz**2 )
              
      ! Gravity
      f_cfl = f_cfl + gravity(1) * odx + gravity(2) * ody &
           &        + gravity(3) * odz
      
      ! Put all together
      tmp = (c_cfl + v_cfl)  + sqrt ( (c_cfl + v_cfl)**2 + four * f_cfl )
      dt  = cfl * two / tmp

      ! Protect against tmp very small
      ! This may happen, for example, when the initial velocity field
      ! is zero for an inviscid flow with no external forcing  
      if ( tmp <= eps ) then
         dt = old_dt
      else 
         dt = min ( dt, old_dt )
      end if
      
   end subroutine compute_new_dt




   !
   ! Computes the term RHS = - div (uu) + div (tau)/rop
   ! along direction "dir"
   !
   ! dir = 1, 2, 3 ( 1=x, 2=y, 3=z ) 
   !
   subroutine compute_fluid_acceleration ( lo, hi, rhs, rlo, rhi, sl, &
        & u, ulo, uhi, v, vlo, vhi, w, wlo, whi, mu, slo, shi, rop,   &
        & dx, dir, order ) bind(C)

      use convection_mod
      use diffusion_mod

      ! Loop bounds
      integer(c_int), intent(in   ) :: lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: rlo(3), rhi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3)
      integer(c_int), intent(in   ) :: vlo(3), vhi(3)
      integer(c_int), intent(in   ) :: wlo(3), whi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)

      ! Scheme order
      integer(c_int), intent(in   ) :: order

      ! Grid 
      real(ar),       intent(in   ) :: dx(3)

      ! Direction
      integer(c_int), intent(in   ) :: dir

      ! Arrays
      real(ar),       intent(in   ) ::                       &
           &   u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           &   v(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
           &   w(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)), &
           &  mu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           & rop(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           &  sl(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))

      real(ar),       intent(inout) ::                       &
           & rhs(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))

      ! Local working arrays
      real(ar)                      ::                        &
           & conv(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3)), &
           & diff(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))

      ! Local variables
      integer                       :: i, j, k
      integer                       :: i0, j0, k0
      real(ar)                      :: orop
      
      ! Compute convection term
      select case ( dir )
      case (1)
         call compute_ugradu_x ( lo, hi, u, ulo, uhi, v, vlo, vhi, &
              & w, wlo, whi, sl, conv, dx, order )  

         call compute_divtau_x ( lo, hi, u, ulo, uhi, v, vlo, vhi, &
              & w, wlo, whi, mu, slo, shi, diff, dx )
      case(2)
         call compute_ugradu_y ( lo, hi, u, ulo, uhi, v, vlo, vhi,  &
              & w, wlo, whi, sl, conv, dx, order )

         call compute_divtau_y ( lo, hi, u, ulo, uhi, v, vlo, vhi, &
              & w, wlo, whi, mu, slo, shi, diff, dx )
      case(3)         
         call compute_ugradu_z ( lo, hi, u, ulo, uhi, v, vlo, vhi,  &
              & w, wlo, whi, sl, conv, dx, order )

         call compute_divtau_z ( lo, hi, u, ulo, uhi, v, vlo, vhi, &
              & w, wlo, whi, mu, slo, shi, diff, dx )
      end select

      !
      !
      !  REMEMBER TO DIVIDE div(tau) by rho !!!!!
      !
      !
      i0 = e_i(dir,1)
      j0 = e_i(dir,2)
      k0 = e_i(dir,3)

      orop = zero

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               orop       = half * ( one/rop(i,j,k) + one/rop(i-i0,j-j0,k-k0) )
               rhs(i,j,k) = -conv(i,j,k) + diff(i,j,k) * orop
            end do
         end do
      end do

   end subroutine compute_fluid_acceleration


   !
   ! Computes  u_i = u_i + C * (beta) * (dphi/dx_i)
   !
   ! u_i  = i-th component of a staggered vector field u.
   !
   ! beta = scalar field defined at u_i location.
   !
   ! phi  = scalar field defined at cell centers
   ! 
   ! C    = real constant
   !
   ! dir  = 1,2,3 indicates x, y, z direction respectively.
   !  
   subroutine add_gradient ( lo, hi, u_i, ulo, uhi, beta,    &
        & phi, slo, shi, dx, c, dir ) bind (C)

      ! Loop bounds
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array bounds
      integer(c_int),  intent(in   ) :: ulo(3), uhi(3)
      integer(c_int),  intent(in   ) :: slo(3), shi(3)

      ! Grid and time spacing
      real(ar),        intent(in   ) :: c, dx(3)

      ! Direction
      integer(c_int),  intent(in   ) :: dir

      ! Arrays
      real(ar),        intent(in   ) ::                        &
           & beta(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),  &
           &  phi(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),        intent(inout) ::                           &
           & u_i(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      ! Local variables
      integer(c_int)                 :: i, j, k, i0, j0, k0
      real(ar)                       :: codx

      i0 = e_i(dir,1)
      j0 = e_i(dir,2)
      k0 = e_i(dir,3)

      codx = c / dx(dir) 
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               u_i(i,j,k) = u_i(i,j,k) + codx * beta(i,j,k) * &
                    &      ( phi(i,j,k) - phi(i-i0,j-j0,k-k0) )


            end do
         end do
      end do

   end subroutine add_gradient



   !
   ! Add forcing (acceleration) terms to velocity component u_i
   ! 
   subroutine add_forcing ( lo, hi, u_i, ulo, uhi, ro_g, slo, shi, &
        & domlo, domhi, dx, dt, dir )  bind(C)
      
      use constant, only: gravity
      use bc      , only: delp_x, delp_y, delp_z 
      use scales,   only: p_scale
      
      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3)

      ! Direction
      integer(c_int), intent(in   ) :: dir

      ! Time step width
      real(ar),       intent(in   ) :: dt

      ! Domain bounds
      integer(c_int), intent(in   ) :: domlo(3), domhi(3)

      ! Grid
      real(ar),       intent(in   ) :: dx(3)
      
      ! Arrays
      real(ar),       intent(in   ) :: &
           ro_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),       intent(inout) :: &
           u_i(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      
      ! Local variables
      integer(c_int)                :: i, j, k 
      real(ar)                      :: odx(3), orog, acc 


      ! 1/dx
      odx = one / dx
      
      select case (dir)
      case(1)                   !X direction

         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  acc = zero

                  ! Pressure drop at boundaries if specified
                  orog = half * ( one/ro_g(i,j,k) + one/ro_g(i-1,j,k) )
                  if ( i == domlo(1) .or. i == domhi(1)+1 ) then
                     acc =  p_scale * delp_x * orog * odx(dir)
                  end if

                  acc = acc + gravity(dir)

                  u_i(i,j,k) = u_i(i,j,k) + dt * acc

               end do
            end do
         end do

      case(2)                   !y direction

         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  acc = zero

                  ! Pressure drop at boundaries if specified
                  orog = half * ( one/ro_g(i,j,k) + one/ro_g(i,j-1,k) )
                  if ( j == domlo(2) .or. j == domhi(2)+1 ) then
                     acc =  p_scale * delp_y * orog * odx(dir)
                  end if

                  acc = acc + gravity(dir)

                  u_i(i,j,k) = u_i(i,j,k) + dt * acc

               end do
            end do
         end do

      case(3)                   !Z direction

         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  acc = zero

                  ! Pressure drop at boundaries if specified
                  orog = half * ( one/ro_g(i,j,k) + one/ro_g(i,j,k-1) )
                  if ( k == domlo(3) .or. k == domhi(3)+1 ) then
                     acc =  p_scale * delp_z * orog * odx(dir)
                  end if

                  acc = acc + gravity(dir)

                  u_i(i,j,k) = u_i(i,j,k) + dt * acc

               end do
            end do
         end do


      case default
         stop "projection_mod: add_forcing: argument dir must be either 1,2, or 3"
      end select



   end subroutine add_forcing


   !
   ! Compute the coefficients of the PPE, i.e. 1 / ro_g = eps_g/rho_g,
   ! at the faces of the pressure cells along the "dir"-axis.
   ! 
   subroutine compute_oro_g ( lo, hi, oro_g, alo, ahi, &
        rop_g, slo, shi, ep_g, dir )  bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: alo(3),ahi(3)

      ! Direction
      integer(c_int), intent(in   ) :: dir

      ! Arrays
      real(ar),       intent(in   ) :: &
           rop_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),       intent(  out) :: &
           oro_g(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

      integer      :: i, j, k, i0, j0, k0

      i0 = e_i(dir,1)
      j0 = e_i(dir,2)
      k0 = e_i(dir,3)


      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               oro_g(i,j,k) = half * ( ep_g(i,j,k) / rop_g(i,j,k) + &
                    & ep_g(i-i0,j-j0,k-k0) / rop_g(i-i0,j-j0,k-k0) )
            end do
         end do
      end do

   end subroutine compute_oro_g


   !
   ! Set the boundary condition for PPE 
   ! 
   subroutine set_ppe_bc ( bc_lo, bc_hi, domlo, domhi, bct_ilo, bct_ihi, &
        & bct_jlo, bct_jhi, bct_klo, bct_khi, singular )  bind(C)
      
      use amrex_lo_bctypes_module
      use bc
      
      ! Array of global BC types
      integer(c_int), intent(  out) :: bc_lo(3), bc_hi(3)

      ! Domain bounds
      integer(c_int), intent(in   ) :: domlo(3), domhi(3)

      ! Whether the system is singular or not
      integer(c_int), intent(  out) :: singular
      
      ! Arrays of point-by-point BC types 
      integer(c_int), intent(in   ), target  ::  &
           & bct_ilo(domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2), &
           & bct_ihi(domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2), &
           & bct_jlo(domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2), &
           & bct_jhi(domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2), &
           & bct_klo(domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2), &
           & bct_khi(domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2)

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

            ! Do not considere the edges: they may cause problems
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
   ! Subroutine to impose dirichlet's BCs to rhs of PPE 
   !
   subroutine set_poisson_solver_bcs ( lo, hi, rhs, slo, shi, phi, &
        & oro_x, ulo, uhi, oro_y, vlo, vhi, oro_z, wlo, whi,       &
        & bc_lo, bc_hi, dx, domlo, domhi ) bind(C)

      use amrex_lo_bctypes_module
      
      ! Loop bounds
      integer(c_int), intent(in   ) :: lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3)
      integer(c_int), intent(in   ) :: vlo(3), vhi(3)
      integer(c_int), intent(in   ) :: wlo(3), whi(3)

      ! Boundary conditions
      integer(c_int), intent(in   ) :: bc_lo(3), bc_hi(3)

      ! Grid bounds
      integer(c_int), intent(in   ) :: domlo(3), domhi(3)

      ! Grid spacings
      real(ar),       intent(in   ) :: dx(3)

      ! Arrays
      real(ar),       intent(in   ) ::                         &
           & oro_x(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           & oro_y(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
           & oro_z(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      
      real(ar),       intent(inout) ::                       &
           & rhs(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           &   phi(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      
      ! Local variables
      integer(c_int)                :: i, j, k
      integer(c_int)                :: nlft, nbot, ndwn
      integer(c_int)                :: nrgt, ntop, nup
      real(ar)                      :: odx2, ody2, odz2


      nlft = max(0,domlo(1)-slo(1))
      nbot = max(0,domlo(2)-slo(2))
      ndwn = max(0,domlo(3)-slo(3))

      nrgt = max(0,shi(1)-domhi(1))
      ntop = max(0,shi(2)-domhi(2))
      nup  = max(0,shi(3)-domhi(3))
      
      odx2 = one / dx(1)**2
      ody2 = one / dx(2)**2
      odz2 = one / dx(3)**2

      ! West boundary
      if ( (bc_lo(1) == amrex_lo_dirichlet) .and. (nlft > 0) ) then
         i = domlo(1) 
         do k = slo(3), shi(3)
            do j = slo(2), shi(2)
               print *, "PHI AT BOUNDARY ", i-1,j,k,phi(i-1,j,k)
               rhs(i,j,k)   = rhs(i,j,k) - 270.0D0 * odx2 !phi(i-1,j,k) * odx2 * half !* oro_x(i,j,k)
               phi(i-1,j,k) = zero
            end do
         end do
      end if

      ! ! East boundary
      ! if ( (bc_hi(1) == amrex_lo_dirichlet) .and. (nrgt > 0) ) then
      !    i = domhi(1)
      !    do k = slo(3), shi(3)
      !       do j = slo(2), shi(2)
      !          rhs(i,j,k)   = rhs(i,j,k) - phi(i+1,j,k) * odx2 * oro_x(i+1,j,k)
      !          phi(i+1,j,k) = zero
      !       end do
      !    end do
      ! end if

      ! ! South boundary
      ! if ( (bc_lo(2) == amrex_lo_dirichlet) .and. (nbot > 0) ) then
      !    j = domlo(2)
      !    do k = slo(3), shi(3)
      !       do i = slo(1), shi(1)
      !          rhs(i,j,k)   = rhs(i,j,k) - phi(i,j-1,k) * ody2 * oro_y(i,j,k)
      !          phi(i,j-1,k) = zero
      !       end do
      !    end do
      ! end if

      ! ! North boundary
      ! if ( (bc_hi(2) == amrex_lo_dirichlet) .and. (ntop > 0) ) then
      !    j = domhi(2)
      !    do k = slo(3), shi(3)
      !       do i = slo(1), shi(1)
      !          rhs(i,j,k)   = rhs(i,j,k) - phi(i,j+1,k) * ody2 * oro_y(i,j+1,k)
      !          phi(i,j+1,k) = zero
      !       end do
      !    end do
      ! end if

      ! ! Bottom boundary
      ! if ( (bc_lo(3) == amrex_lo_dirichlet) .and. (ndwn > 0) ) then
      !    k = domlo(3)
      !    do j = slo(2), shi(2)
      !       do i = slo(1), shi(1)
      !          rhs(i,j,k)   = rhs(i,j,k) - phi(i,j,k-1) * odz2 * oro_z(i,j,k)
      !          phi(i,j,k-1) = zero 
      !       end do
      !    end do
      ! end if

      ! ! Top boundary
      ! if ( (bc_hi(3) == amrex_lo_dirichlet) .and. (nup > 0) ) then
      !    k = domhi(3)
      !    do j = slo(2), shi(2)
      !       do i = slo(1), shi(1)
      !          rhs(i,j,k)   = rhs(i,j,k) - phi(i,j,k+1) * odz2 * oro_z(i,j,k+1)
      !          phi(i,j,k+1) = zero 
      !       end do
      !    end do
      ! end if
 
      
   end subroutine set_poisson_solver_bcs





   !
   ! Set the value of the auxiliary function PHI before the solution of the
   ! poisson equation. The pressure MUST have the correct boundary conditions
   ! 
   !
   ! WARNING: this routine MUST be called before solving the Poisson equation
   !
   subroutine set_phi ( lo, hi, phi, slo, shi, pg, scale, singular,  domlo, domhi ) bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) :: lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: slo(3), shi(3)

      ! Scaling factor
      real(ar),       intent(in   ) :: scale

      ! Whether the sistem is singular
      integer(c_int), intent(in   ) :: singular

      ! Grid bounds
      integer(c_int), intent(in   ) :: domlo(3), domhi(3)

      ! Arrays
      real(ar),       intent(in   ) :: &
           & pg(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),       intent(inout) :: &
           & phi(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Local variables
      integer(c_int)                :: is, ie, js, je, ks, ke


      is = lo(1)
      ie = hi(1)
      js = lo(2)
      je = hi(2)
      ks = lo(3)
      ke = hi(3)

      
      ! First set the interior points to the current value of pressure
      ! scaled by "scale"
      if ( singular == 0 ) then 
         phi(is:ie,js:je,ks:ke) = scale * pg(is:ie,js:je,ks:ke)
      else
         phi(is:ie,js:je,ks:ke) = zero
         return 
      end if
         
      ! Next, assign the boundary values.
      ! Here we assume that all the boundaries are Dirichlet's
      ! and set the boundary nodes accordingly.
      ! The Poisson's solver will ignore all the values set when
      ! the boundary is not of Dirichlet's type.
      ! The boundary cell will hold the value at the boundary location, NOT the
      ! value at the boundary CELL!!!
      if ( domlo(1) > slo(1) )  then 
         phi(domlo(1)-1,js:je,ks:ke) =  scale * half *  &
              & ( pg(domlo(1)-1,js:je,ks:ke) + pg(domlo(1),js:je,ks:ke) )
      end if 

      if ( domhi(1) < shi(1) )  then 
         phi(domhi(1)+1,js:je,ks:ke) =  scale * half *  &
              & ( pg(domhi(1)+1,js:je,ks:ke) + pg(domhi(1),js:je,ks:ke) )
      end if 

      if ( domlo(2) > slo(2) )  then 
         phi(is:ie,domlo(2)-1,ks:ke) =  scale * half *  &
              & ( pg(is:ie,domlo(2)-1,ks:ke) + pg(is:ie,domlo(2),ks:ke) )
      end if 

      if ( domhi(2) < shi(2) )  then 
         phi(is:ie,domhi(2)+1,ks:ke) =  scale * half *  &
              & ( pg(is:ie,domhi(2)+1,ks:ke) + pg(is:ie,domhi(2),ks:ke) )
      end if 

      if ( domlo(3) > slo(3) )  then 
         phi(is:ie,js:je,domlo(3)-1) =  scale * half *  &
              & ( pg(is:ie,js:je,domlo(3)-1) + pg(is:ie,js:je,domlo(3)) )
      end if 

      if ( domhi(3) < shi(3) )  then 
         phi(is:ie,js:je,domhi(3)+1) =  scale * half *  &
              & ( pg(is:ie,js:je,domhi(3)+1) + pg(is:ie,js:je,domhi(3)) )
      end if 
      
   end subroutine set_phi

   
end module projection_mod
