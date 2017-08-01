!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_drag_fluid                                         !
!                                                                      !
!  Purpose: This routine is called before the FLUID solve.             !
!  It calculates the source terms for the center coefficients and RHS  !
!  for the momentum equations. It also saves the drag coefficient for  !
!  each particle.                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine calc_drag_fluid ( slo, shi, ulo, uhi, vlo, vhi, wlo, whi,     &
     np, ep_g, ro_g, u_g, v_g, w_g, mu_g, f_gs, rhs, particles, dx, dy, dz, use_pic ) &
     bind(C, name="calc_drag_fluid")

   use amrex_fort_module,  only : c_real => amrex_real
   use iso_c_binding ,     only: c_int
   use des_drag_gp_module, only: des_drag_gp
   use particle_mod,      only: particle_t

   implicit none

   integer(c_int), intent(in   ) :: slo(3),shi(3)
   integer(c_int), intent(in   ) :: ulo(3),uhi(3)
   integer(c_int), intent(in   ) :: vlo(3),vhi(3)
   integer(c_int), intent(in   ) :: wlo(3),whi(3)
   integer(c_int), intent(in   ) :: np, use_pic

   real(c_real), intent(in   ) :: &
        ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
        ro_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
        mu_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
        u_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
        v_g(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
        w_g(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

   real(c_real), intent(in   ) :: dx, dy, dz

   type(particle_t), intent(inout) :: particles(np)
   real(c_real),     intent(out  ) :: &
        f_gs(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
        rhs(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

   ! Local variables
   !---------------------------------------------------------------------//
   ! Loop counters: Particle, fluid cell, neighbor cells
   integer :: p, i, j, k
   real(c_real) :: velfp(3), velp(3), beta, f_gp
   real(c_real) :: odx, ody, odz, ovol
   !......................................................................!

   odx = 1.0d0/dx
   ody = 1.0d0/dy
   odz = 1.0d0/dz

   ovol = 1.0d0/(dx*dy*dz)

   ! Calculate the gas phase forces acting on each particle.

   do p = 1, np

      i = floor(particles(p) % pos(1)*odx)
      j = floor(particles(p) % pos(2)*ody)
      k = floor(particles(p) % pos(3)*odz)

      ! Gas volume fraction, velocity, at the particle's position.
      velfp(1) = 0.5d0*(u_g(i,j,k) + u_g(i+1,j,k))
      velfp(2) = 0.5d0*(v_g(i,j,k) + v_g(i,j+1,k))
      velfp(3) = 0.5d0*(w_g(i,j,k) + w_g(i,j,k+1))

      velp(:) = particles(p) % vel

      ! Calculate drag coefficient, beta
      call des_drag_gp(slo, shi, p, velp, velfp, ep_g(i,j,k),  &
           ro_g, mu_g, beta, i, j, k, particles(p) % radius,   &
           particles(p) % volume, particles(p) % density,      &
           particles(p) % phase )

      particles(p) % drag(1) = beta

      if (use_pic .eq. 0) then
         f_gp = beta*ovol
         f_gs(i,j,k  )  = f_gs(i,j,k)  + f_gp
          rhs(i,j,k,:) = rhs(i,j,k,:) + f_gp*velp(:)
      end if

   end do

end subroutine calc_drag_fluid


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_drag_particle                                      !
!                                                                      !
!  Purpose: This routine must be called after calc_drag_fluid because  !
!  it assumes you have already computed beta and stored it in          !
!  particles(p)%drag(1).   Here we compute the actual drag on the      !
!  particle.                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine calc_drag_particle( slo, shi, ulo, uhi, vlo, vhi, wlo, whi, &
     np, p_g, u_g, v_g, w_g, particles, dx, dy, dz, xlen, ylen, zlen )&
     bind(C, name="calc_drag_particle")

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int
   use des_drag_gp_module, only: des_drag_gp
   ! Specified pressure drop
   use bc, only: delp_x, delp_y, delp_z
   use particle_mod, only: particle_t

   implicit none

   integer(c_int), intent(in   ) :: slo(3),shi(3), &
        ulo(3),uhi(3), vlo(3),vhi(3), wlo(3),whi(3), np

   real(c_real), intent(in   ) :: &
        p_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
        u_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
        v_g(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
        w_g(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

   real(c_real),     intent(in   ) :: dx, dy, dz, xlen, ylen, zlen
   type(particle_t), intent(inout) :: particles(np)


   ! Local variables
   !---------------------------------------------------------------------//
   ! Loop counters: Particle, fluid cell, neighbor cells
   integer :: p, i, j, k
   real(c_real) :: velfp(3), velp(3), cpg(3), gradpg(3)
   real(c_real) :: beta(np)
   real(c_real) :: odx, ody, odz
   !......................................................................!

   odx = 1.0d0/dx
   ody = 1.0d0/dy
   odz = 1.0d0/dz

   do p = 1, np
      beta(p) = particles(p) % drag(1)
   end do

   cpg(1) = delp_x/xlen
   cpg(2) = delp_y/ylen
   cpg(3) = delp_z/zlen

   ! Calculate the gas phase forces acting on each particle.

   do p = 1, np

      i = floor(particles(p) % pos(1)*odx)
      j = floor(particles(p) % pos(2)*ody)
      k = floor(particles(p) % pos(3)*odz)

      velfp(1) = 0.5d0*(u_g(i,j,k) + u_g(i+1,j,k))
      velfp(2) = 0.5d0*(v_g(i,j,k) + v_g(i,j+1,k))
      velfp(3) = 0.5d0*(w_g(i,j,k) + w_g(i,j,k+1))

      gradpg(1) = 0.5d0 * odx * ( p_g(i+1,j,k) - p_g(i-1,j,k) )

      gradpg(2) = 0.5d0 * ody * ( p_g(i,j+1,k) - p_g(i,j-1,k) )

      gradpg(3) = 0.5d0 * odz * ( p_g(i,j,k+1) - p_g(i,j,k-1) )

      particles(p) % drag = beta(p)*(velfp - particles(p) % vel) + &
           (cpg(:) - gradpg(:)) * particles(p) % volume

   enddo

end subroutine calc_drag_particle
