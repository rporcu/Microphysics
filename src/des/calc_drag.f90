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
subroutine calc_drag_fluid(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, &
  max_pip, ep_g, ro_g, u_g, v_g, w_g, mu_g, f_gs, rhs, &
  pphase, pstate, pvol, ppos, pvel, pradius, drag, dx, dy, dz)&
  bind(C, name="calc_drag_fluid")

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int

  use des_drag_gp_module, only: des_drag_gp

  implicit none

  integer(c_int), intent(in   ) :: slo(3),shi(3)
  integer(c_int), intent(in   ) :: ulo(3),uhi(3)
  integer(c_int), intent(in   ) :: vlo(3),vhi(3)
  integer(c_int), intent(in   ) :: wlo(3),whi(3)
  integer(c_int), intent(in   ) :: max_pip

  integer(c_int), intent(in   ) :: pstate(max_pip)
  integer(c_int), intent(in   ) :: pphase(max_pip)

  real(c_real), intent(in   ) :: &
    ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
    ro_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
    mu_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
     u_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
     v_g(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
     w_g(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

  real(c_real), intent(in   ) :: dx, dy, dz, &
    pvol(max_pip),  ppos(max_pip,3),         &
    pradius(max_pip),  pvel(max_pip,3)


  real(c_real), intent(out  ) :: drag(max_pip,3),    &
    f_gs(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
     rhs(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

! Local variables
!---------------------------------------------------------------------//
! Loop counters: Particle, fluid cell, neighbor cells
  integer :: np, i, j, k
  real(c_real) :: velfp(3), velp(3), beta, f_gp
  real(c_real) :: odx, ody, odz, ovol
!......................................................................!

  odx = 1.0d0/dx
  ody = 1.0d0/dy
  odz = 1.0d0/dz

  ovol = 1.0d0/(dx*dy*dz)

! Calculate the gas phase forces acting on each particle.

  do np=1,max_pip

     i = floor(ppos(np,1)*odx)
     j = floor(ppos(np,2)*ody)
     k = floor(ppos(np,3)*odz)

     ! Gas volume fraction, velocity, at the particle's position.
     velfp(1) = 0.5d0*(u_g(i,j,k) + u_g(i+1,j,k))
     velfp(2) = 0.5d0*(v_g(i,j,k) + v_g(i,j+1,k))
     velfp(3) = 0.5d0*(w_g(i,j,k) + w_g(i,j,k+1))

     velp(:) = pvel(np,:)

     ! Calculate drag coefficient, beta
     call des_drag_gp(slo, shi, np, velp, velfp, ep_g(i,j,k), &
       ro_g, mu_g, beta, i, j, k, pradius(np), pvol(np), pphase(np))

     f_gp = beta*ovol

     f_gs(i,j,k)  = f_gs(i,j,k)  + f_gp
     rhs(i,j,k,:) = rhs(i,j,k,:) + f_gp*velp(:)

     drag(np,1) = beta

  enddo

end subroutine calc_drag_fluid

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_drag_particle                                      !
!                                                                      !
!  Purpose: This routine is called before the FLUID solve.             !
!  It calculates the source terms for the center coefficients and RHS  !
!  for the momentum equations. It also saves the drag coefficient for  !
!  each particle.                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine calc_drag_particle(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, &
  max_pip, p_g, u_g, v_g, w_g, pvol, ppos, pvel, drag, &
  dx, dy, dz, xlen, ylen, zlen) bind(C, name="calc_drag_particle")

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int

  use des_drag_gp_module, only: des_drag_gp
  ! Flags for cyclic BC with pressure drop
  use bc, only: cyclic_x_pd, cyclic_y_pd, cyclic_z_pd
  ! Specified pressure drop
  use bc, only: delp_x, delp_y, delp_z

  implicit none

  integer(c_int), intent(in   ) :: slo(3),shi(3), &
    ulo(3),uhi(3), vlo(3),vhi(3), wlo(3),whi(3), max_pip

  real(c_real), intent(in   ) :: &
    p_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
    u_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
    v_g(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
    w_g(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

  real(c_real), intent(in   ) :: dx, dy, dz, xlen, ylen, zlen, &
    pvol(max_pip), ppos(max_pip,3), pvel(max_pip,3)

  real(c_real), intent(inout) :: drag(max_pip,3)

! Local variables
!---------------------------------------------------------------------//
! Loop counters: Particle, fluid cell, neighbor cells
  integer :: np, i, j, k
  real(c_real) :: velfp(3), velp(3), cpg(3), gradpg(3)
  real(c_real) :: beta(max_pip)
  real(c_real) :: odx, ody, odz
!......................................................................!

  odx = 1.0d0/dx
  ody = 1.0d0/dy
  odz = 1.0d0/dz

  beta = drag(:,1)

  cpg(1) = merge(delp_x/xlen, 0.0d0, cyclic_x_pd)
  cpg(2) = merge(delp_y/ylen, 0.0d0, cyclic_y_pd)
  cpg(3) = merge(delp_z/zlen, 0.0d0, cyclic_z_pd)

! Calculate the gas phase forces acting on each particle.

  do np=1,max_pip

     i = floor(ppos(np,1)*odx)
     j = floor(ppos(np,2)*ody)
     k = floor(ppos(np,3)*odz)

     velfp(1) = 0.5d0*(u_g(i,j,k) + u_g(i+1,j,k))
     velfp(2) = 0.5d0*(v_g(i,j,k) + v_g(i,j+1,k))
     velfp(3) = 0.5d0*(w_g(i,j,k) + w_g(i,j,k+1))

     gradpg(1) = odx*(0.5d0*(p_g(i,j,k) + p_g(i+1,j,k)) - &
                      0.5d0*(p_g(i,j,k) + p_g(i-1,j,k)))

     gradpg(2) = ody*(0.5d0*(p_g(i,j,k) + p_g(i,j+1,k)) - &
                      0.5d0*(p_g(i,j,k) + p_g(i,j-1,k)))

     gradpg(3) = odz*(0.5d0*(p_g(i,j,k) + p_g(i,j,k+1)) - &
                      0.5d0*(p_g(i,j,k) + p_g(i,j,k-1)))

     drag(np,:) = beta(np)*(velfp - pvel(np,:)) + &
                     (cpg(:) - gradpg(:))*pvol(np)

  enddo

end subroutine calc_drag_particle
