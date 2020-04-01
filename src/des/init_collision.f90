!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  module name: init_collision                                         !
!                                                                      !
!  Purpose: DES - allocating DES arrays                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine init_collision (mmax,                  &
     &                     min_dp_in, min_ro_in,  &
     &                     max_dp_in, max_ro_in,  &
     &                     avg_dp_in, avg_ro_in,  &
     &                     tcoll_ratio,           &
     &                     etan_out,  etan_w_out, &
     &                     etat_out, etat_w_out, &
     &                     neighborhood ) &
     bind(C, name="init_collision")

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding,     only: c_int

  use param,             only: dim_m
  use param,             only: zero

  use discretelement,    only: dp_max, dp_min, dp_avg
  use discretelement,    only: ro_max, ro_min, ro_avg

  use discretelement,    only: des_etan, des_etan_wall
  use discretelement,    only: des_etat, des_etat_wall

  implicit none

  integer(c_int), intent(in) :: mmax

  real(rt), intent(in) :: min_dp_in(dim_m), min_ro_in(dim_m)
  real(rt), intent(in) :: max_dp_in(dim_m), max_ro_in(dim_m)
  real(rt), intent(in) :: avg_dp_in(dim_m), avg_ro_in(dim_m)
  real(rt), intent(in) :: tcoll_ratio

  real(rt), intent(inout) :: etan_out(dim_m, dim_m), etan_w_out(dim_m)
  real(rt), intent(inout) :: etat_out(dim_m, dim_m), etat_w_out(dim_m)
  real(rt), intent(inout) :: neighborhood

  real(rt)             :: d_p0(dim_m),   ro_s0(dim_m)

  integer :: ptype, i, j

  if(mmax == 0) return

  ! Work around for cases that have no particles.
  do ptype = 1, mmax

     ! d_p0(ptype)  = merge( min_dp_in(ptype),  100.0d-6, min_dp_in(ptype) > zero)
     ! ro_s0(ptype) = merge( max_ro_in(ptype), 1000.0d+0, max_ro_in(ptype) > zero)

     d_p0(ptype)  = merge( avg_dp_in(ptype),  100.0d-6, avg_dp_in(ptype) > zero)
     ro_s0(ptype) = merge( avg_ro_in(ptype), 1000.0d+0, avg_ro_in(ptype) > zero)

     dp_max(ptype) = max_dp_in(ptype)
     dp_min(ptype) = min_dp_in(ptype)
     dp_avg(ptype) = avg_dp_in(ptype)

     ro_max(ptype) = max_ro_in(ptype)
     ro_min(ptype) = min_ro_in(ptype)
     ro_avg(ptype) = avg_ro_in(ptype)

  enddo

  call init_collision_lsd(mmax)

  ! convert from Fortran to C ordering here
  do i = 1, dim_m
     do j = 1, dim_m
        etan_out(i, j) = des_etan(j, i)
        etat_out(i, j) = des_etat(j, i)
     end do
  end do

  etan_w_out = des_etan_wall
  etat_w_out = des_etat_wall

  neighborhood = (3.0d0*maxval(dp_max(1:mmax)/2.0d0))**2

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: init_collision_lsd                                      !
!  Purpose: Check user input data for DES collision calculations.      !
!                                                                      !
!  References:                                                         !
!   - Schafer et al., J. Phys. I France, 1996, 6, 5-20 (see page 7&13) !
!   -  Van der Hoef et al., Advances in Chemical Engineering, 2006, 31,!
!      65-149 (pages 94-95)                                            !
!   - Silbert et al., Physical Review E, 2001, 64, 051302 1-14 (page 5)!
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine init_collision_lsd(mmax)

  use amrex_constants_module, only: M_PI
  use discretelement, only: kn, kn_w, &
      & des_etan, des_etan_wall, des_etat, des_etat_wall,        &
      & des_en_input, des_en_wall_input, &
      & des_etat_fac, des_etat_w_fac, dtsolid

  integer, intent(in   ) :: mmax

  integer  :: m, l, lc
  real(rt) :: tcoll, tcoll_tmp
  real(rt) :: mass_m, mass_l, mass_eff
  real(rt) :: en

  tcoll = 1.0d0

  lc = 0
  do m = 1, mmax

     mass_m = (M_PI/6.d0)*(d_p0(m)**3)*ro_s0(m)

     ! Particle-Particle Collision Parameters
     do l = m, mmax
        lc = lc+1

        en = des_en_input(lc)

        ! Calculate masses used for collision calculations.
        mass_l = (M_PI/6.d0)*(d_p0(l)**3)*ro_s0(l)
        mass_eff = mass_m*mass_l/(mass_m+mass_l)

        ! Calculate the M-L normal and tangential damping coefficients.
        if(abs(en) > zero) then
           des_etan(m,l) = 2.0d0*sqrt(kn*mass_eff) * abs(log(en))
           des_etan(m,l) = des_etan(m,l)/sqrt(M_PI*M_PI + (log(en)**2))
        else
           des_etan(m,l) = 2.0d0*sqrt(kn*mass_eff)
        endif
        des_etat(m,l) = des_etat_fac*des_etan(m,l)

        ! Store the entries in the symmetric matrix.
        des_etan(l,m) = des_etan(m,l)
        des_etat(l,m) = des_etat(m,l)

        ! Calculate the collision time scale.
        tcoll_tmp = M_PI/sqrt(kn/mass_eff -                          &
             ((des_etan(m,l)/mass_eff)**2)/4.d0)
        tcoll = min(tcoll_tmp, tcoll)
     end do

     ! Particle-Wall Collision Parameters
     en = des_en_wall_input(m)
     mass_eff = mass_m

     ! Calculate the M-Wall normal and tangential damping coefficients.
     if(abs(en) > zero) then
        des_etan_wall(m) = 2.d0*sqrt(kn_w*mass_eff)*abs(log(en))
        des_etan_wall(m) = des_etan_wall(m)/sqrt(M_PI*M_PI+(log(en))**2)
     else
        des_etan_wall(m) = 2.d0*sqrt(kn_w*mass_eff)
     endif
     des_etat_wall(m) = des_etat_w_fac*des_etan_wall(m)

     ! Calculate the collision time scale.
     tcoll_tmp = M_PI/sqrt(kn_w/mass_eff -                           &
       ((des_etan_wall(m)/mass_eff)**2.d0)/4.d0)
     tcoll = min(tcoll_tmp, tcoll)
  enddo

  ! Store the smalled calculated collision time scale.
  dtsolid = tcoll/tcoll_ratio

  end subroutine init_collision_lsd

end subroutine init_collision
