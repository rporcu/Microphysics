!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  module name: init_collision                                         !
!                                                                      !
!  Purpose: DES - allocating DES arrays                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine init_collision(min_dp_in, min_ro_in, &
     &                    max_dp_in, max_ro_in, &
     &                    avg_dp_in, avg_ro_in, &
     &                    tcoll_ratio)    &
     bind(C, name="init_collision")

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding,     only: c_int
  use param,             only: zero, dim_m
  use constant,          only: mmax
  use discretelement,    only: dp_max, dp_min, dp_avg
  use discretelement,    only: ro_max, ro_min, ro_avg

  implicit none

  real(rt), intent(in) :: min_dp_in(dim_m), min_ro_in(dim_m)
  real(rt), intent(in) :: max_dp_in(dim_m), max_ro_in(dim_m)
  real(rt), intent(in) :: avg_dp_in(dim_m), avg_ro_in(dim_m)
  real(rt), intent(in) :: tcoll_ratio

  real(rt)             :: d_p0(dim_m),   ro_s0(dim_m)

  integer :: ptype

  ! Work around for cases that have no particles.
  do ptype =1, mmax

     d_p0(ptype)  = merge( min_dp_in(ptype),  100.0d-6, min_dp_in(ptype) > zero)
     ro_s0(ptype) = merge( max_ro_in(ptype), 1000.0d+0, max_ro_in(ptype) > zero)

     dp_max(ptype) = max_dp_in(ptype)
     dp_min(ptype) = min_dp_in(ptype)
     dp_avg(ptype) = avg_dp_in(ptype)

     ro_max(ptype) = max_ro_in(ptype)
     ro_min(ptype) = min_ro_in(ptype)
     ro_avg(ptype) = avg_ro_in(ptype)

  enddo

  call init_collision_lsd

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
  subroutine init_collision_lsd

  use constant,       only: mmax
  use amrex_constants_module, only: M_PI
  use discretelement, only: kn, kn_w, kt, kt_w, kt_fac, kt_w_fac, &
      & des_etan, des_etan_wall, des_etat, des_etat_wall,        &
      & des_en_input, des_en_wall_input, &
      & des_etat_fac, des_etat_w_fac, dtsolid

  integer      :: m, l, lc
  real(rt) :: tcoll, tcoll_tmp
  real(rt) :: mass_m, mass_l, mass_eff
  real(rt) :: en

  tcoll = 1.0d0

  ! Calculate the particle-particle tangential spring factor.
  kt = kt_fac*kn

  ! Calculate the particle-wall tangential spring factor.
  kt_w = kt_w_fac*kn_w

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
