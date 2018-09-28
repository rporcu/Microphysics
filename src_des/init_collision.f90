!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  module name: init_collision                                         !
!                                                                      !
!  Purpose: DES - allocating DES arrays                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine init_collision(d_p_in, ro_s_in)&
     bind(C, name="init_collision")

  use discretelement,    only: des_coll_model_enum, lsd, hertzian

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding,     only: c_int
  use param,             only: zero, dim_m
  use constant,          only: mmax

  implicit none

  real(rt), intent(in) :: d_p_in(dim_m), ro_s_in(dim_m)
  real(rt)             :: d_p0(dim_m),   ro_s0(dim_m)

  integer :: ptype

   d_p0 =  d_p_in
  ro_s0 = ro_s_in

  ! Work around for cases that have no particles.
  do ptype =1, mmax

     if( d_p_in(ptype) == zero)  d_p0(ptype) =  100.0d-6
     if(ro_s_in(ptype) == zero) ro_s0(ptype) = 1000.0d0

  enddo

  select case (des_coll_model_enum)
     case(lsd)     ; call init_collision_lsd
     case(hertzian); call init_collision_hertz
  end select

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

  use constant,       only: pi, mmax
  use discretelement, only: kn, kn_w, kt, kt_w, kt_fac, kt_w_fac, &
      & des_etan, des_etan_wall, des_etat, des_etat_wall,        &
      & des_en_input, des_en_wall_input, &
      & des_etat_fac, des_etat_w_fac, dtsolid
!     & des_et_input, des_et_wall_input

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

     mass_m = (pi/6.d0)*(d_p0(m)**3)*ro_s0(m)

     ! Particle-Particle Collision Parameters
     do l = m, mmax
        lc = lc+1

        en = des_en_input(lc)

        ! Calculate masses used for collision calculations.
        mass_l = (pi/6.d0)*(d_p0(l)**3)*ro_s0(l)
        mass_eff = mass_m*mass_l/(mass_m+mass_l)

        ! Calculate the M-L normal and tangential damping coefficients.
        if(abs(en) > zero) then
           des_etan(m,l) = 2.0d0*sqrt(kn*mass_eff) * abs(log(en))
           des_etan(m,l) = des_etan(m,l)/sqrt(pi*pi + (log(en)**2))
        else
           des_etan(m,l) = 2.0d0*sqrt(kn*mass_eff)
        endif
        des_etat(m,l) = des_etat_fac*des_etan(m,l)

        ! Store the entries in the symmetric matrix.
        des_etan(l,m) = des_etan(m,l)
        des_etat(l,m) = des_etat(m,l)

        ! Calculate the collision time scale.
        tcoll_tmp = pi/sqrt(kn/mass_eff -                          &
          ((des_etan(m,l)/mass_eff)**2)/4.d0)
        tcoll = min(tcoll_tmp, tcoll)
     end do

     ! Particle-Wall Collision Parameters
     en = des_en_wall_input(m)
     mass_eff = mass_m

     ! Calculate the M-Wall normal and tangential damping coefficients.
     if(abs(en) > zero) then
        des_etan_wall(m) = 2.d0*sqrt(kn_w*mass_eff)*abs(log(en))
        des_etan_wall(m) = des_etan_wall(m)/sqrt(pi*pi+(log(en))**2)
     else
        des_etan_wall(m) = 2.d0*sqrt(kn_w*mass_eff)
     endif
     des_etat_wall(m) = des_etat_w_fac*des_etan_wall(m)

     ! Calculate the collision time scale.
     tcoll_tmp = pi/sqrt(kn_w/mass_eff -                           &
       ((des_etan_wall(m)/mass_eff)**2.d0)/4.d0)
     tcoll = min(tcoll_tmp, tcoll)
  enddo

  ! Store the smalled calculated collision time scale.
  dtsolid = tcoll/50.d0

  end subroutine init_collision_lsd

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: init_collision_hertz                                    !
!                                                                      !
!  References:                                                         !
!   - Schafer et al., J. Phys. I France, 1996, 6, 5-20 (see page 7&13) !
!   -  Van der Hoef et al., Advances in Chemical Engineering, 2006, 31,!
!      65-149 (pages 94-95)                                            !
!   - Silbert et al., Physical Review E, 2001, 64, 051302 1-14 (page 5)!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine init_collision_hertz

  use constant,       only: mmax, pi
  use param,          only: dim_m
  use discretelement, only: des_en_input, des_en_wall_input,   &
                              des_et_input, des_et_wall_input, &
                              des_etan, des_etan_wall,         &
                              des_etat, des_etat_wall,         &
                              hert_kn, hert_kwn,               &
                              hert_kt, hert_kwt,               &
                              e_young, ew_young,               &
                              v_poisson, vw_poisson,           &
                              dtsolid
!                             des_etat_fac, des_etat_w_fac

    integer           :: m, l, lc
    real(rt)      :: tcoll, tcoll_tmp
    ! Particle and effective mass.
    real(rt)      :: mass_m, mass_l, mass_eff
    ! Effective physical quantities. Radius, Youngs, Shear
    real(rt)      :: r_eff, e_eff, g_mod_eff, red_mass_eff
    ! Alias for coefficient restitution
    real(rt)      :: en, et
    ! Shear modules for particles and wall
    real(rt)      :: g_mod(dim_m), g_mod_wall

    tcoll = 1.0d0

    do m=1,mmax
       g_mod(m) = 0.5d0*e_young(m)/(1.d0+v_poisson(m))
    enddo
    g_mod_wall = 0.5d0*ew_young/(1.d0+vw_poisson)

    lc = 0
    do m=1,mmax
       ! Calculate the mass of a phase M particle.
       mass_m = (pi/6.d0)*(d_p0(m)**3)*ro_s0(m)

       do l=m,mmax
          lc = lc+1

          en = des_en_input(lc)
          et = des_et_input(lc)

          ! Calculate masses used for collision calculations.
          mass_l = (pi/6.d0)*(d_p0(l)**3)*ro_s0(l)
          mass_eff = (mass_m*mass_l)/(mass_m+mass_l)
          red_mass_eff = (2.d0/7.d0)*mass_eff
          ! Calculate the effective radius, Youngs modulus, and shear modulus.
          r_eff = 0.5d0*(d_p0(m)*d_p0(l)/ (d_p0(m) + d_p0(l)))
          e_eff = e_young(m)*e_young(l) /                            &
               (e_young(m)*(1.d0 - v_poisson(l)**2) +                &
               e_young(l)*(1.d0 - v_poisson(m)**2))
          g_mod_eff = g_mod(m)*g_mod(l)/                             &
               (g_mod(m)*(2.d0 - v_poisson(l)) +                     &
               g_mod(l)*(2.d0 - v_poisson(m)))

          ! Calculate the spring properties and store in symmetric matrix format.
          hert_kn(m,l)=(4.d0/3.d0)*sqrt(r_eff)*e_eff
          hert_kt(m,l)= 8.d0*sqrt(r_eff)*g_mod_eff

          hert_kn(l,m) = hert_kn(m,l)
          hert_kt(l,m) = hert_kt(m,l)

          ! Calculate the normal coefficient.
          if(abs(en) > zero) then
             des_etan(m,l) = 2.d0*sqrt(hert_kn(m,l)*mass_eff)*abs(log(en))
             des_etan(m,l) = des_etan(m,l)/sqrt(pi*pi + (log(en))**2)
          else
             des_etan(m,l) = 2.d0*sqrt(hert_kn(m,l)*mass_eff)
          endif
          des_etan(l,m) = des_etan(m,l)

          ! Calculate the tangential coefficients.
          if(abs(et) > zero) then
             des_etat(m,l) = 2.d0*sqrt(hert_kt(m,l)*red_mass_eff)*abs(log(et))
             des_etat(m,l) = des_etat(m,l)/ sqrt(pi*pi + (log(et))**2)
          else
             des_etat(m,l) = 2.d0*sqrt(hert_kt(m,l)*red_mass_eff)
          endif
          des_etat(l,m) = des_etat(m,l)

          tcoll_tmp = pi/sqrt(hert_kn(m,l)/mass_eff -                &
               ((des_etan(m,l)/mass_eff)**2)/4.d0)
          tcoll = min(tcoll_tmp, tcoll)
       enddo

       en = des_en_wall_input(m)
       et = des_et_wall_input(m)

       ! Calculate masses used for collision calculations.
       mass_eff = mass_m
       red_mass_eff = (2.d0/7.d0)*mass_eff
       ! Calculate the effective radius, Youngs modulus, and shear modulus.
       r_eff = 0.5d0*d_p0(m)
       e_eff = e_young(m)*ew_young /                                 &
            (e_young(m)*(1.d0-vw_poisson**2) +                       &
            ew_young  *(1.d0-v_poisson(m)**2))
       g_mod_eff = g_mod(m)*g_mod_wall /                             &
            (g_mod(m)*(2.d0 - vw_poisson) +                          &
            g_mod_wall*(2.d0 - v_poisson(m)))

       ! calculate the spring properties.
       hert_kwn(m) = (4.d0/3.d0)*sqrt(r_eff)*e_eff
       hert_kwt(m) = 8.0*sqrt(r_eff)*g_mod_eff

       ! Calculate the tangential coefficients.
       if(abs(en) > zero) then
          des_etan_wall(m) = 2.d0*sqrt(hert_kwn(m)*mass_eff)*&
               abs(log(en))
          des_etan_wall(m) = des_etan_wall(m)/&
               sqrt(pi*pi + (log(en))**2)
       else
          des_etan_wall(m) = 2.d0*sqrt(hert_kwn(m)*mass_eff)
       endif

       if(abs(et) > zero) then
          des_etat_wall(m) = 2.d0*sqrt(hert_kwt(m)*red_mass_eff)*    &
               abs(log(et))
          des_etat_wall(m) = des_etat_wall(m)/sqrt(pi*pi+(log(et))**2)
       else
          des_etat_wall(m) = 2.d0*sqrt(hert_kwt(m)*red_mass_eff)
       endif

       ! Calculate the collision time scale.
       tcoll_tmp = pi/sqrt(hert_kwn(m)/mass_eff -                    &
            ((des_etan_wall(m)/mass_eff)**2)/4.d0)
       tcoll = min(tcoll_tmp, tcoll)
    enddo

    dtsolid = tcoll/50.d0

  end subroutine init_collision_hertz

end subroutine init_collision


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! Procedure Name: sum_particle_props                                   !
!                                                                      !
! Purpose: Sum diameter, density and count number of particles. This   !
! is used to calculate the average particle diameter and density which !
! go into calculating the particle collision time scale.               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine sum_particle_props ( np, particles, sum_np, sum_dp, sum_ro) &
     bind(C, name="sum_particle_props")

  use amrex_fort_module, only: rt => amrex_real
  use iso_c_binding ,    only: c_int

  use param, only: dim_m
  use particle_mod

  integer(c_int),   intent(in   ) :: np
  type(particle_t), intent(in   ) :: particles(np)

  real(rt), intent(inout) :: sum_np(dim_m)
  real(rt), intent(inout) :: sum_dp(dim_m)
  real(rt), intent(inout) :: sum_ro(dim_m)

  integer :: p, m

  do p=1, np
     m = particles(p) % phase
     sum_np(m) = sum_np(m) + 1.0d0
     sum_dp(m) = sum_dp(m) + particles(p) % radius * 2.0d0
     sum_ro(m) = sum_ro(m) + particles(p) % density
  enddo

end subroutine sum_particle_props
