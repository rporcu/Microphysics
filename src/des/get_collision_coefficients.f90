subroutine set_lsd_collision_coefficients ( mew_in,      mew_w_in,     &
     &                                      kn_in,       kn_w_in,      &
     &                                      kt_in,       kt_w_in,      &
     &                                      en_in,       en_w_in,      &
     &                                      kt_fac_in,   kt_w_fac_in,  &
     &                                      eta_fac_in,  eta_w_fac_in) &
     &                                      bind(C, name="set_lsd_collision_coefficients")

   use amrex_fort_module, only: rt => amrex_real
   use iso_c_binding,     only: c_int

   use param,             only: dim_m  ! Maximal number of particle types

   use discretelement,    only: mew,          mew_w          ! Friction coeffs
   use discretelement,    only: kn,           kn_w           ! Spring constants (normal)
   use discretelement,    only: kt,           kt_w           ! Spring constants (tangential)
   use discretelement,    only: des_en_input, des_en_wall_input ! Restitution coeffs
   use discretelement,    only: kt_fac,       kt_w_fac          ! tan/norm spring factor
   use discretelement,    only: des_etat_fac, des_etat_w_fac    ! tan/norm damping factor


   implicit none

   integer, parameter :: dim_lm = dim_m+dim_m*(dim_m-1)/2;

   real(rt),       intent(in   ) :: mew_in,        mew_w_in
   real(rt),       intent(in   ) :: kn_in,         kn_w_in
   real(rt),       intent(in   ) :: kt_in,         kt_w_in
   real(rt),       intent(in   ) :: en_in(dim_lm), en_w_in(dim_m)
   real(rt),       intent(in   ) :: kt_fac_in,     kt_w_fac_in
   real(rt),       intent(in   ) :: eta_fac_in,    eta_w_fac_in

   ! Copy inputs into fortran module variables. This is a temporary
   ! work around until we completely use c++.

   mew   = mew_in
   mew_w = mew_w_in

   kn   = kn_in
   kn_w = kn_w_in

   kt   = kt_in
   kt_w = kt_w_in

   des_en_input(:dim_lm)     = en_in(:dim_lm)
   des_en_wall_input(:dim_m) = en_w_in(:dim_m)

   kt_fac   = kt_fac_in
   kt_w_fac = kt_w_fac_in

   des_etat_fac   = eta_fac_in
   des_etat_w_fac = eta_w_fac_in

end subroutine set_lsd_collision_coefficients
