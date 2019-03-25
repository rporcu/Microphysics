
subroutine get_collision_coefficients ( nphase_out, kt_out, kt_w_out, kn_out, kn_w_out, &
                                        etan_out, etan_w_out, etat_out, etat_w_out ) &
  bind(C, name="get_collision_coefficients")

   use amrex_fort_module,               only: rt => amrex_real
   use iso_c_binding,                   only: c_int
   use param,                           only: dim_m
   use constant,                        only: mmax
   use discretelement,                   only: kn, kn_w, kt, kt_w, kt_fac, kt_w_fac, &
                                              des_etan, des_etan_wall, des_etat, des_etat_wall
                                              
   implicit none

   real(rt),       intent(inout) :: kt_out, kt_w_out, kn_out, kn_w_out
   real(rt),       intent(inout) :: etan_out(dim_m, dim_m), etan_w_out(dim_m), &
                                    etat_out(dim_m, dim_m), etat_w_out(dim_m)
   integer(c_int), intent(inout) :: nphase_out

   nphase_out = mmax
   etan_out = des_etan
   etan_w_out = des_etan_wall
   etat_out = des_etat
   etat_w_out = des_etat_wall

end subroutine get_collision_coefficients
