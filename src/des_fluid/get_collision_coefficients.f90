subroutine get_collision_model( flag_out ) bind(C, name="get_collision_model")

  use iso_c_binding,                 only: c_int
  use discretelement,                only: des_coll_model_enum

  implicit none

  integer(c_int), intent(inout) :: flag_out

  flag_out = des_coll_model_enum

end subroutine get_collision_model

subroutine get_lsd_collision_coefficients ( nphase_out, kt_out, kt_w_out, kn_out, kn_w_out, &
     mew_out, mew_w_out, etan_out, etan_w_out, etat_out, & 
     etat_w_out ) &
     bind(C, name="get_lsd_collision_coefficients")

   use amrex_fort_module,               only: rt => amrex_real
   use iso_c_binding,                   only: c_int
   use param,                           only: dim_m
   use constant,                        only: mmax
   use discretelement,                   only: kn, kn_w, kt, kt_w, kt_fac, kt_w_fac, mew, mew_w, &
                                              des_etan, des_etan_wall, des_etat, des_etat_wall
                                              
   implicit none

   real(rt),       intent(inout) :: kt_out, kt_w_out, kn_out, kn_w_out, mew_out, mew_w_out
   real(rt),       intent(inout) :: etan_out(dim_m, dim_m), etan_w_out(dim_m), &
                                    etat_out(dim_m, dim_m), etat_w_out(dim_m)
   integer(c_int), intent(inout) :: nphase_out

   integer :: i, j

   nphase_out = mmax
   kt_out = kt
   kt_w_out = kt_w
   kn_out = kn
   kn_w_out = kn_w
   mew_out = mew
   mew_w_out = mew_w

   ! convert from Fortran to C ordering here
   do i = 1, dim_m
      do j = 1, dim_m
         print *, des_etan(i, j)
         etan_out(i, j) = des_etan(j, i)
         etat_out(i, j) = des_etat(j, i)
      end do
   end do

   etan_w_out = des_etan_wall
   etat_w_out = des_etat_wall

end subroutine get_lsd_collision_coefficients
