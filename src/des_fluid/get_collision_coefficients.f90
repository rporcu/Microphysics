subroutine get_lsd_collision_coefficients ( nphase_out, kt_out, kt_w_out, kn_out, kn_w_out, &
     mew_out, mew_w_out, etan_out, etan_w_out, etat_out, &
     etat_w_out, neighborhood ) &
     bind(C, name="get_lsd_collision_coefficients")

   use amrex_fort_module, only: rt => amrex_real
   use iso_c_binding,     only: c_int
   use param,             only: dim_m
   use constant,          only: mmax
   use discretelement,    only: kn, kn_w, kt, kt_w, mew, mew_w, &
        &                       des_etan, des_etan_wall, des_etat, des_etat_wall, &
        &                       dp_max


   implicit none

   real(rt),       intent(inout) :: kt_out, kt_w_out, kn_out, kn_w_out, mew_out, mew_w_out
   real(rt),       intent(inout) :: etan_out(dim_m, dim_m), etan_w_out(dim_m), &
        &                           etat_out(dim_m, dim_m), etat_w_out(dim_m), &
        &                           neighborhood

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
         etan_out(i, j) = des_etan(j, i)
         etat_out(i, j) = des_etat(j, i)
      end do
   end do

   etan_w_out = des_etan_wall
   etat_w_out = des_etat_wall

   neighborhood = (3.0d0*maxval(dp_max(1:mmax)))**2

end subroutine get_lsd_collision_coefficients
