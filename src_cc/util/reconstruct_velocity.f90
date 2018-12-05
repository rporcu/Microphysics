!
! This subroutines reconstructs the velocity field in the covered
! cells within a narrow band around the EB walls
!
! Author: Michele Rosso, LBL
!
! Date:   November 15, 2018
!
subroutine reconstruct_velocity ( vel_out, volo, vohi,       &
 &                                 vel_in, vilo, vihi,       &
 &                                phi, phlo, phhi, n_refine, &
 &                                flags, flo, fhi,           &
 &                                x0, dx  ) bind(C)

   use amrex_error_module,      only: amrex_abort
   use interpolation_m,         only: trilinear_interp, interp_stencil_is_valid
   use amrex_fort_module,       only: rt => amrex_real
   use iso_c_binding,           only: c_int
   use particle_mod,            only: particle_t
   use param,                   only: half, one, two
   use amrex_ebcellflag_module, only: is_covered_cell
   use amrex_eb_levelset_module,only: amrex_eb_interp_levelset, &
    &                                 amrex_eb_normal_levelset

   implicit none

   ! Array bounds
   integer(c_int), intent(in   ) :: vilo(3), vihi(3)
   integer(c_int), intent(in   ) :: volo(3), vohi(3)
   integer(c_int), intent(in   ) :: phlo(3), phhi(3)
   integer(c_int), intent(in   ) ::  flo(3),  fhi(3)

   ! Grid info
   real(rt),       intent(in   ) :: dx(3)
   real(rt),       intent(in   ) :: x0(3)

   ! Level-set refinement level
   integer,        intent(in   ), value  :: n_refine

   ! Arrays
   real(rt),       intent(in   ) :: &
    &  vel_in(vilo(1):vihi(1),vilo(2):vihi(2),vilo(3):vihi(3), 1:3), &
    &     phi(phlo(1):phhi(1),phlo(2):phhi(2),phlo(3):phhi(3))  ! Phi is nodal

   integer(c_int), intent(in   ) :: &
    & flags(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

   real(rt),       intent(inout) :: &
    &  vel_out(volo(1):vohi(1),volo(2):vohi(2),volo(3):vohi(3), 1:3)

   ! Loop indeces
   integer             :: i, j, k

   ! Width of narrow band 
   real(rt), parameter :: band_width    = two
   real(rt)            :: phi_threshold 

   ! Amout of "correction" to go from mirror point to interpolation point
   real(rt), parameter :: eps = 0.1_rt

   ! Coordinates, level set, and normal of cell center
   real(rt)            :: x_cc(3), phi_cc, norm_cc(3)

   ! ! Coordinates, level set, and normal of mirror point
   ! real(rt)            :: x_m(3), phi_m, norm_m(3)

   ! Coordinates, level set, normal  and velocity at interpolation point
   real(rt)            :: x_i(3), phi_i, norm_i(3), vel_i(3)

   ! Other local variables
   integer, parameter  :: max_iter = 20 ! Max number of iterations to find mirror point
   integer             :: iter 
   real(rt)            :: odx(3)

   odx = one / dx
   phi_threshold = band_width * maxval(dx)


   ! We can avoid filling the first and last layer of cells since
   ! vel_in will be used for trilinear interpolation only and this require only
   ! 1 layer of ghost nodes ( and in EB land, the velocity field has multiple layers of
   ! ghost nodes)
   do k = volo(3)+1,vohi(3)-1
      do j = volo(2)+1, vohi(2)-1
         do i = volo(1)+1, vohi(1)-1

            vel_out(i,j,k,:) = vel_in(i,j,k,:)

            ! Only reconstruct the velocity value if a covered cell is
            ! inside a band of width band_width.
            ! To check this, we look up the values of phi for the cell.
            ! Note that phi is nodal, so we check all the 8 corners
            ! of the cell
            if ( is_covered_cell(flags(i,j,k))                          .and. &
             &   minval(abs(phi(i:i+1,j:j+1,k:k+1))) <= phi_threshold ) then

               ! Coordinates of cell center
               x_cc = ( real([i,j,k],rt) + half ) * dx

               ! Get phi at cell center
               call amrex_eb_interp_levelset(x_cc, x0, n_refine, phi, phlo, phhi, dx, phi_cc)

               ! Skip innermost loop if cell is too far from zero level set
               if (abs(phi_cc) > band_width*maxval(dx)) cycle
               
               ! Get normal at cell center
               call amrex_eb_normal_levelset(x_cc, x0, n_refine, phi, phlo, phhi, dx, norm_cc)

               ! Find location of interpolation point by iteration if necessary
               x_i  = x_cc + two * abs(phi_cc) * norm_cc

               iter = 0
               find_xi: do

                  if ( interp_stencil_is_valid(x_i, x0, dx, flags, flo, fhi) ) exit find_xi

                  ! Get normal at interpolation point                  
                  call amrex_eb_normal_levelset(x_i, x0, n_refine, phi, phlo, phhi, dx, norm_i)

                  x_i = x_i + maxval(dx) * norm_i

                  iter = iter + 1

                  if ( iter > max_iter ) &
                   call amrex_abort("reconstruct_velocity: cannot find mirror point")
                  
               end do find_xi

               ! Get phi at interpolation point
               call amrex_eb_interp_levelset(x_i, x0, n_refine, phi, phlo, phhi, dx, phi_i)
               
               ! Compute interpolated velocity at x_i
               vel_i = trilinear_interp(vel_in, vilo, vihi, 3, x_i, x0, dx)

               ! Since interpolation point is only slightly shifted with respect to
               ! the mirror point, we approximate vel at mirror point with vel_i and
               ! then use linear interpolation between x_m and x_c
               vel_out(i,j,k,:) = vel_i * phi_cc / phi_i 

            end if

         end do
      end do
   end do


end subroutine reconstruct_velocity
