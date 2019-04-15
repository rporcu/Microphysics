module mfix_particle_pic_module

   use iso_c_binding
   use amrex_fort_module,               only: amrex_real, amrex_particle_real
   use param,                           only: one, zero, half
   use amrex_ebcellflag_module,         only: is_covered_cell

   implicit none

   private

   public :: mfix_deposit_cic_eb
   public :: mfix_multi_deposit_cic_eb

contains

   subroutine mfix_deposit_cic_eb ( particles, ns, np,                  &
        &                           vol   , vlo, vhi,                   &
        &                           vratio, slo, shi,                   &
        &                           flags,  flo, fhi,                   &
        &                           plo, dx, particle_comp ) bind(C)

      real(amrex_particle_real), intent(in   )  :: particles(ns,np)
      integer, value,            intent(in   )  :: ns, np
      integer,                   intent(in   )  :: vlo(3), vhi(3)
      integer,                   intent(in   )  :: slo(3), shi(3)
      integer,                   intent(in   )  :: flo(3), fhi(3)

      real(amrex_real),          intent(in   )  ::           &
           &  vratio(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3))

      real(amrex_real),          intent(inout)  ::           &
           &  vol(vlo(1):vhi(1), vlo(2):vhi(2), vlo(3):vhi(3))

      integer(c_int),             intent(in   ) ::  &
           & flags(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

      real(amrex_real),          intent(in   )  :: plo(3), dx(3)
      integer,                   intent(in   )  :: particle_comp


      ! Local variables
      integer                       :: i, j, k, ii, jj, kk, n
      real(amrex_real)              :: wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
      real(amrex_real)              :: lx, ly, lz, pvol
      real(amrex_real)              :: inv_dx(3)
      real(amrex_real)              :: weights(-1:0,-1:0,-1:0)
      real(amrex_real)              :: regular_cell_volume, this_cell_volume

      inv_dx = one/dx

      regular_cell_volume = dx(1)*dx(2)*dx(3)

      do n = 1, np

         lx = (particles(1, n) - plo(1))*inv_dx(1) + half
         ly = (particles(2, n) - plo(2))*inv_dx(2) + half
         lz = (particles(3, n) - plo(3))*inv_dx(3) + half

         i = floor(lx)
         j = floor(ly)
         k = floor(lz)

         ! Compute deposition weights: use a combination of
         ! trilinear weighting and volume ratio weighting.
         wx_hi = lx - i
         wy_hi = ly - j
         wz_hi = lz - k

         wx_lo = one - wx_hi
         wy_lo = one - wy_hi
         wz_lo = one - wz_hi

         weights(-1,-1,-1) = vratio(i-1, j-1, k-1) * wx_lo * wy_lo * wz_lo
         weights(-1,-1, 0) = vratio(i-1, j-1, k  ) * wx_lo * wy_lo * wz_hi
         weights(-1, 0,-1) = vratio(i-1, j,   k-1) * wx_lo * wy_hi * wz_lo
         weights(-1, 0, 0) = vratio(i-1, j,   k  ) * wx_lo * wy_hi * wz_hi
         weights( 0,-1,-1) = vratio(i,   j-1, k-1) * wx_hi * wy_lo * wz_lo
         weights( 0,-1, 0) = vratio(i,   j-1, k  ) * wx_hi * wy_lo * wz_hi
         weights( 0, 0,-1) = vratio(i,   j,   k-1) * wx_hi * wy_hi * wz_lo
         weights( 0, 0, 0) = vratio(i,   j,   k  ) * wx_hi * wy_hi * wz_hi

         ! Normalize so that weights sums up to one
         weights = weights / sum(weights)

         pvol = particles(particle_comp,n)

         ! Perform deposition
         do kk = -1, 0
            do jj = -1, 0
               do ii = -1, 0
                  if (is_covered_cell(flags(i+ii, j+jj, k+kk))) cycle
                  this_cell_volume = vratio(i+ii, j+jj, k+kk) * regular_cell_volume
                  vol(i+ii, j+jj, k+kk) = vol(i+ii, j+jj, k+kk) + &
                       &                     weights(ii,jj,kk) * pvol / this_cell_volume
               end do
            end do
         end do

      end do

   end subroutine mfix_deposit_cic_eb




   subroutine mfix_multi_deposit_cic_eb (particles, ns, np, &
        mf_x, mf_u, lo, hi, &
        vratio, vlo, vhi,   &
        flags,  flo, fhi,   &
        plo, dx, beta_comp, vel_comp) bind(c)

      real(amrex_particle_real), intent(in   )  :: particles(ns,np)
      integer, value,            intent(in   )  :: ns, np
      integer,                   intent(in   )  ::  lo(3),  hi(3)
      integer,                   intent(in   )  :: vlo(3), vhi(3)
      integer,                   intent(in   )  :: flo(3), fhi(3)

      real(amrex_real),          intent(in   )  ::           &
           &  vratio(vlo(1):vhi(1), vlo(2):vhi(2), vlo(3):vhi(3))

      integer(c_int),             intent(in   ) ::  &
           & flags(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

      real(amrex_real),          intent(inout)  ::           &
           &  mf_x(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)),   &
           &  mf_u(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),3)

      real(amrex_real),          intent(in   )  :: plo(3), dx(3)
      integer,                   intent(in   )  :: beta_comp, vel_comp

      ! Local variables
      integer           :: i, j, k, ii, jj, kk, n, idir
      real(amrex_real)  :: wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
      real(amrex_real)  :: lx, ly, lz, pbeta, pvel(3), weights(-1:0,-1:0,-1:0)
      real(amrex_real)  :: inv_dx(3)

      real(amrex_real)  :: weight_vol, regular_cell_volume, this_cell_volume

      regular_cell_volume = dx(1)*dx(2)*dx(3)

      inv_dx = one/dx

      do n = 1, np

         lx = (particles(1, n) - plo(1))*inv_dx(1) + half
         ly = (particles(2, n) - plo(2))*inv_dx(2) + half
         lz = (particles(3, n) - plo(3))*inv_dx(3) + half

         i = floor(lx)
         j = floor(ly)
         k = floor(lz)

         wx_hi = lx - i
         wy_hi = ly - j
         wz_hi = lz - k

         wx_lo = one - wx_hi
         wy_lo = one - wy_hi
         wz_lo = one - wz_hi

         ! Compute deposition weights: use a combination of
         ! trilinear weighting and volume ratio weighting.
         weights(-1,-1,-1) = vratio(i-1, j-1, k-1) * wx_lo * wy_lo * wz_lo
         weights(-1,-1, 0) = vratio(i-1, j-1, k  ) * wx_lo * wy_lo * wz_hi
         weights(-1, 0,-1) = vratio(i-1, j,   k-1) * wx_lo * wy_hi * wz_lo
         weights(-1, 0, 0) = vratio(i-1, j,   k  ) * wx_lo * wy_hi * wz_hi
         weights( 0,-1,-1) = vratio(i,   j-1, k-1) * wx_hi * wy_lo * wz_lo
         weights( 0,-1, 0) = vratio(i,   j-1, k  ) * wx_hi * wy_lo * wz_hi
         weights( 0, 0,-1) = vratio(i,   j,   k-1) * wx_hi * wy_hi * wz_lo
         weights( 0, 0, 0) = vratio(i,   j,   k  ) * wx_hi * wy_hi * wz_hi

         ! Normalize so that weights sums up to one
         weights = weights / sum(weights)

         pbeta   = particles(beta_comp,n)
         pvel(1) = particles(vel_comp  ,n) * pbeta
         pvel(2) = particles(vel_comp+1,n) * pbeta
         pvel(3) = particles(vel_comp+2,n) * pbeta


         ! The projection method uses drag to update u, not (cell_vol * u),
         ! so we must divide by the cell volume here so that later when we
         ! divide by bulk density it all works.

         ! Perform deposition
         do kk = -1, 0
            do jj = -1, 0
               do ii = -1, 0

                  if (is_covered_cell(flags(i+ii, j+jj, k+kk))) cycle

                  this_cell_volume = vratio(i+ii, j+jj, k+kk) * regular_cell_volume

                  weight_vol = weights(ii,jj,kk) / this_cell_volume


                  mf_x(i+ii, j+jj, k+kk) = mf_x(i+ii, j+jj, k+kk) + &
                       &                     weight_vol * pbeta

                  do idir = 1, 3
                     mf_u(i+ii, j+jj, k+kk, idir) = mf_u(i+ii, j+jj, k+kk, idir) + &
                          &                     weight_vol * pvel(idir)
                  end do

               end do
            end do
         end do

      end do

   end subroutine mfix_multi_deposit_cic_eb

end module mfix_particle_pic_module
