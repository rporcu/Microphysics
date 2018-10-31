module mfix_particle_pic_module

   use iso_c_binding
   use amrex_fort_module, only : amrex_real, amrex_particle_real
   use amrex_paralleldescriptor_module, only : amrex_pd_ioprocessor

   implicit none

   private

   public :: mfix_deposit_cic

contains

   subroutine mfix_deposit_cic(particles, ns, np, nc, vol, lo, hi, plo, dx, particle_comp) &
        bind(c,name='mfix_deposit_cic')
      integer, value                :: ns, np, nc
      real(amrex_particle_real)     :: particles(ns,np)
      integer                       :: lo(3)
      integer                       :: hi(3)
      real(amrex_real)              :: vol(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),nc)
      real(amrex_real)              :: plo(3)
      real(amrex_real)              :: dx(3)
      integer                       :: particle_comp

      integer i, j, k, n
      real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
      real(amrex_real) lx, ly, lz, pvol
      real(amrex_real) inv_dx(3), oovol
      inv_dx = 1.0d0/dx

      oovol = 1.0d0/(dx(1)*dx(2)*dx(3))

      do n = 1, np

         lx = (particles(1, n) - plo(1))*inv_dx(1) + 0.5d0
         ly = (particles(2, n) - plo(2))*inv_dx(2) + 0.5d0
         lz = (particles(3, n) - plo(3))*inv_dx(3) + 0.5d0

         i = floor(lx)
         j = floor(ly)
         k = floor(lz)

         wx_hi = lx - i
         wy_hi = ly - j
         wz_hi = lz - k

         wx_lo = 1.0d0 - wx_hi
         wy_lo = 1.0d0 - wy_hi
         wz_lo = 1.0d0 - wz_hi

         pvol = particles(particle_comp,n) * oovol

         vol(i-1, j-1, k-1, 1) = vol(i-1, j-1, k-1, 1) + wx_lo*wy_lo*wz_lo*pvol
         vol(i-1, j-1, k  , 1) = vol(i-1, j-1, k  , 1) + wx_lo*wy_lo*wz_hi*pvol
         vol(i-1, j,   k-1, 1) = vol(i-1, j,   k-1, 1) + wx_lo*wy_hi*wz_lo*pvol
         vol(i-1, j,   k  , 1) = vol(i-1, j,   k,   1) + wx_lo*wy_hi*wz_hi*pvol
         vol(i,   j-1, k-1, 1) = vol(i,   j-1, k-1, 1) + wx_hi*wy_lo*wz_lo*pvol
         vol(i,   j-1, k  , 1) = vol(i,   j-1, k  , 1) + wx_hi*wy_lo*wz_hi*pvol
         vol(i,   j,   k-1, 1) = vol(i,   j,   k-1, 1) + wx_hi*wy_hi*wz_lo*pvol
         vol(i,   j,   k  , 1) = vol(i,   j,   k  , 1) + wx_hi*wy_hi*wz_hi*pvol

      end do

   end subroutine mfix_deposit_cic


   subroutine mfix_deposit_cic_eb ( particles, ns, np, nc, vol, lo, hi, vratio, slo, shi, &
        &                           plo, dx, particle_comp ) bind(C)

      use param, only: one, zero, half
      
      integer, value                :: ns, np, nc
      real(amrex_particle_real)     :: particles(ns,np)
      integer                       :: lo(3), slo(3)
      integer                       :: hi(3), shi(3)
      real(amrex_real)              :: vol(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),nc)
      real(amrex_real)              :: vratio(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3))
      real(amrex_real)              :: plo(3)
      real(amrex_real)              :: dx(3)
      integer                       :: particle_comp
      integer                       :: i, j, k, ii, jj, kk, n
      real(amrex_real)              :: wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
      real(amrex_real)              :: lx, ly, lz, pvol
      real(amrex_real)              :: inv_dx(3), oovol
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


         do kk = -1, 0
            do jj = -1, 0
               do ii = -1, 0
                  this_cell_volume = vratio(i+ii, j+jj, k+kk) * regular_cell_volume
                  vol(i+ii, j+jj, k+kk, 1) = vol(i+ii, j+jj, k+kk, 1) + &
                       &                     weights(ii,jj,kk) * pvol / this_cell_volume
               end do
            end do
         end do

      end do

   end subroutine mfix_deposit_cic_eb


   


   subroutine mfix_multi_deposit_cic(particles, ns, np, mf_x, mf_u, &
        lo, hi, plo, dx, beta_comp, vel_comp) &
        bind(c,name='mfix_multi_deposit_cic')


      integer, value                :: ns, np
      real(amrex_particle_real)     :: particles(ns,np)
      integer                       :: lo(3), hi(3)
      real(amrex_real)              :: mf_x(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
      real(amrex_real)              :: mf_u(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),3)
      real(amrex_real)              :: plo(3)
      real(amrex_real)              :: dx(3)
      integer                       :: beta_comp, vel_comp

      integer i, j, k, n, idir
      real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
      real(amrex_real) lx, ly, lz, pbeta, pvel(3)
      real(amrex_real) inv_dx(3)

      logical, parameter :: no_interpolation = .false.

      if(no_interpolation .and. amrex_pd_ioprocessor()) &
           write(*,*) 'WARNING: No interpolation of drag force grid deposition'

      inv_dx = 1.0d0/dx

      if (no_interpolation) then

         do n = 1, np

            lx = (particles(1, n) - plo(1))*inv_dx(1)
            ly = (particles(2, n) - plo(2))*inv_dx(2)
            lz = (particles(3, n) - plo(3))*inv_dx(3)

            i = floor(lx)
            j = floor(ly)
            k = floor(lz)

            pbeta   = particles(beta_comp,n)

            mf_x(i,j,k)   = mf_x(i,j,k)   + pbeta
            mf_u(i,j,k,1) = mf_u(i,j,k,1) + pbeta * particles(vel_comp  ,n)
            mf_u(i,j,k,2) = mf_u(i,j,k,2) + pbeta * particles(vel_comp+1,n)
            mf_u(i,j,k,3) = mf_u(i,j,k,3) + pbeta * particles(vel_comp+2,n)

         end do

      else

         do n = 1, np

            lx = (particles(1, n) - plo(1))*inv_dx(1) + 0.5d0
            ly = (particles(2, n) - plo(2))*inv_dx(2) + 0.5d0
            lz = (particles(3, n) - plo(3))*inv_dx(3) + 0.5d0

            i = floor(lx)
            j = floor(ly)
            k = floor(lz)

            wx_hi = lx - i
            wy_hi = ly - j
            wz_hi = lz - k

            wx_lo = 1.0d0 - wx_hi
            wy_lo = 1.0d0 - wy_hi
            wz_lo = 1.0d0 - wz_hi

            pbeta   = particles(beta_comp,n)
            pvel(1) = particles(vel_comp  ,n) * pbeta
            pvel(2) = particles(vel_comp+1,n) * pbeta
            pvel(3) = particles(vel_comp+2,n) * pbeta

            mf_x(i-1, j-1, k-1) = mf_x(i-1, j-1, k-1) + wx_lo*wy_lo*wz_lo*pbeta
            mf_x(i-1, j-1, k  ) = mf_x(i-1, j-1, k  ) + wx_lo*wy_lo*wz_hi*pbeta
            mf_x(i-1, j,   k-1) = mf_x(i-1, j,   k-1) + wx_lo*wy_hi*wz_lo*pbeta
            mf_x(i-1, j,   k  ) = mf_x(i-1, j,   k  ) + wx_lo*wy_hi*wz_hi*pbeta
            mf_x(i  , j-1, k-1) = mf_x(i  , j-1, k-1) + wx_hi*wy_lo*wz_lo*pbeta
            mf_x(i  , j-1, k  ) = mf_x(i  , j-1, k  ) + wx_hi*wy_lo*wz_hi*pbeta
            mf_x(i  , j,   k-1) = mf_x(i  , j,   k-1) + wx_hi*wy_hi*wz_lo*pbeta
            mf_x(i  , j,   k  ) = mf_x(i  , j,   k  ) + wx_hi*wy_hi*wz_hi*pbeta

            do idir = 1, 3
               mf_u(i-1, j-1, k-1,idir) = mf_u(i-1, j-1, k-1,idir) + wx_lo*wy_lo*wz_lo*pvel(idir)
               mf_u(i-1, j-1, k  ,idir) = mf_u(i-1, j-1, k  ,idir) + wx_lo*wy_lo*wz_hi*pvel(idir)
               mf_u(i-1, j,   k-1,idir) = mf_u(i-1, j,   k-1,idir) + wx_lo*wy_hi*wz_lo*pvel(idir)
               mf_u(i-1, j,   k  ,idir) = mf_u(i-1, j,   k  ,idir) + wx_lo*wy_hi*wz_hi*pvel(idir)
               mf_u(i  , j-1, k-1,idir) = mf_u(i  , j-1, k-1,idir) + wx_hi*wy_lo*wz_lo*pvel(idir)
               mf_u(i  , j-1, k  ,idir) = mf_u(i  , j-1, k  ,idir) + wx_hi*wy_lo*wz_hi*pvel(idir)
               mf_u(i  , j,   k-1,idir) = mf_u(i  , j,   k-1,idir) + wx_hi*wy_hi*wz_lo*pvel(idir)
               mf_u(i  , j,   k  ,idir) = mf_u(i  , j,   k  ,idir) + wx_hi*wy_hi*wz_hi*pvel(idir)
            end do

         end do
      endif


   end subroutine mfix_multi_deposit_cic

end module mfix_particle_pic_module
