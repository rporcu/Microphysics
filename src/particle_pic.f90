module mfix_particle_pic_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real, amrex_particle_real

  implicit none

  private

  public :: mfix_deposit_cic, mfix_interpolate_cic

contains

  subroutine mfix_deposit_cic(particles, ns, np, nc, vol, lo, hi, plo, dx) &
       bind(c,name='mfix_deposit_cic')
    integer, value                :: ns, np, nc
    real(amrex_particle_real)     :: particles(ns,np)
    integer                       :: lo(3)
    integer                       :: hi(3)
    real(amrex_real)              :: vol(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),nc)
    real(amrex_real)              :: plo(3)
    real(amrex_real)              :: dx(3)

    integer i, j, k, n, comp
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

       pvol = particles(5,n) * oovol

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

  subroutine mfix_interpolate_cic(particles, ns, np, acc, lo, hi, ncomp, plo, dx) &
       bind(c,name='mfix_interpolate_cic')
    integer, value                :: ns, np, ncomp
    real(amrex_particle_real)     :: particles(ns,np)
    integer                       :: lo(3)
    integer                       :: hi(3)
    real(amrex_real)              :: acc(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), ncomp)
    real(amrex_real)              :: plo(3)
    real(amrex_real)              :: dx(3)
    real(amrex_real)              :: acceleration(ncomp)

    integer i, j, k, n, nc
    real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
    real(amrex_real) lx, ly, lz
    real(amrex_real) inv_dx(3)
    inv_dx = 1.0d0/dx

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

       do nc = 1, ncomp
          acceleration(nc) = wx_lo*wy_lo*wz_lo*acc(i-1, j-1, k-1, nc) + &
                             wx_lo*wy_lo*wz_hi*acc(i-1, j-1, k  , nc) + &
                             wx_lo*wy_hi*wz_lo*acc(i-1, j,   k-1, nc) + &
                             wx_lo*wy_hi*wz_hi*acc(i-1, j,   k  , nc) + &
                             wx_hi*wy_lo*wz_lo*acc(i,   j-1, k-1, nc) + &
                             wx_hi*wy_lo*wz_hi*acc(i,   j-1, k  , nc) + &
                             wx_hi*wy_hi*wz_lo*acc(i,   j,   k-1, nc) + &
                             wx_hi*wy_hi*wz_hi*acc(i,   j,   k  , nc)
       
          if (abs(acceleration(nc) - 5.d0) .ge. 1.0d-9) then
             print *, particles(1, n), particles(2, n), particles(3, n)
          end if

       end do
    end do

  end subroutine mfix_interpolate_cic

end module mfix_particle_pic_module
