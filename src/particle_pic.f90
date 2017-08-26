module mfix_particle_pic_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real, amrex_particle_real

  implicit none

  private

  public :: mfix_deposit_cic, mfix_interpolate_cic

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

  subroutine mfix_multi_deposit_cic(particles, ns, np, nc, mf, lo, hi, plo, dx, beta_comp, vel_comp) &
       bind(c,name='mfix_multi_deposit_cic')
    integer, value                :: ns, np, nc
    real(amrex_particle_real)     :: particles(ns,np)
    integer                       :: lo(3)
    integer                       :: hi(3)
    real(amrex_real)              :: mf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),nc)
    real(amrex_real)              :: plo(3)
    real(amrex_real)              :: dx(3)
    integer                       :: beta_comp, vel_comp

    integer i, j, k, n, idim
    integer ii, jj, kk
    real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
    real(amrex_real) wu_lo, wv_lo, ww_lo, wu_hi, wv_hi, ww_hi
    real(amrex_real) lx, ly, lz, pbeta, pvel(3)
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

       pbeta   = particles(beta_comp,n) * oovol
       pvel(1) = particles(vel_comp  ,n)
       pvel(2) = particles(vel_comp+1,n)
       pvel(3) = particles(vel_comp+2,n)

       mf(i-1, j-1, k-1, 1) = mf(i-1, j-1, k-1, 1) + wx_lo*wy_lo*wz_lo*pbeta
       mf(i-1, j-1, k  , 1) = mf(i-1, j-1, k,   1) + wx_lo*wy_lo*wz_hi*pbeta
       mf(i-1, j,   k-1, 1) = mf(i-1, j,   k-1, 1) + wx_lo*wy_hi*wz_lo*pbeta
       mf(i-1, j,   k  , 1) = mf(i-1, j,   k,   1) + wx_lo*wy_hi*wz_hi*pbeta
       mf(i,   j-1, k-1, 1) = mf(i,   j-1, k-1, 1) + wx_hi*wy_lo*wz_lo*pbeta
       mf(i,   j-1, k  , 1) = mf(i,   j-1, k,   1) + wx_hi*wy_lo*wz_hi*pbeta
       mf(i,   j,   k-1, 1) = mf(i,   j,   k-1, 1) + wx_hi*wy_hi*wz_lo*pbeta
       mf(i,   j,   k  , 1) = mf(i,   j,   k,   1) + wx_hi*wy_hi*wz_hi*pbeta


       lx = (particles(1, n) - plo(1))*inv_dx(1)
       ii = floor(lx)
       wu_hi = lx - ii
       wu_lo = 1.0d0 - wu_hi

       idim = 1
       mf(ii,   j-1, k-1, idim+1) = mf(ii,   j-1, k-1, idim+1) + wu_lo*wy_lo*wz_lo*pbeta*pvel(idim)
       mf(ii,   j-1, k  , idim+1) = mf(ii,   j-1, k,   idim+1) + wu_lo*wy_lo*wz_hi*pbeta*pvel(idim)
       mf(ii,   j,   k-1, idim+1) = mf(ii,   j,   k-1, idim+1) + wu_lo*wy_hi*wz_lo*pbeta*pvel(idim)
       mf(ii,   j,   k  , idim+1) = mf(ii,   j,   k,   idim+1) + wu_lo*wy_hi*wz_hi*pbeta*pvel(idim)
       mf(ii+1, j-1, k-1, idim+1) = mf(ii+1, j-1, k-1, idim+1) + wu_hi*wy_lo*wz_lo*pbeta*pvel(idim)
       mf(ii+1, j-1, k  , idim+1) = mf(ii+1, j-1, k,   idim+1) + wu_hi*wy_lo*wz_hi*pbeta*pvel(idim)
       mf(ii+1, j,   k-1, idim+1) = mf(ii+1, j,   k-1, idim+1) + wu_hi*wy_hi*wz_lo*pbeta*pvel(idim)
       mf(ii+1, j,   k  , idim+1) = mf(ii+1, j,   k,   idim+1) + wu_hi*wy_hi*wz_hi*pbeta*pvel(idim)


       ly = (particles(2, n) - plo(2))*inv_dx(2)
       jj = floor(ly)
       wv_hi = ly - jj
       wv_lo = 1.0d0 - wv_hi

       idim = 2
       mf(i-1, jj,   k-1, idim+1) = mf(i-1, jj,   k-1, idim+1) + wx_lo*wv_lo*wz_lo*pbeta*pvel(idim)
       mf(i-1, jj,   k  , idim+1) = mf(i-1, jj,   k,   idim+1) + wx_lo*wv_lo*wz_hi*pbeta*pvel(idim)
       mf(i-1, jj+1, k-1, idim+1) = mf(i-1, jj+1, k-1, idim+1) + wx_lo*wv_hi*wz_lo*pbeta*pvel(idim)
       mf(i-1, jj+1, k  , idim+1) = mf(i-1, jj+1, k,   idim+1) + wx_lo*wv_hi*wz_hi*pbeta*pvel(idim)
       mf(i,   jj,   k-1, idim+1) = mf(i,   jj,   k-1, idim+1) + wx_hi*wv_lo*wz_lo*pbeta*pvel(idim)
       mf(i,   jj,   k  , idim+1) = mf(i,   jj,   k,   idim+1) + wx_hi*wv_lo*wz_hi*pbeta*pvel(idim)
       mf(i,   jj+1, k-1, idim+1) = mf(i,   jj+1, k-1, idim+1) + wx_hi*wv_hi*wz_lo*pbeta*pvel(idim)
       mf(i,   jj+1, k  , idim+1) = mf(i,   jj+1, k,   idim+1) + wx_hi*wv_hi*wz_hi*pbeta*pvel(idim)


       lz = (particles(3, n) - plo(3))*inv_dx(3)
       kk = floor(lz)
       ww_hi = lz - kk
       ww_lo = 1.0d0 - ww_hi

       idim = 3
       mf(i-1, j-1, kk,   idim+1) = mf(i-1, j-1, kk,   idim+1) + wx_lo*wy_lo*ww_lo*pbeta*pvel(idim)
       mf(i-1, j-1, kk+1, idim+1) = mf(i-1, j-1, kk+1, idim+1) + wx_lo*wy_lo*ww_hi*pbeta*pvel(idim)
       mf(i-1, j,   kk,   idim+1) = mf(i-1, j,   kk,   idim+1) + wx_lo*wy_hi*ww_lo*pbeta*pvel(idim)
       mf(i-1, j,   kk+1, idim+1) = mf(i-1, j,   kk+1, idim+1) + wx_lo*wy_hi*ww_hi*pbeta*pvel(idim)
       mf(i,   j-1, kk,   idim+1) = mf(i,   j-1, kk,   idim+1) + wx_hi*wy_lo*ww_lo*pbeta*pvel(idim)
       mf(i,   j-1, kk+1, idim+1) = mf(i,   j-1, kk+1, idim+1) + wx_hi*wy_lo*ww_hi*pbeta*pvel(idim)
       mf(i,   j,   kk,   idim+1) = mf(i,   j,   kk,   idim+1) + wx_hi*wy_hi*ww_lo*pbeta*pvel(idim)
       mf(i,   j,   kk+1, idim+1) = mf(i,   j,   kk+1, idim+1) + wx_hi*wy_hi*ww_hi*pbeta*pvel(idim)

    end do

  end subroutine mfix_multi_deposit_cic

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
