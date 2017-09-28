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

  subroutine mfix_multi_deposit_cic(particles, ns, np, mf_x, mf_y, mf_z, mf_u, mf_v, mf_w, &
                                    lo_x, hi_x, lo_y, hi_y, lo_z, hi_z, plo, dx, beta_comp, vel_comp) &
       bind(c,name='mfix_multi_deposit_cic')
    integer, value                :: ns, np
    real(amrex_particle_real)     :: particles(ns,np)
    integer                       :: lo_x(3), lo_y(3), lo_z(3)
    integer                       :: hi_x(3), hi_y(3), hi_z(3)
    real(amrex_real)              :: mf_x(lo_x(1):hi_x(1), lo_x(2):hi_x(2), lo_x(3):hi_x(3))
    real(amrex_real)              :: mf_y(lo_y(1):hi_y(1), lo_y(2):hi_y(2), lo_y(3):hi_y(3))
    real(amrex_real)              :: mf_z(lo_z(1):hi_z(1), lo_z(2):hi_z(2), lo_z(3):hi_z(3))
    real(amrex_real)              :: mf_u(lo_x(1):hi_x(1), lo_x(2):hi_x(2), lo_x(3):hi_x(3))
    real(amrex_real)              :: mf_v(lo_y(1):hi_y(1), lo_y(2):hi_y(2), lo_y(3):hi_y(3))
    real(amrex_real)              :: mf_w(lo_z(1):hi_z(1), lo_z(2):hi_z(2), lo_z(3):hi_z(3))
    real(amrex_real)              :: plo(3)
    real(amrex_real)              :: dx(3)
    integer                       :: beta_comp, vel_comp

    integer i, j, k, n
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

       pbeta   = particles(beta_comp,n)
       pvel(1) = particles(vel_comp  ,n) * pbeta
       pvel(2) = particles(vel_comp+1,n) * pbeta
       pvel(3) = particles(vel_comp+2,n) * pbeta

       lx = (particles(1, n) - plo(1))*inv_dx(1)
       ii = floor(lx)
       wu_hi = lx - ii
       wu_lo = 1.0d0 - wu_hi

       mf_x(ii,   j-1, k-1) = mf_x(ii,   j-1, k-1) + wu_lo*wy_lo*wz_lo*pbeta
       mf_x(ii,   j-1, k  ) = mf_x(ii,   j-1, k  ) + wu_lo*wy_lo*wz_hi*pbeta
       mf_x(ii,   j,   k-1) = mf_x(ii,   j,   k-1) + wu_lo*wy_hi*wz_lo*pbeta
       mf_x(ii,   j,   k  ) = mf_x(ii,   j,   k  ) + wu_lo*wy_hi*wz_hi*pbeta
       mf_x(ii+1, j-1, k-1) = mf_x(ii+1, j-1, k-1) + wu_hi*wy_lo*wz_lo*pbeta
       mf_x(ii+1, j-1, k  ) = mf_x(ii+1, j-1, k  ) + wu_hi*wy_lo*wz_hi*pbeta
       mf_x(ii+1, j,   k-1) = mf_x(ii+1, j,   k-1) + wu_hi*wy_hi*wz_lo*pbeta
       mf_x(ii+1, j,   k  ) = mf_x(ii+1, j,   k  ) + wu_hi*wy_hi*wz_hi*pbeta

       mf_u(ii,   j-1, k-1) = mf_u(ii,   j-1, k-1) + wu_lo*wy_lo*wz_lo*pvel(1)
       mf_u(ii,   j-1, k  ) = mf_u(ii,   j-1, k  ) + wu_lo*wy_lo*wz_hi*pvel(1)
       mf_u(ii,   j,   k-1) = mf_u(ii,   j,   k-1) + wu_lo*wy_hi*wz_lo*pvel(1)
       mf_u(ii,   j,   k  ) = mf_u(ii,   j,   k  ) + wu_lo*wy_hi*wz_hi*pvel(1)
       mf_u(ii+1, j-1, k-1) = mf_u(ii+1, j-1, k-1) + wu_hi*wy_lo*wz_lo*pvel(1)
       mf_u(ii+1, j-1, k  ) = mf_u(ii+1, j-1, k  ) + wu_hi*wy_lo*wz_hi*pvel(1)
       mf_u(ii+1, j,   k-1) = mf_u(ii+1, j,   k-1) + wu_hi*wy_hi*wz_lo*pvel(1)
       mf_u(ii+1, j,   k  ) = mf_u(ii+1, j,   k  ) + wu_hi*wy_hi*wz_hi*pvel(1)

       ly = (particles(2, n) - plo(2))*inv_dx(2)
       jj = floor(ly)
       wv_hi = ly - jj
       wv_lo = 1.0d0 - wv_hi

       mf_y(i-1, jj,   k-1) = mf_y(i-1, jj,   k-1) + wx_lo*wv_lo*wz_lo*pbeta
       mf_y(i-1, jj,   k  ) = mf_y(i-1, jj,   k  ) + wx_lo*wv_lo*wz_hi*pbeta
       mf_y(i-1, jj+1, k-1) = mf_y(i-1, jj+1, k-1) + wx_lo*wv_hi*wz_lo*pbeta
       mf_y(i-1, jj+1, k  ) = mf_y(i-1, jj+1, k  ) + wx_lo*wv_hi*wz_hi*pbeta
       mf_y(i,   jj,   k-1) = mf_y(i,   jj,   k-1) + wx_hi*wv_lo*wz_lo*pbeta
       mf_y(i,   jj,   k  ) = mf_y(i,   jj,   k  ) + wx_hi*wv_lo*wz_hi*pbeta
       mf_y(i,   jj+1, k-1) = mf_y(i,   jj+1, k-1) + wx_hi*wv_hi*wz_lo*pbeta
       mf_y(i,   jj+1, k  ) = mf_y(i,   jj+1, k  ) + wx_hi*wv_hi*wz_hi*pbeta

       mf_v(i-1, jj,   k-1) = mf_v(i-1, jj,   k-1) + wx_lo*wv_lo*wz_lo*pvel(2)
       mf_v(i-1, jj,   k  ) = mf_v(i-1, jj,   k  ) + wx_lo*wv_lo*wz_hi*pvel(2)
       mf_v(i-1, jj+1, k-1) = mf_v(i-1, jj+1, k-1) + wx_lo*wv_hi*wz_lo*pvel(2)
       mf_v(i-1, jj+1, k  ) = mf_v(i-1, jj+1, k  ) + wx_lo*wv_hi*wz_hi*pvel(2)
       mf_v(i,   jj,   k-1) = mf_v(i,   jj,   k-1) + wx_hi*wv_lo*wz_lo*pvel(2)
       mf_v(i,   jj,   k  ) = mf_v(i,   jj,   k  ) + wx_hi*wv_lo*wz_hi*pvel(2)
       mf_v(i,   jj+1, k-1) = mf_v(i,   jj+1, k-1) + wx_hi*wv_hi*wz_lo*pvel(2)
       mf_v(i,   jj+1, k  ) = mf_v(i,   jj+1, k  ) + wx_hi*wv_hi*wz_hi*pvel(2)

       lz = (particles(3, n) - plo(3))*inv_dx(3)
       kk = floor(lz)
       ww_hi = lz - kk
       ww_lo = 1.0d0 - ww_hi

       mf_z(i-1, j-1, kk  ) = mf_z(i-1, j-1, kk  ) + wx_lo*wy_lo*ww_lo*pbeta
       mf_z(i-1, j-1, kk+1) = mf_z(i-1, j-1, kk+1) + wx_lo*wy_lo*ww_hi*pbeta
       mf_z(i-1, j,   kk  ) = mf_z(i-1, j,   kk  ) + wx_lo*wy_hi*ww_lo*pbeta
       mf_z(i-1, j,   kk+1) = mf_z(i-1, j,   kk+1) + wx_lo*wy_hi*ww_hi*pbeta
       mf_z(i,   j-1, kk  ) = mf_z(i,   j-1, kk  ) + wx_hi*wy_lo*ww_lo*pbeta
       mf_z(i,   j-1, kk+1) = mf_z(i,   j-1, kk+1) + wx_hi*wy_lo*ww_hi*pbeta
       mf_z(i,   j,   kk  ) = mf_z(i,   j,   kk  ) + wx_hi*wy_hi*ww_lo*pbeta
       mf_z(i,   j,   kk+1) = mf_z(i,   j,   kk+1) + wx_hi*wy_hi*ww_hi*pbeta

       mf_w(i-1, j-1, kk  ) = mf_w(i-1, j-1, kk  ) + wx_lo*wy_lo*ww_lo*pvel(3)
       mf_w(i-1, j-1, kk+1) = mf_w(i-1, j-1, kk+1) + wx_lo*wy_lo*ww_hi*pvel(3)
       mf_w(i-1, j,   kk  ) = mf_w(i-1, j,   kk  ) + wx_lo*wy_hi*ww_lo*pvel(3)
       mf_w(i-1, j,   kk+1) = mf_w(i-1, j,   kk+1) + wx_lo*wy_hi*ww_hi*pvel(3)
       mf_w(i,   j-1, kk  ) = mf_w(i,   j-1, kk  ) + wx_hi*wy_lo*ww_lo*pvel(3)
       mf_w(i,   j-1, kk+1) = mf_w(i,   j-1, kk+1) + wx_hi*wy_lo*ww_hi*pvel(3)
       mf_w(i,   j,   kk  ) = mf_w(i,   j,   kk  ) + wx_hi*wy_hi*ww_lo*pvel(3)
       mf_w(i,   j,   kk+1) = mf_w(i,   j,   kk+1) + wx_hi*wy_hi*ww_hi*pvel(3)

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
