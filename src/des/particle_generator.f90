module par_gen_module

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int

  implicit none

  real(c_real),   allocatable :: rdata(:,:)
  integer(c_int), allocatable :: idata(:,:)

  !< Position............... 1,2,3
  !< Radius................. 4
  !< Density................ 5
  !< Linear velocity........ 6,7,8
  integer, parameter :: nr =  8

  !< Type................... 1
  integer, parameter :: ni =  1

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: GENERATE_PARTICLE_CONFIG                                !
!                                                                      !
!  Purpose: Generate particle configuration based on maximum particle  !
!           radius and filling from top to bottom within specified     !
!           bounds                                                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine particle_generator(pc, lo, hi, dx, dy, dz) &
  bind(C, name="mfix_particle_generator")

  use ic, only: dim_ic, ic_defined
  use ic, only: ic_ep_s, ic_ep_g

  use ic, only: ic_x_e, ic_y_n, ic_z_t
  use ic, only: ic_x_w, ic_y_s, ic_z_b

  use ic, only: ic_u_s, ic_v_s, ic_w_s

  use ic, only: ic_dp_dist, ic_ro_s_dist
  use ic, only: ic_dp_mean, ic_ro_s_mean
  use ic, only: ic_dp_std,  ic_ro_s_std
  use ic, only: ic_dp_min,  ic_ro_s_min
  use ic, only: ic_dp_max,  ic_ro_s_max

  use param, only: undefined, is_defined

  use calc_cell_module, only: calc_cell_ic
  use discretelement, only: particle_types
  use constant, only: pi

  implicit none

  integer(c_int), intent(inout) :: pc
  integer(c_int), intent(in   ) :: lo(3),hi(3)
  real(c_real),   intent(in   ) :: dx, dy, dz

  ! local index for initial condition
  integer :: icv

  ! indices
  integer :: i_w, i_e
  integer :: j_s, j_n
  integer :: k_b, k_t

  integer :: np, type, i,j,k, imax, jmax, kmax, pmax, init_pc
  real(c_real) :: ic_vol, type_vol, acc_vol
  real(c_real) :: ic_dlo(3), ic_dhi(3)
  real(c_real) :: max_dp
  real(c_real) :: rp(3)
  real(c_real) :: x(3)

  real(c_real), allocatable :: dp(:), ro_s(:)

  init_pc = pc

  !  Set the initial conditions.
  do icv = 1, dim_ic
     if (ic_defined(icv) .and. abs(ic_ep_g(icv)-1.0d0)>epsilon(0.0d0)) then

        call calc_cell_ic(dx, dy, dz, &
             ic_x_w(icv), ic_y_s(icv), ic_z_b(icv), &
             ic_x_e(icv), ic_y_n(icv), ic_z_t(icv), &
             i_w, i_e, j_s, j_n, k_b, k_t)

        ! Start/end of IC domain bounds
        ic_dlo(1) = (max(lo(1), i_w)    ) * dx
        ic_dlo(2) = (max(lo(2), j_s)    ) * dy
        ic_dlo(3) = (max(lo(3), k_b)    ) * dz
        ic_dhi(1) = (min(hi(1), i_e) + 1) * dx
        ic_dhi(2) = (min(hi(2), j_n) + 1) * dy
        ic_dhi(3) = (min(hi(3), k_t) + 1) * dz

        ! physical volume of local piece of IC region
        ic_vol = (ic_dhi(1) - ic_dlo(1)) * &
                 (ic_dhi(2) - ic_dlo(2)) * &
                 (ic_dhi(3) - ic_dlo(3))

        do type=1, particle_types

           acc_vol = 0.0d0
           type_vol = ic_ep_s(icv,type) * ic_vol
           if(type_vol > 0.0d0) then

              if(is_defined(ic_dp_max(icv,type))) then
                 max_dp = ic_dp_max(icv,type)
              else
                 max_dp = ic_dp_mean(icv,type)
              endif

              ic_dlo = ic_dlo + 0.5d0*max_dp

              imax = floor((ic_dhi(1)-ic_dlo(1))/max_dp)
              jmax = floor((ic_dhi(2)-ic_dlo(2))/max_dp)
              kmax = floor((ic_dhi(3)-ic_dlo(3))/max_dp)
              pmax = imax*jmax*kmax

              rp(1) = 0.5d0 * (ic_dhi(1)-ic_dlo(1))/imax
              rp(2) = 0.5d0 * (ic_dhi(2)-ic_dlo(2))/jmax
              rp(3) = 0.5d0 * (ic_dhi(3)-ic_dlo(3))/kmax

              call grow_pdata(pmax)
              allocate(dp(pmax))

              ! Setup particle diameters
              if(ic_dp_dist(icv,type) == 'NORMAL') then
                 call nor_rno(dp, ic_dp_mean(icv,type), ic_dp_std(icv,type), &
                      ic_dp_min(icv,type), ic_dp_max(icv,type))

              else if(ic_dp_dist(icv,type) == 'UNIFORM') then
                 call uni_rno(dp, ic_dp_min(icv,type), ic_dp_max(icv,type))
              else
                 dp = ic_dp_mean(icv,type)
              endif

              np = 0
              klp: do k = 0, kmax-1

                 x(3) = (2.0d0*sqrt(6.0d0)/3.0d0 * k)

                 do j = 0, jmax-1

                    x(2) = (sqrt(3.0d0)*(j + (1.0d0/3.0d0)*mod(k,2)))

                    do i = 0, imax-1

                       x(1) = (2.0d0*i + (mod(j+k,2)))

                       np = np + 1
                       pc = pc + 1

                       rdata(pc,1:3) = ic_dlo + rp*x
                       rdata(pc,4) = 0.5d0*dp(np)

                       acc_vol = acc_vol + (pi/6.0d0)*dp(np)**3
                       if(acc_vol >= type_vol) exit klp

                    enddo
                 enddo
              enddo klp
              if(allocated(dp)) deallocate(dp)

              ! At this point, np particles were created. Now given
              ! then some properties

              allocate(ro_s(np))
              if(ic_ro_s_dist(icv,type) == 'NORMAL') then
                 call nor_rno(ro_s, ic_ro_s_mean(icv,type), ic_ro_s_std(icv,type), &
                      ic_ro_s_min(icv,type), ic_ro_s_max(icv,type))

              else if(ic_ro_s_dist(icv,type) == 'UNIFORM') then
                 call uni_rno(ro_s, ic_ro_s_min(icv,type), ic_ro_s_max(icv,type))
              else
                 ro_s = ic_ro_s_mean(icv,type)
              endif


              do i=1,np

                 !< Density................ 5
                 rdata(pc-np+i,5) = ro_s(i)

                 !< Linear velocity........ 6,7,8
                 rdata(pc-np+i,6) = ic_u_s(icv,type)
                 rdata(pc-np+i,7) = ic_v_s(icv,type)
                 rdata(pc-np+i,8) = ic_w_s(icv,type)

                 !< Type................... 1
                 idata(pc-np+i,1) = type

              enddo
              if(allocated(ro_s)) deallocate(ro_s)

           endif
        enddo

     endif
  enddo

  return
end subroutine particle_generator


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: GENERATE_PARTICLE_CONFIG                                !
!                                                                      !
!  Purpose: Generate particle configuration based on maximum particle  !
!           radius and filling from top to bottom within specified     !
!           bounds                                                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine mfix_particle_generator_prop(nrp, particles) &
     bind(C, name="mfix_particle_generator_prop")

  use particle_mod
  use constant, only: pi

  integer(c_int),   intent(in   ) :: nrp
  type(particle_t), intent(inout) :: particles(nrp)

  integer :: p

  real(c_real)   :: rad, rho
  real(c_real)   :: vol, mass, omoi

  do p = 1, nrp

     particles(p) % pos(1:3) = rdata(p,1:3)
     particles(p) % vel(1:3) = rdata(p,6:8)

     rad  = rdata(p,4)
     rho  = rdata(p,5)

     vol  = (4.0d0/3.0d0)*pi*rad**3
     mass = vol * rho
     omoi = 2.5d0/(mass * rad**2)

     particles(p) % radius   = rad
     particles(p) % density  = rho

     particles(p) % volume   = vol
     particles(p) % mass     = mass
     particles(p) % omoi     = omoi

     particles(p) % omega    = 0.0d0
     particles(p) % drag     = 0.0d0

     particles(p) % phase    = idata(p,1)
     particles(p) % state    = 1

  end do

  if(allocated(rdata)) deallocate(rdata)
  if(allocated(idata)) deallocate(idata)

end subroutine mfix_particle_generator_prop

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                     !
!                                                                     !
!                                                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine nor_rno(dp, mean, sigma, dp_min, dp_max)

  implicit none

  real(c_real), intent(inout) :: dp(:)
  real(c_real), intent(in   ) :: mean, sigma, dp_min, dp_max

  ! Local variables
  !-----------------------------------------------
  real(c_real) :: lmean, lvariance, lsigma
  real(c_real) :: x(2), w, dp1, dp2
  integer i, nsize
  logical :: debug = .false.
  !-----------------------------------------------

  nsize = size(dp(:))
  ! call init_random_seed

  i=1
  do while(i<= ceiling(real(nsize/2.0)))
     w=1.0
     do while(w>=1.0)
        call random_number(x)
        x = 2.0 * x - 1.0
        w = x(1)**2 + x(2)**2
     end do

     w = sqrt( (-2.0 * log( w ) ) / w )

     dp1 = x(1) * w * sigma + mean
     dp2 = x(2) * w * sigma + mean

     if(dp1 >= dp_min .and. dp2 >= dp_min .and. &
        dp1 <= dp_max .and. dp2 <= dp_max) then
        if(2*i -1 >=     1) dp(2*i-1) = dp1
        if(2*i    <= nsize) dp(2*i  ) = dp2
        i= i+1
     endif
  end do


  if(debug) then
     lmean = sum(dp(:))/nsize

     lvariance = 0.0
     do i = 1, nsize
        lvariance = lvariance + (dp(i)-lmean)**2
     end do

     lvariance = lvariance/nsize
     lsigma = sqrt(lvariance)

     write(*,*) '   '
     write(*,1000) ! Divider
     write(*,1010) ! Header
     write(*,1000) ! Divider
     write(*,1020) '  Mean     ', mean, lmean
     write(*,1000) ! Divider
     write(*,1020) '  Sigma    ', sigma, lsigma
     write(*,1000) ! Divider

  endif

1000 format(4x,'+',11('-'),2('+',17('-')),'+')
1010 format(4x,'|',11(' '),'|    Specified    |    Computed     |')
1020 format(4x,'|',A11,2('|',1x,es15.6,1x),'|')


  return
end subroutine nor_rno


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                     !
!                                                                     !
!                                                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine uni_rno(dp, dp_min, dp_max)

  implicit none

  real(c_real), intent(inout) :: dp(:)
  real(c_real), intent(in   ) :: dp_min, dp_max

  integer :: nsize, lc
  real(c_real) :: lscale

  ! call init_random_seed
  call random_number(dp)

  lscale = dp_max - dp_min

  nsize = size(dp(:))
  do lc = 1, nsize
     dp(lc) = dp_min + lscale*dp(lc)

  enddo

  return

end subroutine uni_rno


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                     !
!                                                                     !
!                                                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine init_random_seed

  implicit none

  !-----------------------------------------------
  ! local variables
  !-----------------------------------------------
  integer              :: isize,idate(8)
  integer,allocatable  :: iseed(:)
        !-----------------------------------------------

  call date_and_time(values=idate)
  call random_seed(size=isize)
  allocate( iseed(isize) )
  call random_seed(get=iseed)
  iseed = iseed * (idate(8)-500) ! idate(8) contains millisecond
  call random_seed(put=iseed)

  if(allocated(iseed)) deallocate( iseed )

end subroutine init_random_seed



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                     !
!                                                                     !
!                                                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine grow_pdata(gsize)
  implicit none

  integer, intent(in) :: gsize

  integer :: csize, nsize

  real(c_real),   allocatable :: rtmp(:,:)
  integer(c_int), allocatable :: itmp(:,:)

  ! Increase real data
  if(.not.(allocated(rdata))) then
     allocate(rdata(gsize,nr))
  else

     csize = size(rdata,1)
     nsize = csize + gsize
     allocate(rtmp(nsize,nr))

     rtmp(1:csize,:) = rdata(1:csize,:)
     call move_alloc(rtmp,rdata)
  endif

  ! Increase integer data
  if(.not.(allocated(idata))) then
     allocate(idata(gsize,ni))
  else

     csize = size(idata,1)
     nsize = csize + gsize
     allocate(itmp(nsize,ni))

     itmp(1:csize,:) = idata(1:csize,:)
     call move_alloc(itmp,idata)
  endif

end subroutine grow_pdata


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                     !
!                                                                     !
!                                                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine particle_write(nrp, particles) &
     bind(C, name="mfix_particle_write")

  use particle_mod

  implicit none

  integer(c_int),   intent(in) :: nrp
  type(particle_t), intent(in) :: particles(nrp)

  integer :: lc1

  open(unit=100, file='test.vtp', status='unknown')

! Write the necessary header information for a PolyData file type
  write(100,"(A)")'<?xml version="1.0"?>'
  write(100,"(2A)") '<VTKFile type="PolyData"',&
       ' version="0.1" byte_order="LittleEndian">'
  write(100,"(3x,A)") '<PolyData>'

! Write Piece tag and identify the number of particles in the system.
  write(100,"(6x,a,i10.10,a,a)") &
       '<Piece NumberOfPoints="',nrp, '" NumberOfVerts="0" ', &
       'NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">'

  write(100,"(9x,a)")'<PointData>'


  write(100,"(12x,a)") '<DataArray type="Float32" Name="radius" &
       &NumberOfComponents="1" format="ascii">'
  do lc1 = 1, nrp
     write (100,"(15x,es13.6)") real(particles(lc1) % radius)
  end do
  write(100,"(12x,a)") '</DataArray>'

  write(100,"(12x,a)") '<DataArray type="Float32" Name="density" &
       &NumberOfComponents="1" format="ascii">'
  do lc1 = 1, nrp
     write (100,"(15x,es13.6)") real(particles(lc1) % density)
  end do
  write(100,"(12x,a)") '</DataArray>'


  write(100,"( 9x,a)") '</PointData>'

  write(100,"(9x,a)") '<Points>'
  write(100,"(12x,a,a)") '<DataArray type="Float32" ',&
       'Name="Position" NumberOfComponents="3" format="ascii">'
  do lc1 = 1,nrp
     write (100,"(15x,3(es13.6,3x))") real(particles(lc1) % pos)
  enddo
  write(100,"(12x,a,/9x,a)")'</DataArray>','</Points>'

! Write tags for data not included (vtp format style)
  write(100,"(9x,a,/9x,a,/9x,a,/9x,a)")'<Verts></Verts>',&
       '<Lines></Lines>','<Strips></Strips>','<Polys></Polys>'
  write(100,"(6x,a,/3x,a,/a)")&
       '</Piece>','</PolyData>','</VTKFile>'

  close(100)

return
end subroutine particle_write

end module par_gen_module
