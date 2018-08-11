module par_gen_module

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int

  implicit none

  real(rt), allocatable :: rdata(:,:)
  integer,      allocatable :: idata(:,:)

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
  use ic, only: ic_ep_s, ic_ep_g, ic_pack

  use ic, only: ic_x_e, ic_y_n, ic_z_t
  use ic, only: ic_x_w, ic_y_s, ic_z_b

  use ic, only: ic_u_s, ic_v_s, ic_w_s

  use ic, only: ic_dp_dist, ic_ro_s_dist
  use ic, only: ic_dp_mean, ic_ro_s_mean
  use ic, only: ic_dp_std,  ic_ro_s_std
  use ic, only: ic_dp_min,  ic_ro_s_min
  use ic, only: ic_dp_max,  ic_ro_s_max

  use param, only: is_defined

  use discretelement, only: particle_types
  use constant, only: pi

  implicit none

  integer(c_int), intent(inout) :: pc
  integer(c_int), intent(in   ) :: lo(3),hi(3)
  real(rt),   intent(in   ) :: dx, dy, dz


  real(rt), parameter :: sqrt3 = sqrt(3.0)
  real(rt), parameter :: sqrt6o3x2 = 2.0*sqrt(6.0)/3.0

  ! local index for initial condition
  integer :: icv

  ! indices
  integer  :: p, i_w, i_e, j_s, j_n, k_b, k_t
  integer  :: np, type, init_pc

  real(rt) :: ic_vol, type_vol, acc_vol, pvol
  real(rt) :: ic_dlo(3), ic_dhi(3)
  real(rt) :: mean_dp, max_dp, max_rp
  real(rt) :: pos(3)

  real(rt), allocatable :: dp(:), ro_s(:)

  init_pc = pc

  ! Get the IC index
  do icv = 1, dim_ic
     if (ic_defined(icv) .and. abs(ic_ep_g(icv)-1.0d0)>epsilon(0.0d0)) exit
  enddo

  ! Get the solids type index
  do type=1, particle_types
     if(ic_ep_s(icv,type) > epsilon(0.d0)) exit
  enddo

  select case(trim(ic_pack(icv)))
  case('HCP'   ); call hex_close_pack(icv, type, lo, hi, np, pc, dx, dy, dz)
  case('RANDOM'); call random_fill(icv, type, lo, hi, np, pc, dx, dy, dz)
  case DEFAULT
     write(*,*) "Unknown particle generator fill type"
     stop 1000
  end select

  ! No more work.
  if(np == 0) return

  allocate(dp(np))
  allocate(ro_s(np))

  ! Setup particle diameters
  if(ic_dp_dist(icv,type) == 'NORMAL') then
     call nor_rno(dp, ic_dp_mean(icv,type), ic_dp_std(icv,type), &
          ic_dp_min(icv,type), ic_dp_max(icv,type))

  else if(ic_dp_dist(icv,type) == 'UNIFORM') then
     call uni_rno(dp, ic_dp_min(icv,type), ic_dp_max(icv,type))
  else
     dp = ic_dp_mean(icv,type)
  endif

  if(ic_ro_s_dist(icv,type) == 'NORMAL') then
     call nor_rno(ro_s, ic_ro_s_mean(icv,type), ic_ro_s_std(icv,type), &
          ic_ro_s_min(icv,type), ic_ro_s_max(icv,type))

  else if(ic_ro_s_dist(icv,type) == 'UNIFORM') then
     call uni_rno(ro_s, ic_ro_s_min(icv,type), ic_ro_s_max(icv,type))
  else
     ro_s = ic_ro_s_mean(icv,type)
  endif

  pc = init_pc
  acc_vol = 0.0d0
  nplp: do p=1,np

     pvol = (pi/6.0d0)*dp(p)**3

     pc = pc + 1

     !< Radius................. 4
     rdata(pc,4) = 0.5d0*dp(p)

     !< Density................ 5
     rdata(pc,5) = ro_s(p)

     !< Linear velocity........ 6,7,8
     rdata(pc,6) = ic_u_s(icv,type)
     rdata(pc,7) = ic_v_s(icv,type)
     rdata(pc,8) = ic_w_s(icv,type)

     !< Type................... 1
     idata(pc,1) = type

  enddo nplp

  if(allocated(ro_s)) deallocate(ro_s)
  if(allocated(dp  )) deallocate(dp)


  return
end subroutine particle_generator

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: hex_close_pack                                          !
!                                                                      !
!  Purpose: Generate initial solids packing based on hexagonal close   !
!           packing of mono-sized spheres.                             !
!                                                                      !
!  TODO: * generalize fill direction to follow gravity.                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine hex_close_pack(icv, type, lo, hi, np, pc, dx, dy, dz)

  use ic, only: dim_ic, ic_defined
  use ic, only: ic_ep_s, ic_ep_g

  use ic, only: ic_x_e, ic_y_n, ic_z_t
  use ic, only: ic_x_w, ic_y_s, ic_z_b

  use ic, only: ic_dp_dist, ic_ro_s_dist
  use ic, only: ic_dp_mean, ic_ro_s_mean
  use ic, only: ic_dp_std,  ic_ro_s_std
  use ic, only: ic_dp_min,  ic_ro_s_min
  use ic, only: ic_dp_max,  ic_ro_s_max

  use param, only: is_defined

  use calc_cell_module, only: calc_cell_ic
  use discretelement, only: particle_types
  use constant, only: pi

  implicit none

  integer(c_int), intent(in   ) :: icv, type, lo(3), hi(3)
  integer(c_int), intent(inout) :: np, pc
  real(rt),       intent(in   ) :: dx, dy, dz

  real(rt), parameter :: sqrt3 = sqrt(3.0)
  real(rt), parameter :: sqrt6o3x2 = 2.0*sqrt(6.0)/3.0

  ! indices
  integer :: i_w, i_e
  integer :: j_s, j_n
  integer :: k_b, k_t

  integer  :: i,j,k, init_pc
  real(rt) :: ic_vol, type_vol, acc_vol, pvol
  real(rt) :: ic_dlo(3), ic_dhi(3)
  real(rt) :: mean_dp, max_dp, max_rp
  real(rt) :: pos(3)

  integer :: seed, max_seed(3), seed_lo(3), seed_hi(3)


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

! physical volume of IC region
  ic_vol = (ic_x_e(icv) - ic_x_w(icv)) * &
           (ic_y_n(icv) - ic_y_s(icv)) * &
           (ic_z_t(icv) - ic_z_b(icv))

  ! Spacing is based on maximum particle size
  if(is_defined(ic_dp_max(icv,type))) then
     max_dp = ic_dp_max(icv,type)
  else
     max_dp = ic_dp_mean(icv,type)
  endif
  max_rp = 0.5d0 * max_dp

  ! Particle count is based on mean particle size
  seed = ic_vol * ic_ep_s(icv,type) / &
       ((pi/6.0d0)*ic_dp_mean(icv,type)**3)

  ! Total to seed over the whole IC region
  max_seed(1) = int((ic_x_e(icv) - ic_x_w(icv) - max_dp)/max_dp)
  max_seed(3) = int((ic_z_t(icv) - ic_z_b(icv) - max_dp)/(sqrt3*max_rp))
  max_seed(2) = int(seed / (max_seed(1)*max_seed(3)))

  ! local grid seed loop hi/lo
  seed_lo(1) = nint((ic_dlo(1) - i_w*dx) / max_dp)
  seed_lo(3) = nint((ic_dlo(3) - k_b*dz) / (sqrt3 * max_rp))
  seed_lo(2) = nint((ic_dlo(2) - j_s*dy) / ((sqrt6o3x2) * max_rp))

  seed_hi(1) = nint((ic_dhi(1) - i_w*dx) /  max_dp - seed_lo(1)*max_dp)
  seed_hi(3) = nint((ic_dhi(3) - k_b*dz) / (sqrt3 * max_rp) - seed_lo(1)*max_dp)
  seed_hi(2) = nint((ic_dhi(2) - j_s*dy) / ((sqrt6o3x2) * max_rp) - seed_lo(1)*max_dp)

  seed_hi(1) = min(max_seed(1), seed_hi(1)-1)
  seed_hi(3) = min(max_seed(3), seed_hi(3)-1)
  seed_hi(2) = min(max_seed(2), seed_hi(2)-1)

  pos = -1.0d20
  np = 0

  do j=seed_lo(2), seed_hi(2)

     pos(2) = j_s*dy + max_rp*(1.0d0 + j*sqrt6o3x2)

     do k=seed_lo(3), seed_hi(3)

        pos(3) = k_b*dz + max_rp*(1.0d0 + sqrt3*(k+(mod(j,2)/3.0d0)))

        do i=seed_lo(1), seed_hi(1)

           pos(1) = i_w*dx + max_rp* (1.0d0 + 2.0d0*i + mod(j+k,2))

           np = np + 1 ! local to type
           pc = pc + 1 ! local to routine

           call grow_pdata(pc)

           rdata(pc,1:3) = pos

        enddo
     enddo
  enddo

  return
end subroutine hex_close_pack

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: random fill                                             !
!                                                                      !
!  Purpose: Generate initial solids packing based on randomly placing  !
!           particles in the ic region.                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine random_fill(icv, type, lo, hi, np, pc, dx, dy, dz)

  use ic, only: dim_ic, ic_defined
  use ic, only: ic_ep_s, ic_ep_g

  use ic, only: ic_x_e, ic_y_n, ic_z_t
  use ic, only: ic_x_w, ic_y_s, ic_z_b

  use ic, only: ic_dp_dist, ic_ro_s_dist
  use ic, only: ic_dp_mean, ic_ro_s_mean
  use ic, only: ic_dp_std,  ic_ro_s_std
  use ic, only: ic_dp_min,  ic_ro_s_min
  use ic, only: ic_dp_max,  ic_ro_s_max

  use param, only: is_defined

  use calc_cell_module, only: calc_cell_ic
  use discretelement, only: particle_types
  use constant, only: pi

  implicit none

  integer(c_int), intent(in   ) :: icv, type, lo(3), hi(3)
  integer(c_int), intent(inout) :: np, pc
  real(rt),       intent(in   ) :: dx, dy, dz

  real(rt), parameter :: sqrt3 = sqrt(3.0)
  real(rt), parameter :: sqrt6o3x2 = 2.0*sqrt(6.0)/3.0

  ! indices
  integer :: i_w, i_e
  integer :: j_s, j_n
  integer :: k_b, k_t

  integer  :: i,j,k,l,ll,ii,jj,kk
  real(rt) :: ic_vol, type_vol, acc_vol, pvol
  real(rt) :: ic_dlo(3), ic_dhi(3), ic_len(3)
  real(rt) :: mean_dp, max_dp, max_rp
  real(rt) :: pos(3), Oodx(3), rand3(3), mindist, dist

  integer, allocatable :: pinc(:,:,:), pbin(:,:,:,:), tbin(:,:,:,:)

  integer :: seed, overlaps, fails, nb, ob

  integer, parameter :: maxfails = 1000

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

! physical volume of IC region intersecting this grid
  ic_vol = (ic_dhi(1) - ic_dlo(1)) * &
           (ic_dhi(2) - ic_dlo(2)) * &
           (ic_dhi(3) - ic_dlo(3))

  ! Spacing is based on maximum particle size
  if(is_defined(ic_dp_max(icv,type))) then
     max_dp = ic_dp_max(icv,type)
  else
     max_dp = ic_dp_mean(icv,type)
  endif
  max_rp = 0.5d0 * max_dp

  ! Particle count is based on mean particle size
  seed = ic_vol * ic_ep_s(icv,type) / &
       ((pi/6.0d0)*ic_dp_mean(icv,type)**3)

  ic_len = ic_dhi - ic_dlo - max_dp
  ic_dlo = ic_dlo + max_rp

  pos = -1.0d20
  np = 0

  Oodx(1) = 1.0_rt/dx
  Oodx(2) = 1.0_rt/dy
  Oodx(3) = 1.0_rt/dz

  mindist = (1.05d0*max_dp)**2

  allocate(pinc(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
  allocate(pbin(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),8))

  np = 0
  fails = 0
  pinc = 0

  do while (np < seed .and. fails < maxfails)

     do

        call random_number(rand3)
        pos = ic_dlo + ic_len*rand3(:)

        ! Grid containing the new particle
        i = floor(pos(1)*Oodx(1) + 0.5_rt)
        j = floor(pos(2)*Oodx(2) + 0.5_rt)
        k = floor(pos(3)*Oodx(3) + 0.5_rt)

        ! Local grid search for collisions.
        overlaps=0
        do kk=max(lo(3),k-1), min(k+1,hi(3))
           do jj=max(lo(2),j-1), min(j+1,hi(2))
              do ii=max(lo(1),i-1), min(i+1,hi(1))

                 do l=1, pinc(ii,jj,kk)

                    ll = pbin(ii,jj,kk,l)

                    dist=(rdata(ll,1) - pos(1))**2 + &
                         (rdata(ll,2) - pos(2))**2 + &
                         (rdata(ll,3) - pos(3))**2

                    if(dist < mindist) overlaps = overlaps+1

                 enddo
              enddo
           enddo
        enddo

        if(overlaps == 0) then

           np = np + 1 ! local to call
           pc = pc + 1 ! local to grid

           call grow_pdata(pc)

           rdata(pc,1:3) = pos

           pinc(i,j,k) = pinc(i,j,k) + 1
           pbin(i,j,k,pinc(i,j,k)) = np

           ob = ubound(pbin,4)

           if(pinc(i,j,k) + 1 >= ob) then
              nb = ob + 2
              allocate(tbin(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),nb))
              tbin(:,:,:,1:ob) = pbin(:,:,:,1:ob)
              call move_alloc(tbin, pbin)
           endif

           exit
        else
           fails = fails + 1
        endif
     enddo

     if((mod(np, seed/10) == 0 .and. np < seed*0.95) .or. np==seed) &
          write(*,"(2x,'Seeded: ',I9,3x,'(',f5.0,'%)')") np,100*dble(np)/seed

  enddo


  return
end subroutine random_fill

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

  real(rt)   :: rad, rho
  real(rt)   :: vol, mass, omoi

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

  real(rt), intent(inout) :: dp(:)
  real(rt), intent(in   ) :: mean, sigma, dp_min, dp_max

  ! Local variables
  !-----------------------------------------------
  real(rt) :: lmean, lvariance, lsigma
  real(rt) :: x(2), w, dp1, dp2
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

  real(rt), intent(inout) :: dp(:)
  real(rt), intent(in   ) :: dp_min, dp_max

  integer :: nsize, lc
  real(rt) :: lscale

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

  real(rt),   allocatable :: rtmp(:,:)
  integer(c_int), allocatable :: itmp(:,:)

  ! Increase real data
  if(.not.(allocated(rdata))) then
     allocate(rdata(max(gsize,1024),nr))
  else

     csize = size(rdata,1)
     if(gsize >= csize) then
        nsize = max(2*csize, gsize)
        allocate(rtmp(nsize,nr))
        rtmp(1:csize,:) = rdata(1:csize,:)
        call move_alloc(rtmp,rdata)
     endif
  endif

  ! Increase integer data
  if(.not.(allocated(idata))) then
     allocate(idata(max(gsize,1024),ni))
  else

     csize = size(idata,1)
     if(gsize >= csize) then
        nsize = max(2*csize, gsize)
        allocate(itmp(nsize,ni))
        itmp(1:csize,:) = idata(1:csize,:)
        call move_alloc(itmp,idata)
     endif
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



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                     !
!                                                                     !
!                                                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine rm_wall_collisions ( particles, nrp, flag, fglo, fghi, &
     bcent, blo, bhi, apx, axlo, axhi, apy, aylo, ayhi, apz, azlo, azhi, dx) &
     bind(C, name="rm_wall_collisions")

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding    , only: c_int


  use param        , only: small_number, zero, one
  use particle_mod,  only: particle_t

  use amrex_ebcellflag_module, only : is_regular_cell, is_covered_cell, is_single_valued_cell
  use amrex_ebcellflag_module, only : get_neighbor_cells

  implicit none

  integer, intent(in) ::  nrp

  type(particle_t), intent(inout), target :: particles(nrp)

  integer, dimension(3), intent(in) :: axlo, aylo, azlo, fglo, blo
  integer, dimension(3), intent(in) :: axhi, ayhi, azhi, fghi, bhi

  integer,      intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))
  real(rt), intent(in) :: bcent(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),3)
  real(rt), intent(in) :: apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
  real(rt), intent(in) :: apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
  real(rt), intent(in) :: apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
  real(rt), intent(in) :: dx(3)

  type(particle_t), pointer :: p

  ! facet barycenter (bcent) in global coordinates
  real(rt) :: eb_cent(3)

  integer :: ll, ii, jj, kk, i, j, k

  real(rt), parameter :: fudge = 1.0d0 - 1.0d-8
  real(rt) :: inv_dx(3)

  real(rt) :: adx, ady, adz, apnorminv

  integer :: nbr(-1:1,-1:1,-1:1)

  real(rt) :: dist, normal(3)

  ! inverse cell size: used to convert positions to cell indices
  !   -> dx is a vector, fortran is amazing!
    inv_dx = 1.0d0 / dx

    ! itterate over particles
    do ll = 1, nrp
       ! get current particle
       p => particles(ll)

       ! cell indices for position corresponding to (lx, ly, lz)
       i = floor(p%pos(1)*inv_dx(1))
       j = floor(p%pos(2)*inv_dx(2))
       k = floor(p%pos(3)*inv_dx(3))

       if(is_covered_cell(flag(i,j,k))) then
          p%id = -1
       else

          ! ignore disconnected cells
          ! -> one reason for this is the get accurate wall normals...
          call get_neighbor_cells(flag(i,j,k),nbr)

          do kk = k-1, k+1
             do jj = j-1, j+1
                do ii = i-1, i+1
                   ! only consider cells that contain EB's
                   if ( (nbr(ii-i, jj-j, kk-k) == 1) .and. &
                        (.not. is_regular_cell(flag(ii, jj, kk))) ) then

                      adx = apx(ii,jj,kk) - apx(ii+1,jj  , kk  )
                      ady = apy(ii,jj,kk) - apy(ii,  jj+1, kk  )
                      adz = apz(ii,jj,kk) - apz(ii,  jj  , kk+1)

                      apnorminv = 1.d0 / sqrt(adx**2 + ady**2 + adz**2)

                      normal(1) = adx * apnorminv   ! pointing to the wall
                      normal(2) = ady * apnorminv
                      normal(3) = adz * apnorminv

                      ! convert bcent to global coordinate system centered at plo
                      eb_cent(:) = ( bcent(ii, jj, kk, :) + &
                           (/dble(ii), dble(jj), dble(kk)/) + &
                           (/0.5d0, 0.5d0, 0.5d0/) )*dx(:)

                      ! Distance to closest point on EB
                      dist = dot_product( p%pos(:) - eb_cent(:), -normal(:))

                      if (dist .lt. p%radius) p%id = -1
                   end if
                end do
             end do
          end do

       endif

    end do ! loop over particles
  end subroutine rm_wall_collisions
end module par_gen_module
