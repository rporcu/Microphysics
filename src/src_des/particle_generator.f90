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
      use ic, only: ic_ep_s, ic_ep_g, ic_pack_type

      use ic, only: ic_u_s, ic_v_s, ic_w_s

      use ic, only: ic_dp_dist, ic_ro_s_dist
      use ic, only: ic_dp_mean, ic_ro_s_mean
      use ic, only: ic_dp_std,  ic_ro_s_std
      use ic, only: ic_dp_min,  ic_ro_s_min
      use ic, only: ic_dp_max,  ic_ro_s_max

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
      integer  :: p
      integer  :: np, type, init_pc

      real(rt) :: pvol

      real(rt), allocatable :: dp(:), ro_s(:)

      init_pc = pc

      ! Get the IC index
      do icv = 1, dim_ic
         if (ic_defined(icv) .and. abs(ic_ep_g(icv)-1.0d0)>epsilon(0.0d0)) exit
      enddo

      if(icv > dim_ic) return

      ! Get the solids type index
      do type=1, particle_types
         if(ic_ep_s(icv,type) > epsilon(0.d0)) exit
      enddo

      select case(trim(ic_pack_type(icv)))
          case('HCP'   ); call hex_close_pack(icv, type, lo, hi, np, pc, dx, dy, dz)
          case('RANDOM'); call random_fill(icv, type, lo, hi, np, pc, dx, dy, dz, .false.)
          case('PSEUDO_RANDOM'); call random_fill(icv, type, lo, hi, np, pc, dx, dy, dz, .true. )
          case('ONEPER'); call one_per_fill(icv, type, lo, hi, np, pc, dx, dy, dz)
          case('EIGHTPER'); call eight_per_fill(icv, type, lo, hi, np, pc, dx, dy, dz)
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
      use ic, only: ic_ep_s

      use ic, only: ic_x_e, ic_y_n, ic_z_t
      use ic, only: ic_x_w, ic_y_s, ic_z_b

      use ic, only: ic_dp_mean
      use ic, only: ic_dp_max

      use param, only: is_defined

      use calc_cell_module, only: calc_cell_ic
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

      integer  :: i,j,k
      real(rt) :: ic_vol
      real(rt) :: ic_dlo(3), ic_dhi(3)
      real(rt) :: max_dp, max_rp
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
   !  Subroutine: one_per_fill                                            !
   !                                                                      !
   !  Purpose: Generate initial solids packing based on putting one       !
   !           per cell                                                   !
   !                                                                      !
   !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine one_per_fill(icv, type, lo, hi, np, pc, dx, dy, dz)

      use ic, only: dim_ic, ic_defined
      use ic, only: ic_ep_s

      use ic, only: ic_x_e, ic_y_n, ic_z_t
      use ic, only: ic_x_w, ic_y_s, ic_z_b

      use ic, only: ic_dp_mean
      use ic, only: ic_dp_max

      use param, only: is_defined

      use calc_cell_module, only: calc_cell_ic
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

      integer  :: i,j,k
      real(rt) :: ic_vol
      real(rt) :: ic_dlo(3), ic_dhi(3)
      real(rt) :: max_dp, max_rp
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

      ! This routine arbitrarily puts one particle at the center of each cell.
      do j = lo(2), hi(2)
         pos(2) = (j + 0.5d0)*dy
         if (pos(2).ge.ic_y_s(icv) .and. pos(2).le.ic_y_n(icv)) then

            do k = lo(3), hi(3)
                pos(3) = (k + 0.5d0)*dz
                if (pos(3).ge.ic_z_b(icv) .and. pos(3).le.ic_z_t(icv)) then

                   do i = lo(1), hi(1)
                       pos(1) = (i + 0.5d0)*dx
                       if (pos(1).ge.ic_x_w(icv) .and. pos(1).le.ic_x_e(icv)) then

                            np = np + 1 ! local to type
                            pc = pc + 1 ! local to routine

                            call grow_pdata(pc)

                            rdata(pc,1:3) = pos
                       end if
                   enddo

                end if
            enddo
         end if
      enddo

   end subroutine one_per_fill

   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
   !                                                                      !
   !  Subroutine: eight_per_fill                                            !
   !                                                                      !
   !  Purpose: Generate initial solids packing based on putting eight     !
   !           per cell                                                   !
   !                                                                      !
   !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine eight_per_fill(icv, type, lo, hi, np, pc, dx, dy, dz)

      use ic, only: dim_ic, ic_defined
      use ic, only: ic_ep_s

      use ic, only: ic_x_e, ic_y_n, ic_z_t
      use ic, only: ic_x_w, ic_y_s, ic_z_b

      use ic, only: ic_dp_mean
      use ic, only: ic_dp_max

      use param, only: is_defined

      use calc_cell_module, only: calc_cell_ic
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

      integer  :: i,j,k
      real(rt) :: ic_vol
      real(rt) :: ic_dlo(3), ic_dhi(3)
      real(rt) :: max_dp, max_rp
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

      ! This routine arbitrarily puts eight particles per cell, one at the center
      !      of each quadrant

      do j = 2*lo(2), 2*hi(2)+1
         pos(2) = (j + 0.5d0)*dy/2.d0
         if (pos(2).ge.ic_y_s(icv) .and. pos(2).le.ic_y_n(icv)) then

            do k = 2*lo(3), 2*hi(3)+1
                pos(3) = (k + 0.5d0)*dz/2.d0
                if (pos(3).ge.ic_z_b(icv) .and. pos(3).le.ic_z_t(icv)) then

                   do i = 2*lo(1), 2*hi(1)+1
                       pos(1) = (i + 0.5d0)*dx/2.d0
                       if (pos(1).ge.ic_x_w(icv) .and. pos(1).le.ic_x_e(icv)) then

                            np = np + 1 ! local to type
                            pc = pc + 1 ! local to routine

                            call grow_pdata(pc)

                            rdata(pc,1:3) = pos
                       end if
                   enddo

                end if
            enddo
         end if
      enddo

   end subroutine eight_per_fill

   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
   !                                                                      !
   !  Subroutine: random fill                                             !
   !                                                                      !
   !  Purpose: Generate initial solids packing based on randomly placing  !
   !           particles in the ic region.                                !
   !                                                                      !
   !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine random_fill(icv, type, lo, hi, np, pc, dx, dy, dz, fix_seed)

      use ic, only: dim_ic, ic_defined
      use ic, only: ic_ep_s

      use ic, only: ic_x_e, ic_y_n, ic_z_t
      use ic, only: ic_x_w, ic_y_s, ic_z_b

      use ic, only: ic_dp_mean
      use ic, only: ic_dp_max

      use param, only: is_defined

      use calc_cell_module, only: calc_cell_ic
      use constant, only: pi

      implicit none

      integer(c_int), intent(in   ) :: icv, type, lo(3), hi(3)
      integer(c_int), intent(inout) :: np, pc
      logical       , intent(in   ) :: fix_seed
      real(rt),       intent(in   ) :: dx, dy, dz

      real(rt), parameter :: sqrt3 = sqrt(3.0)
      real(rt), parameter :: sqrt6o3x2 = 2.0*sqrt(6.0)/3.0

      ! indices
      integer :: i_w, i_e
      integer :: j_s, j_n
      integer :: k_b, k_t

      integer  :: i,j,k,l,ll,ii,jj,kk
      real(rt) :: ic_vol
      real(rt) :: ic_dlo(3), ic_dhi(3), ic_len(3)
      real(rt) :: max_dp, max_rp
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

      mindist = (1.01d0*max_dp)**2

      allocate(pinc(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      allocate(pbin(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),8))

      np = 0
      fails = 0
      pinc = 0

      if (fix_seed) &
         call init_random_seed(fix_seed)

      do while (np < seed .and. fails < maxfails)

         do

            call random_number(rand3)
            pos = ic_dlo + ic_len*rand3(:)

            ! Grid containing the new particle
            i = floor(pos(1)*Oodx(1))
            j = floor(pos(2)*Oodx(2))
            k = floor(pos(3)*Oodx(3))

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
                  tbin(:,:,:,ob+1:nb) = 0
                  call move_alloc(tbin, pbin)
               endif
               fails = 0

               exit
            else
               fails = fails + 1
            endif
         enddo

         !     if((mod(np, seed/10) == 0 .and. np < seed*0.95) .or. np==seed) &
         !          write(*,"(2x,'Seeded: ',I9,3x,'(',f5.0,'%)')") np,100*dble(np)/seed

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
      ! call init_random_seed(.false.)

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

1000  format(4x,'+',11('-'),2('+',17('-')),'+')
1010  format(4x,'|',11(' '),'|    Specified    |    Computed     |')
1020  format(4x,'|',A11,2('|',1x,es15.6,1x),'|')


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

      ! call init_random_seed(.false.)
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
   subroutine init_random_seed(fix_seed)

      implicit none

      logical, intent(in)  :: fix_seed

      !-----------------------------------------------
      ! local variables
      !-----------------------------------------------
      integer              :: isize,idate(8)
      integer,allocatable  :: iseed(:)
      !-----------------------------------------------

      call random_seed(size=isize)
      allocate( iseed(isize) )
      call random_seed(get=iseed)

      ! Note -- "10" is arbitrary -- we just need something repeatable for
      !     regression testing
      if ( fix_seed ) then
         iseed(:) = 10
      else
         call date_and_time(values=idate)
         iseed = iseed * (idate(8)-500) ! idate(8) contains millisecond
      end if

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
   subroutine rm_wall_collisions ( particles, nrp,           &
    valid,     vlo,  vhi,     &
    phi,       phlo, phhi,    &
    dx,        n_refine     ) &
    bind(C, name="rm_wall_collisions")

      use particle_mod, only: particle_t

      use param       , only: small_number, zero, one

      use amrex_eb_levelset_module, only: amrex_eb_interp_levelset
      use amrex_eb_levelset_module, only: amrex_eb_normal_levelset

      implicit none

      ! ** input varaibles

      type(particle_t), intent(inout), target :: particles(nrp)
      integer,          intent(in   )         :: nrp
      integer,          intent(in   )         :: n_refine
      integer,          intent(in   )         :: vlo(3),  vhi(3)
      integer,          intent(in   )         :: phlo(3), phhi(3)

      integer,          intent(in   )         :: valid( vlo(1): vhi(1),  vlo(2): vhi(2),  vlo(3): vhi(3) )
      real(rt),         intent(in   )         :: phi(  phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3) )
      real(rt),         intent(in   )         :: dx (3)

      ! pos: current particle's position
      !  rp: current particle's radius
      real(rt) :: pos(3), rp

      !    plo: domain lo -- HACK this should get passed down
      real(rt) :: plo(3), inv_dx(3), ls_value

      integer :: p

      inv_dx = 1.0_rt / dx
      plo    = (/ 0.0_rt, 0.0_rt, 0.0_rt /)

      ! itterate over particles
      do p = 1, nrp

         rp  = particles(p)%radius
         pos = particles(p)%pos

         ! interpolates levelset from nodal phi to position pos
         call amrex_eb_interp_levelset(pos, plo, n_refine, &
          phi, phlo, phhi, dx, ls_value);

         if (ls_value < rp) particles(p)%id = -1

      end do ! loop over particles

   end subroutine rm_wall_collisions



   subroutine rm_wall_collisions_eb ( particles,    nrp, &
    &                                 valid,  vlo,  vhi, &
    &                                 phi,   phlo, phhi, &
    &                                 flags,  flo,  fhi, &
    &                                 plo, dx, n_refine ) bind(C)

      use particle_mod,             only: particle_t
      use param,                    only: small_number, zero, one
      use amrex_ebcellflag_module,  only: is_covered_cell
      use amrex_eb_levelset_module, only: amrex_eb_interp_levelset, &
       &                                  amrex_eb_normal_levelset

      ! Particles
      type(particle_t), intent(inout), target :: particles(nrp)
      integer,          intent(in   )         :: nrp

      ! Array bounds
      integer,          intent(in   )         :: vlo(3),  vhi(3)
      integer,          intent(in   )         :: phlo(3), phhi(3)
      integer,          intent(in   )         :: flo(3),  fhi(3)

      ! LS refinement
      integer,          intent(in   )         :: n_refine

      ! Grid spacing
      real(rt),         intent(in   )         :: dx (3)

      ! Coordinates of bottom-east-south corner of domain
      real(rt),         intent(in   )         :: plo(3)

      ! Arrays
      integer,          intent(in   )         :: &
       &    valid( vlo(1):vhi(1), vlo(2):vhi(2), vlo(3):vhi(3) ), &
       &    flags( flo(1):fhi(1), flo(2):fhi(2), flo(3):fhi(3) )

      real(rt),         intent(in   )         :: &
       &      phi(  phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3) )


      ! Local variables
      real(rt) :: odx(3), ls_value
      integer  :: p, ic, jc, kc

      odx = one / dx

      do p = 1, nrp

         associate( rp => particles(p)%radius, pos => particles(p)%pos )

            ! Indeces of the cells containing the particle center
            ic = floor( ( pos(1) - plo(1) ) * odx(1) )
            jc = floor( ( pos(2) - plo(2) ) * odx(2) )
            kc = floor( ( pos(3) - plo(3) ) * odx(3) )

            if (is_covered_cell(flags(ic,jc,kc))) then
               particles(p)%id = -1
            else
               ! interpolates level-set from nodal phi to position pos
               call amrex_eb_interp_levelset( pos, plo, n_refine, &
                &   phi, phlo, phhi, dx, ls_value)
               if (ls_value < rp) particles(p)%id = -1
            end if

         end associate

      end do

   end subroutine rm_wall_collisions_eb

end module par_gen_module
