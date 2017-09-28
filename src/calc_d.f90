module calc_d_mod

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   use param, only: zero, small_number

   ! Flag: Coupled DEM simulation
   use discretelement, only: des_continuum_coupled
   use discretelement, only: des_oneway_coupled

   ! Pressure scale factor
   use scales, only: p_scale
   use bc, only: minf_, nsw_, psw_, fsw_

   implicit none

   contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_d_*                                                !
!  Author: M. Syamlal                                 Date: 21-JUN-96  !
!                                                                      !
!  Purpose: calculate coefficients linking velocity correction to      !
!           pressure correction                                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine calc_d_e(lo, hi, slo, shi, ulo, uhi, alo, ahi, dlo, dhi, &
        d_e, A_m, ep_g, f_gds_u, dx, dy, dz, domlo, domhi, bc_ilo_type, bc_ihi_type)

      integer, intent(in   ) ::  lo(3), hi(3)
      integer, intent(in   ) :: slo(3),shi(3)
      integer, intent(in   ) :: ulo(3),uhi(3)
      integer, intent(in   ) :: alo(3),ahi(3)
      integer, intent(in   ) :: dlo(3),dhi(3)
      integer, intent(in   ) :: domlo(3),domhi(3)

      ! Pressure correction
      real(c_real), intent(  out) :: d_e&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      real(c_real), intent(in   ):: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3), -3:3)

      real(c_real), intent(in   ):: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: f_gds_u&
         (dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))

      integer(c_int), intent(in   ) :: bc_ilo_type&
           (domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_ihi_type&
           (domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2)

      real(c_real), intent(in   ) :: dx, dy, dz

      integer      :: i,j,k, bcv
      real(c_real) :: ayz
      real(c_real) :: Am0, epga
      logical      :: coupled

      coupled = (des_continuum_coupled .and. .not.des_oneway_coupled)

      ayz = dy*dz

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               Am0 = -A_m(i,j,k,0)

               if (abs(am0) > small_number) then

                  epga = ayz*0.5d0*(ep_g(i-1,j,k)+ep_g(i,j,k))
                  if (coupled) Am0 = Am0 + f_gds_u(i,j,k)
                  d_e(i,j,k) = p_scale*epga/am0

               else

                  d_e(i,j,k) = zero

               endif

            enddo
         enddo
      enddo

      ! At left boundary
      if (slo(1) .lt. domlo(1) .and. lo(1).eq.alo(1)) then
         i = alo(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               if(bc_ilo_type(j,k,1) == MINF_ .or. &
                  bc_ilo_type(j,k,1) == NSW_ .or. &
                  bc_ilo_type(j,k,1) == FSW_ .or. &
                  bc_ilo_type(j,k,1) == PSW_) then
                  d_e(i,j,k) =  0.0d0
               endif
            end do
         end do
      endif

      ! At right boundary
      if (shi(1) .gt. domhi(1) .and. hi(1).eq.ahi(1)) then
         i = ahi(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               if(bc_ihi_type(j,k,1) == MINF_ .or. &
                  bc_ihi_type(j,k,1) == NSW_  .or. &
                  bc_ihi_type(j,k,1) == FSW_  .or. &
                  bc_ihi_type(j,k,1) == PSW_) then
                  d_e(i,j,k) = 0.0d0
               endif

            end do
         end do
      endif


   end subroutine calc_d_e

   subroutine calc_d_n(lo, hi, slo, shi, vlo, vhi, alo, ahi, dlo, dhi, &
        d_n, A_m, ep_g, f_gds_v, dx, dy, dz, domlo, domhi, bc_jlo_type, bc_jhi_type)


      integer, intent(in   ) ::  lo(3), hi(3)
      integer, intent(in   ) :: slo(3),shi(3)
      integer, intent(in   ) :: vlo(3),vhi(3)
      integer, intent(in   ) :: alo(3),ahi(3)
      integer, intent(in   ) :: dlo(3),dhi(3)
      integer, intent(in   ) :: domlo(3),domhi(3)

      ! Pressure correction
      real(c_real), intent(  out) :: d_n&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

      real(c_real), intent(in   ):: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3), -3:3)

      real(c_real), intent(in   ):: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: f_gds_v&
         (dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))

      integer(c_int), intent(in   ) :: bc_jlo_type&
           (domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_jhi_type&
           (domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2)

      real(c_real), intent(in   ) :: dx, dy, dz

      integer      :: i,j,k, bcv
      real(c_real) :: axz
      real(c_real) :: Am0, epga
      logical      :: coupled

      coupled = (des_continuum_coupled .and. .not.des_oneway_coupled)

      axz = dx*dz

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               Am0 = -A_m(i,j,k,0)

               if(abs(Am0) > small_number) then

                  epga = axz*0.5d0*(ep_g(i,j-1,k)+ep_g(i,j,k))
                  if (coupled) Am0 = Am0 + f_gds_v(i,j,k)
                  d_n(i,j,k) = p_scale*epga/am0

               else

                  d_n(i,j,k) = zero

               endif

            enddo
         enddo
      enddo

      ! At bottom boundary
      if (slo(2) .lt. domlo(2) .and. lo(2).eq.alo(2)) then
         j = alo(2)
         do k = lo(3), hi(3)
            do i = lo(1), hi(1)
               bcv = bc_jlo_type(i,k,2)
               if(bc_jlo_type(i,k,1) == MINF_ .or. &
                  bc_jlo_type(i,k,1) == NSW_  .or. &
                  bc_jlo_type(i,k,1) == PSW_  .or. &
                  bc_jlo_type(i,k,1) == FSW_) then
                  d_n(i,j,k) =  zero
               endif
            end do
         end do
      endif

      ! At top boundary
      if (shi(2) .gt. domhi(2) .and. hi(2).eq.ahi(2)) then
         j = ahi(2)
         do k = lo(3), hi(3)
            do i = lo(1), hi(1)
               bcv = bc_jhi_type(i,k,2)
               if(bc_jhi_type(i,k,1) == MINF_ .or. &
                  bc_jhi_type(i,k,1) == NSW_  .or. &
                  bc_jhi_type(i,k,1) == PSW_  .or. &
                  bc_jhi_type(i,k,1) == FSW_ ) then
                  d_n(i,j,k) =  zero
               endif
            end do
         end do
      endif

   end subroutine calc_d_n

   subroutine calc_d_t(lo, hi, slo, shi, wlo, whi, alo, ahi, dlo, dhi, &
      d_t, A_m, ep_g, f_gds_w, dx, dy, dz, domlo, domhi, bc_klo_type, bc_khi_type)

      integer, intent(in   ) ::  lo(3), hi(3)
      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) :: wlo(3),whi(3)
      integer     , intent(in   ) :: alo(3),ahi(3)
      integer, intent(in   ) :: dlo(3),dhi(3)
      integer, intent(in   ) :: domlo(3),domhi(3)

      ! Pressure correction
      real(c_real), intent(  out) :: d_t&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ):: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3), -3:3)

      real(c_real), intent(in   ):: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: f_gds_w&
         (dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))

      integer(c_int), intent(in   ) :: bc_klo_type&
           (domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2)
      integer(c_int), intent(in   ) :: bc_khi_type&
           (domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2)

      real(c_real), intent(in   ) :: dx, dy, dz

      integer      :: i,j,k, bcv
      real(c_real) :: axy
      real(c_real) :: Am0, epga
      logical      :: coupled

      coupled = (des_continuum_coupled .and. .not.des_oneway_coupled)

      axy = dx*dy

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              Am0 = -A_m(I,J,K,0)

              if (abs(Am0) > small_number) THEN

                 epga = axy*0.5d0*(ep_g(i,j,k-1)+ep_g(i,j,k))
                 if (coupled) Am0 = Am0 + f_gds_w(i,j,k)
                 d_t(i,j,k) = p_scale*epga/am0

              else

                 d_t(i,j,k) = zero

              endif

           enddo
        enddo
     enddo

      ! At down boundary
     if (slo(3) .lt. domlo(3) .and. lo(3).eq.alo(3)) then
        k = alo(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              bcv = bc_klo_type(i,j,2)
              if (bc_klo_type(i,j,1) == MINF_ .or. &
                   bc_klo_type(i,j,1) == NSW_ .or. &
                   bc_klo_type(i,j,1) == FSW_ .or. &
                   bc_klo_type(i,j,1) == PSW_) then
                 d_t(i,j,k) = 0.0d0
              endif
           end do
        end do
      endif

      ! At up boundary
      if (shi(3) .gt. domhi(3) .and. hi(3).eq.ahi(3)) then
         k = ahi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               bcv = bc_khi_type(i,j,2)
               if(bc_khi_type(i,j,1) == MINF_ .or. &
                  bc_khi_type(i,j,1) == NSW_ .or. &
                  bc_khi_type(i,j,1) == FSW_ .or. &
                  bc_khi_type(i,j,1) == PSW_) then
                  d_t(i,j,k) = 0.0d0
               endif
            end do
         end do
      endif

   end subroutine calc_d_t

   end module calc_d_mod
