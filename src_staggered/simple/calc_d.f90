module calc_d_mod

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   use param, only: zero, small_number

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
     subroutine calc_d_e(lo, hi, slo, shi, ulo, uhi, alo, ahi, d_e, A_m, ep_g, &
                         dy, dz, domlo, domhi, bc_ilo_type, bc_ihi_type, ng)

      integer, intent(in   ) ::  lo(3), hi(3), ng
      integer, intent(in   ) :: slo(3),shi(3)
      integer, intent(in   ) :: ulo(3),uhi(3)
      integer, intent(in   ) :: alo(3),ahi(3)
      integer, intent(in   ) :: domlo(3),domhi(3)

      ! Pressure correction
      real(rt), intent(  out) :: &
           d_e (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      real(rt), intent(in   ) :: &
           A_m (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3), -3:3), &
           ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer(c_int), intent(in   ) :: &
           bc_ilo_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_ihi_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2)

      real(rt), intent(in   ) :: dy, dz

      integer      :: i,j,k
      real(rt) :: const

      const = -p_scale*0.5d0*dy*dz

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               d_e(i,j,k) = const*(ep_g(i-1,j,k)+ep_g(i,j,k))/A_m(i,j,k,0)

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

   subroutine calc_d_n(lo, hi, slo, shi, vlo, vhi, alo, ahi, d_n, A_m, ep_g, &
                       dx, dz, domlo, domhi, bc_jlo_type, bc_jhi_type, ng)


      integer, intent(in   ) ::  lo(3), hi(3), ng
      integer, intent(in   ) :: slo(3),shi(3)
      integer, intent(in   ) :: vlo(3),vhi(3)
      integer, intent(in   ) :: alo(3),ahi(3)
      integer, intent(in   ) :: domlo(3),domhi(3)

      ! Pressure correction
      real(rt), intent(  out) :: &
           d_n (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

      real(rt), intent(in   ):: &
           A_m (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3), -3:3), &
           ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer(c_int), intent(in   ) :: &
           bc_jlo_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_jhi_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2)

      real(rt), intent(in   ) :: dx, dz

      integer      :: i,j,k
      real(rt) :: const

      const = -p_scale*0.5d0*dx*dz

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               d_n(i,j,k) = const*(ep_g(i,j-1,k)+ep_g(i,j,k))/A_m(i,j,k,0)

            enddo
         enddo
      enddo

      ! At bottom boundary
      if (slo(2) .lt. domlo(2) .and. lo(2).eq.alo(2)) then
         j = alo(2)
         do k = lo(3), hi(3)
            do i = lo(1), hi(1)
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

   subroutine calc_d_t(lo, hi, slo, shi, wlo, whi, alo, ahi, d_t, A_m, ep_g, &
                       dx, dy, domlo, domhi, bc_klo_type, bc_khi_type, ng)

      integer, intent(in   ) ::  lo(3), hi(3)
      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) :: wlo(3),whi(3)
      integer     , intent(in   ) :: alo(3),ahi(3)
      integer, intent(in   ) :: domlo(3),domhi(3), ng

      ! Pressure correction
      real(rt), intent(  out) :: &
           d_t (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(rt), intent(in   ):: &
           A_m (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3), -3:3), &
           ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer(c_int), intent(in   ) :: &
           bc_klo_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
           bc_khi_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

      real(rt), intent(in   ) :: dx, dy

      integer      :: i,j,k
      real(rt) :: const

      const = -p_scale*0.5d0*dx*dy

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              d_t(i,j,k) = const*(ep_g(i,j,k-1)+ep_g(i,j,k))/A_m(i,j,k,0)

           enddo
        enddo
     enddo

      ! At down boundary
     if (slo(3) .lt. domlo(3) .and. lo(3).eq.alo(3)) then
        k = alo(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
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
