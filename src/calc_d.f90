module calc_d_mod

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   use param, only: zero, small_number

   ! Flag: Coupled DEM simulation
   use discretelement, only: DES_CONTINUUM_COUPLED
   use discretelement, only: DES_ONEWAY_COUPLED

   ! Pressure scale factor
   use scales, only: P_SCALE

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
   subroutine calc_d_e(slo, shi, ulo, uhi, alo, ahi, d_e, A_m, &
                       ep_g, f_gds, dx, dy, dz, domlo, domhi)

      use bc, only: cyclic_x

      integer, intent(in   ) :: slo(3),shi(3)
      integer, intent(in   ) :: ulo(3),uhi(3)
      integer, intent(in   ) :: alo(3),ahi(3)
      integer, intent(in   ) :: domlo(3),domhi(3)

      ! Pressure correction
      real(c_real), intent(  out) :: d_e&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      real(c_real), intent(in   ):: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3), -3:3)

      real(c_real), intent(in   ):: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: dx, dy, dz

      integer      :: i,j,k
      real(c_real) :: ayz, vol
      real(c_real) :: Am0, epga
      logical      :: coupled
      integer :: llo(3), lhi(3)

      COUPLED = (DES_CONTINUUM_COUPLED .AND. .NOT.DES_ONEWAY_COUPLED)

      ayz = dy*dz
      vol = dx*dy*dz

      llo = alo
      lhi = ahi

      ! if(.not.cyclic_x .and. alo(1) == domlo(1)) llo(1) = alo(1)+1

      do k = llo(3), lhi(3)
         do j = llo(2), lhi(2)
            do i = llo(1), lhi(1)
               Am0 = -A_m(i,j,k,0)
               if (abs(am0) > small_number) then
                  epga = ayz*0.5d0*(ep_g(i-1,j,k)+ep_g(i,j,k))
                  if(coupled) Am0 = Am0 + 0.5d0*vol* &
                     (f_gds(i-1,j,k) + f_gds(i,j,k))
                  d_e(i,j,k) = p_scale*epga/am0
               else
                  d_e(i,j,k) = zero
               endif
            enddo
         enddo
      enddo

   end subroutine calc_d_e

   subroutine calc_d_n(slo, shi, vlo, vhi, alo, ahi, d_n, A_m,&
                       ep_g, f_gds, dx, dy, dz, domlo, domhi)

      use bc, only: cyclic_y

      integer, intent(in   ) :: slo(3),shi(3)
      integer, intent(in   ) :: vlo(3),vhi(3)
      integer, intent(in   ) :: alo(3),ahi(3)
      integer, intent(in   ) :: domlo(3),domhi(3)

      ! Pressure correction
      real(c_real), intent(  out) :: d_n&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

      real(c_real), intent(in   ):: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3), -3:3)

      real(c_real), intent(in   ):: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: dx, dy, dz

      integer      :: i,j,k
      real(c_real) :: axz, vol
      real(c_real) :: Am0, epga
      logical      :: coupled
      integer :: llo(3), lhi(3)

      coupled = (des_continuum_coupled .and. .not.des_oneway_coupled)

      axz = dx*dz
      vol = dx*dy*dz

      llo = alo
      lhi = ahi

      do k = llo(3), lhi(3)
         do j = llo(2), lhi(2)
            do i = llo(1), lhi(1)
               Am0 = -A_m(i,j,k,0)
               if(abs(Am0) > small_number) then
                  epga = axz*0.5d0*(ep_g(i,j-1,k)+ep_g(i,j,k))
                  if(coupled) Am0 = Am0 + 0.5d0*vol* &
                     (f_gds(i,j-1,k) + f_gds(i,j,k))
                  d_n(i,j,k) = p_scale*epga/am0
               else
                  d_n(i,j,k) = zero
               endif

            enddo
         enddo
      enddo

   end subroutine calc_d_n

   subroutine calc_d_t(slo, shi, wlo, whi, alo, ahi, d_t, A_m,&
      ep_g, f_gds, dx, dy, dz, domlo, domhi)

      use bc, only: cyclic_z

      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) :: wlo(3),whi(3)
      integer     , intent(in   ) :: alo(3),ahi(3)
      integer, intent(in   ) :: domlo(3),domhi(3)

      ! Pressure correction
      real(c_real), intent(  out) :: d_t&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ):: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3), -3:3)

      real(c_real), intent(in   ):: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: dx, dy, dz

      integer      :: i,j,k
      real(c_real) :: axy, vol
      real(c_real) :: Am0, epga
      logical      :: coupled
      integer :: llo(3), lhi(3)

      coupled = (des_continuum_coupled .and. .not.des_oneway_coupled)

      axy = dx*dy
      vol = dx*dy*dz

      llo = alo
      lhi = ahi

      ! if(.not.cyclic_z .and. alo(3) == domlo(3)) llo(3) = alo(3)+1

      do k = llo(3), lhi(3)
        do j = llo(2), lhi(2)
           do i = llo(1), lhi(1)
              Am0 = -A_m(I,J,K,0)

              if (abs(Am0) > small_number) THEN

                 epga = axy*0.5d0*(ep_g(i,j,k-1)+ep_g(i,j,k))
                 if(coupled) Am0 = Am0 + 0.5d0*vol* &
                    (f_gds(i,j,k-1) + f_gds(i,j,k))
                 d_t(i,j,k) = p_scale*epga/am0

              else

                 d_t(i,j,k) = zero

              endif

           enddo
        enddo
     enddo

   end subroutine calc_d_t

   end module calc_d_mod
