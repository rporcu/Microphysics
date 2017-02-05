! ::: -----------------------------------------------------------
! ::: This routine is intended to be a generic fill function
! ::: for cell centered data.  It knows how to exrapolate,
! ::: and reflect data and can be used to suppliment problem
! ::: specific fill functions (ie. EXT_DIR).
! :::
! ::: INPUTS/OUTPUTS:
! ::: q        <=  array to fill
! ::: DIMS(q)   => index extent of q array
! ::: domlo,hi  => index extent of problem domain
! ::: dx        => cell spacing
! ::: xlo       => physical location of lower left hand
! :::            corner of q array
! ::: bc	=> array of boundary flags bc(SPACEDIM,lo:hi)
! :::
! ::: NOTE: corner data not used in computing soln but must have
! :::       reasonable values for arithmetic to live
! ::: -----------------------------------------------------------

      subroutine fill_bc(s,slo,shi,flag,vtype) &
         bind(C, name="fill_bc")

      use iso_c_binding , only: c_int
      use bl_fort_module, only: c_real

      use geometry      , only: domlo, domhi

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: vtype

      real(c_real), intent(inout) ::  s&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) ::  flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      integer    nlft, nrgt, nbot, ntop, nup, ndwn
      integer    ilo, ihi, jlo, jhi, klo, khi
      integer    is,  ie,  js,  je,  ks,  ke
      integer    i, j, k

      integer, parameter:: REFLECT_ODD  = -1
      ! integer, parameter:: INT_DIR      =  0
      integer, parameter:: REFLECT_EVEN =  1
      integer, parameter:: FOEXTRAP     =  2
      ! integer, parameter:: EXT_DIR      =  3
      integer, parameter:: HOEXTRAP     =  4

      real(c_real), parameter:: half         = 0.5d0


      integer, parameter :: bc(3,2) = FOEXTRAP

      is = max(slo(1),domlo(1))
      ie = min(shi(1),domhi(1))
      js = max(slo(2),domlo(2))
      je = min(shi(2),domhi(2))
      ks = max(slo(3),domlo(3))
      ke = min(shi(3),domhi(3))

      nlft = max(0,domlo(1)-slo(1))
      nbot = max(0,domlo(2)-slo(2))
      ndwn = max(0,domlo(3)-slo(3))

      nrgt = max(0,shi(1)-domhi(1))
      ntop = max(0,shi(2)-domhi(2))
      nup  = max(0,shi(3)-domhi(3))



      if(vtype /= 0) return

!
!     ::::: first fill sides
!
      if (nlft .gt. 0) then
         ilo = domlo(1)

         if (bc(1,1) .eq. FOEXTRAP) then
            do i = 1, nlft
               s(ilo-i,:,:) = s(ilo,:,:)
            end do
         else if (bc(1,1) .eq. HOEXTRAP) then
            do i = 2, nlft
               s(ilo-i,:,:) = s(ilo,:,:)
            end do
            if (ilo+2 .le. ie) then
                  s(ilo-1,:,:) = (15*s(ilo,:,:) - 10*s(ilo+1,:,:) + &
                       3*s(ilo+2,:,:)) * 0.125d0
            else
                  s(ilo-1,:,:) = half*(3*s(ilo,:,:) - s(ilo+1,:,:))
            end if
         else if (bc(1,1) .eq. REFLECT_EVEN) then
            do i = 1, nlft
               s(ilo-i,:,:) = s(ilo+i-1,:,:)
            end do
         else if (bc(1,1) .eq. REFLECT_ODD) then
            do i = 1, nlft
               s(ilo-i,:,:) = -s(ilo+i-1,:,:)
            end do
         end if
      end if

      if (nrgt .gt. 0) then
         ihi = domhi(1)

         if (bc(1,2) .eq. FOEXTRAP) then
            do i = 1, nrgt
               s(ihi+i,:,:) = s(ihi,:,:)
            end do
         else if (bc(1,2) .eq. HOEXTRAP) then
            do i = 2, nrgt
               s(ihi+i,:,:) = s(ihi,:,:)
            end do
            if (ihi-2 .ge. is) then
                  s(ihi+1,:,:) = (15*s(ihi,:,:) - 10*s(ihi-1,:,:) +  &
                       3*s(ihi-2,:,:)) * 0.125d0
            else
                  s(ihi+1,:,:) = half*(3*s(ihi,:,:) - s(ihi-1,:,:))
            end if
         else if (bc(1,2) .eq. REFLECT_EVEN) then
            do i = 1, nrgt
                  s(ihi+i,:,:) = s(ihi-i+1,:,:)
            end do
        else if (bc(1,2) .eq. REFLECT_ODD) then
            do i = 1, nrgt
               s(ihi+i,:,:) = -s(ihi-i+1,:,:)
            end do
        end if
      end if

      if (nbot .gt. 0) then
         jlo = domlo(2)

         if (bc(2,1) .eq. FOEXTRAP) then
!     THIS IS A HACK SO WE DONT OVERWRITE INFLOW BC'S IN DEM06
            ! do j = 1, nbot
            !    s(:,jlo-j,:) = s(:,jlo,:)
            ! end do
         else if (bc(2,1) .eq. HOEXTRAP) then
            do j = 2, nbot
               s(:,jlo-j,:) = s(:,jlo,:)
            end do
            if (jlo+2 .le. je) then
               s(:,jlo-1,:) = (15.*s(:,jlo,:) - 10.*s(:,jlo+1,:) + &
                                3.*s(:,jlo+2,:)) * 0.125d0
            else
               s(:,jlo-1,:) = half*(3*s(:,jlo,:) - s(:,jlo+1,:))
            end if
         else if (bc(2,1) .eq. REFLECT_EVEN) then
            do j = 1, nbot
               s(:,jlo-j,:) = s(:,jlo+j-1,:)
            end do
         else if (bc(2,1) .eq. REFLECT_ODD) then
            do j = 1, nbot
               s(:,jlo-j,:) = -s(:,jlo+j-1,:)
            end do
         end if
      end if

      if (ntop .gt. 0) then
         jhi = domhi(2)

         if (bc(2,2) .eq. FOEXTRAP) then
!     THIS IS A HACK SO WE DONT OVERWRITE OUTFLOW BC'S IN DEM06
            ! do j = 1, ntop
            !    s(:,jhi+j,:) = s(:,jhi,:)
            ! end do
         else if (bc(2,2) .eq. HOEXTRAP) then
            do j = 2, ntop
               s(:,jhi+j,:) = s(:,jhi,:)
            end do
            if (jhi-2 .ge. js) then
                  s(:,jhi+1,:) = (15*s(:,jhi,:) - 10*s(:,jhi-1,:) + &
                       3*s(:,jhi-2,:)) * 0.125d0
            else
                  s(:,jhi+1,:) = half*(3*s(:,jhi,:) - s(:,jhi-1,:))
            end if
         else if (bc(2,2) .eq. REFLECT_EVEN) then
            do j = 1, ntop
               s(:,jhi+j,:) = s(:,jhi-j+1,:)
            end do
         else if (bc(2,2) .eq. REFLECT_ODD) then
            do j = 1, ntop
               s(:,jhi+j,:) = -s(:,jhi-j+1,:)
            end do
         end if
      end if

      if (ndwn .gt. 0) then
         klo = domlo(3)

         if (bc(3,1) .eq. FOEXTRAP) then
            do k = 1, ndwn
               s(:,:,klo-k) = s(:,:,klo)
            end do
         else if (bc(3,1) .eq. HOEXTRAP) then
            do k = 2, ndwn
               s(:,:,klo-k) = s(:,:,klo)
            end do
            if (klo+2 .le. ke) then
                  s(:,:,klo-1) = (15*s(:,:,klo) - 10*s(:,:,klo+1) + &
                       3*s(:,:,klo+2)) * 0.125d0
         else
                  s(:,:,klo-1) = half*(3*s(:,:,klo) - s(:,:,klo+1))
            end if
         else if (bc(3,1) .eq. REFLECT_EVEN) then
            do k = 1, ndwn
               s(:,:,klo-k) = s(:,:,klo+k-1)
            end do
         else if (bc(3,1) .eq. REFLECT_ODD) then
            do k = 1, ndwn
               s(:,:,klo-k) = -s(:,:,klo+k-1)
            end do
         end if
      end if

      if (nup .gt. 0) then
         khi = domhi(3)

         if (bc(3,2) .eq. FOEXTRAP) then
            do k = 1, nup
               s(:,:,khi+k) = s(:,:,khi)
            end do
         else if (bc(3,2) .eq. HOEXTRAP) then
            do k = 2, nup
               s(:,:,khi+k) = s(:,:,khi)
            end do
            if (khi-2 .ge. ks) then
               s(:,:,khi+1) = (15.*s(:,:,khi) - 10*s(:,:,khi-1) + &
                                3.*s(:,:,khi-2)) * 0.125d0
            else
               s(:,:,khi+1) = half*(3*s(:,:,khi) - s(:,:,khi-1))
            end if
         else if (bc(3,2) .eq. REFLECT_EVEN) then
            do k = 1, nup
               s(:,:,khi+k) = s(:,:,khi-k+1)
            end do
         else if (bc(3,2) .eq. REFLECT_ODD) then
            do k = 1, nup
               s(:,:,khi+k) = -s(:,:,khi-k+1)
            end do
         end if
      end if


!
!    First correct the i-j edges and all corners
!
      if((nlft .gt. 0 .and. bc(1,1) .eq. HOEXTRAP) .and. &
         (nbot .gt. 0 .and. bc(2,1) .eq. HOEXTRAP) ) then
         if (jlo+2 .le. je) then
            do k = slo(3), shi(3)
               s(ilo-1,jlo-1,k) = half * 0.125d0 * &
                    (15*s(ilo-1,jlo,k) - 10*s(ilo-1,jlo+1,k) + &
                    3*s(ilo-1,jlo+2,k))
            end do
         else
            do k = slo(3), shi(3)
               s(ilo-1,jlo-1,k) = half * half * &
                    (3*s(ilo-1,jlo,k) - s(ilo-1,jlo+1,k))
            end do
         end if

         if (ilo+2 .le. ie) then
            do k = slo(3), shi(3)
               s(ilo-1,jlo-1,k) = s(ilo-1,jlo-1,k) + half * 0.125d0 * &
                    (15*s(ilo,jlo-1,k) - 10*s(ilo+1,jlo-1,k) + &
                    3*s(ilo+2,jlo-1,k))
            end do
         else
            do k = slo(3), shi(3)
               s(ilo-1,jlo-1,k) = s(ilo-1,jlo-1,k) + half * half * &
                    (3*s(ilo,jlo-1,k) - s(ilo+1,jlo-1,k))
            end do
         end if

         if (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) then
            if (klo+2 .le. ke) then
               s(ilo-1,jlo-1,klo-1) = 0.125d0 * ( &
                    (15*s(ilo-1,jlo-1,klo) - 10*s(ilo-1,jlo-1,klo+1) + &
                    3*s(ilo-1,jlo-1,klo+2)) )
            else
               s(ilo-1,jlo-1,klo-1) = half *  &
                    (3*s(ilo-1,jlo-1,klo) - s(ilo-1,jlo-1,klo+1))
            end if
         end if

         if (nup .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) then
            if (khi-2 .ge. ks) then
               s(ilo-1,jlo-1,khi+1) = 0.125d0 * ( &
                    (15*s(ilo-1,jlo-1,khi) - 10*s(ilo-1,jlo-1,khi-1) + &
                    3*s(ilo-1,jlo-1,khi-2)) )
            else
               s(ilo-1,jlo-1,khi+1) = half *  &
                    (3*s(ilo-1,jlo-1,khi) - s(ilo-1,jlo-1,khi-1))
            end if
         end if

      end if
!
! ****************************************************************************
!
      if ((nlft .gt. 0 .and. bc(1,1) .eq. HOEXTRAP) .and. &
          (ntop .gt. 0 .and. bc(2,2) .eq. HOEXTRAP) ) then
         if (jhi-2 .ge. js) then
            do k = slo(3), shi(3)
               s(ilo-1,jhi+1,k) = half * 0.125d0 * &
                    (15*s(ilo-1,jhi,k) - 10*s(ilo-1,jhi-1,k) + &
                    3*s(ilo-1,jhi-2,k))
            end do
         else
            do k = slo(3), shi(3)
               s(ilo-1,jhi+1,k) = half * half * &
                    (3*s(ilo-1,jhi,k) - s(ilo-1,jhi-1,k))
            end do
         end if

         if (ilo+2 .le. ie) then
            do k = slo(3), shi(3)
               s(ilo-1,jhi+1,k) = s(ilo-1,jhi+1,k) + half * 0.125d0 * &
                    (15*s(ilo,jhi+1,k) - 10*s(ilo+1,jhi+1,k) + &
                    3*s(ilo+2,jhi+1,k))
            end do
         else
            do k = slo(3), shi(3)
               s(ilo-1,jhi+1,k) = s(ilo-1,jhi+1,k) + half * half * &
                    (3*s(ilo,jhi+1,k) - s(ilo+1,jhi+1,k))
            end do
         end if

         if (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) then
            if (klo+2 .le. ke) then
               s(ilo-1,jhi+1,klo-1) = 0.125d0 * ( &
                    (15*s(ilo-1,jhi+1,klo) - 10*s(ilo-1,jhi+1,klo+1) + &
                    3*s(ilo-1,jhi+1,klo+2)) )
            else
               s(ilo-1,jhi+1,klo-1) = half * &
                    (3*s(ilo-1,jhi+1,klo) - s(ilo-1,jhi+1,klo+1))
            end if
         end if

         if (nup .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) then
            if (khi-2 .ge. ks) then
               s(ilo-1,jhi+1,khi+1) = 0.125d0 * ( &
                    (15*s(ilo-1,jhi+1,khi) - 10*s(ilo-1,jhi+1,khi-1) + &
                    3*s(ilo-1,jhi+1,khi-2)) )
            else
               s(ilo-1,jhi+1,khi+1) = half * &
                    (3*s(ilo-1,jhi+1,khi) - s(ilo-1,jhi+1,khi-1))
            end if
         end if

      end if
!
! ****************************************************************************
!
      if ((nrgt .gt. 0 .and. bc(1,2) .eq. HOEXTRAP) .and. &
          (nbot .gt. 0 .and. bc(2,1) .eq. HOEXTRAP) ) then
         if (jlo+2 .le. je) then
            do k = slo(3), shi(3)
               s(ihi+1,jlo-1,k) = half * 0.125d0 * &
                    (15*s(ihi+1,jlo,k) - 10*s(ihi+1,jlo+1,k) + &
                    3*s(ihi+1,jlo+2,k))
            end do
         else
            do k = slo(3), shi(3)
               s(ihi+1,jlo-1,k) = half * half * &
                    (3*s(ihi+1,jlo,k) - s(ihi+1,jlo+1,k))
            end do
         end if

         if (ihi-2 .ge. is) then
            do k = slo(3), shi(3)
               s(ihi+1,jlo-1,k) = s(ihi+1,jlo-1,k) + half * 0.125d0 * &
                    (15*s(ihi,jlo-1,k) - 10*s(ihi-1,jlo-1,k) + &
                    3*s(ihi-2,jlo-1,k))
            end do
         else
            do k = slo(3), shi(3)
               s(ihi+1,jlo-1,k) = s(ihi+1,jlo-1,k) + half * half * &
                    (3*s(ihi,jlo-1,k) - s(ihi-1,jlo-1,k))
            end do
         end if

         if (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) then
            if (klo+2 .le. ke) then
               s(ihi+1,jlo-1,klo-1) = 0.125d0 * &
                    (15*s(ihi+1,jlo-1,klo) - 10*s(ihi+1,jlo-1,klo+1) + &
                    3*s(ihi+1,jlo-1,klo+2))
            else
               s(ihi+1,jlo-1,klo-1) = half * &
                    (3*s(ihi+1,jlo-1,klo) - s(ihi+1,jlo-1,klo+1))
            end if
         end if

         if (nup .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) then
            if (khi-2 .ge. ks) then
               s(ihi+1,jlo-1,khi+1) = 0.125d0 * &
                    (15*s(ihi+1,jlo-1,khi) - 10*s(ihi+1,jlo-1,khi-1) + &
                    3*s(ihi+1,jlo-1,khi-2))
            else
               s(ihi+1,jlo-1,khi+1) = half * &
                    (3*s(ihi+1,jlo-1,khi) - s(ihi+1,jlo-1,khi-1))
            end if
         end if

      end if
!
! ****************************************************************************
!
      if ((nrgt .gt. 0 .and. bc(1,2) .eq. HOEXTRAP) .and. &
          (ntop .gt. 0 .and. bc(2,2) .eq. HOEXTRAP) ) then
         if (jhi-2 .ge. js) then
            do k = slo(3), shi(3)
               s(ihi+1,jhi+1,k) = half * 0.125d0 * &
                    (15*s(ihi+1,jhi,k) - 10*s(ihi+1,jhi-1,k) + &
                    3*s(ihi+1,jhi-2,k))
            end do
         else
            do k = slo(3), shi(3)
               s(ihi+1,jhi+1,k) = half * half * &
                    (3*s(ihi+1,jhi,k) - s(ihi+1,jhi-1,k))
            end do
         end if

         if (ihi-2 .ge. is) then
            do k = slo(3), shi(3)
               s(ihi+1,jhi+1,k) = s(ihi+1,jhi+1,k) + half * 0.125d0 * &
                    (15*s(ihi,jhi+1,k) - 10*s(ihi-1,jhi+1,k) + &
                    3*s(ihi-2,jhi+1,k))
            end do
         else
            do k = slo(3), shi(3)
               s(ihi+1,jhi+1,k) = s(ihi+1,jhi+1,k) + half * half * &
                    (3*s(ihi,jhi+1,k) - s(ihi-1,jhi+1,k))
            end do
         end if

         if (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) then
            if (klo+2 .le. ke) then
               s(ihi+1,jhi+1,klo-1) = 0.125d0 * &
                    (15*s(ihi+1,jhi+1,klo) - 10*s(ihi+1,jhi+1,klo+1) + &
                    3*s(ihi+1,jhi+1,klo+2))
            else
               s(ihi+1,jhi+1,klo-1) = half * &
                    (3*s(ihi+1,jhi+1,klo) - s(ihi+1,jhi+1,klo+1))
            end if
         end if

         if (nup .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) then
            if (khi-2 .ge. ks) then
               s(ihi+1,jhi+1,khi+1) = 0.125d0 * &
                    (15*s(ihi+1,jhi+1,khi) - 10*s(ihi+1,jhi+1,khi-1) + &
                    3*s(ihi+1,jhi+1,khi-2))
            else
               s(ihi+1,jhi+1,khi+1) = half * &
                    (3*s(ihi+1,jhi+1,khi) - s(ihi+1,jhi+1,khi-1))
            end if
         end if

      end if
!
!    Next correct the i-k edges
!
      if ((nlft .gt. 0 .and. bc(1,1) .eq. HOEXTRAP) .and. &
          (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) ) then
         if (klo+2 .le. ke) then
            do j = slo(2), shi(2)
               s(ilo-1,j,klo-1) = half * 0.125d0 * &
                    (15*s(ilo-1,j,klo) - 10*s(ilo-1,j,klo+1) +  &
                    3*s(ilo-1,j,klo+2))
            end do
         else
            do j = slo(2), shi(2)
               s(ilo-1,j,klo-1) = half * half * &
                    (3*s(ilo-1,j,klo) - s(ilo-1,j,klo+1))
            end do
         end if

         if (ilo+2 .le. ie) then
            do j = slo(2), shi(2)
               s(ilo-1,j,klo-1) = s(ilo-1,j,klo-1) + half * 0.125d0 * &
                    (15*s(ilo,j,klo-1) - 10*s(ilo+1,j,klo-1) + &
                    3*s(ilo+2,j,klo-1))
            end do
         else
            do j = slo(2), shi(2)
               s(ilo-1,j,klo-1) = s(ilo-1,j,klo-1) + half * half * &
                    (3*s(ilo,j,klo-1) - s(ilo+1,j,klo-1))
            end do
         end if
      end if
!
! ****************************************************************************
!
      if ((nlft .gt. 0 .and. bc(1,1) .eq. HOEXTRAP) .and. &
          (nup  .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) ) then
         if (khi-2 .ge. ks) then
            do j = slo(2), shi(2)
               s(ilo-1,j,khi+1) = half * 0.125d0 * &
                    (15*s(ilo-1,j,khi) - 10*s(ilo-1,j,khi-1) + &
                    3*s(ilo-1,j,khi-2))
            end do
         else
            do j = slo(2), shi(2)
               s(ilo-1,j,khi+1) = half * half * &
                    (3*s(ilo-1,j,khi) - s(ilo-1,j,khi-1))
            end do
         end if

         if (ilo+2 .le. ie) then
            do j = slo(2), shi(2)
               s(ilo-1,j,khi+1) = s(ilo-1,j,khi+1) + half * 0.125d0 * &
                    (15*s(ilo,j,khi+1) - 10*s(ilo+1,j,khi+1) + &
                    3*s(ilo+2,j,khi+1))
            end do
         else
            do j = slo(2), shi(2)
               s(ilo-1,j,khi+1) = s(ilo-1,j,khi+1) + half * half * &
                    (3*s(ilo,j,khi+1) - s(ilo+1,j,khi+1))
            end do
         end if
      end if
!
! ****************************************************************************
!
      if ((nrgt .gt. 0 .and. bc(1,2) .eq. HOEXTRAP) .and. &
          (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) ) then
         if (klo+2 .le. ke) then
            do j = slo(2), shi(2)
               s(ihi+1,j,klo-1) = half * 0.125d0 * &
                    (15*s(ihi+1,j,klo) - 10*s(ihi+1,j,klo+1) +  &
                    3*s(ihi+1,j,klo+2))
            end do
         else
            do j = slo(2), shi(2)
               s(ihi+1,j,klo-1) = half * half *  &
                    (3*s(ihi+1,j,klo) - s(ihi+1,j,klo+1))
            end do
         end if

         if (ihi-2 .ge. is) then
            do j = slo(2), shi(2)
               s(ihi+1,j,klo-1) = s(ihi+1,j,klo-1) + half * 0.125d0 * &
                    (15*s(ihi,j,klo-1) - 10*s(ihi-1,j,klo-1) +  &
                    3*s(ihi-2,j,klo-1))
            end do
         else
            do j = slo(2), shi(2)
               s(ihi+1,j,klo-1) = s(ihi+1,j,klo-1) + half * half * &
                    (3*s(ihi,j,klo-1) - s(ihi-1,j,klo-1))
            end do
         end if
      end if
!
! ****************************************************************************
!
      if ((nrgt .gt. 0 .and. bc(1,2) .eq. HOEXTRAP) .and. &
          (nup  .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) ) then
         if (khi-2 .ge. ks) then
            do j = slo(2), shi(2)
               s(ihi+1,j,khi+1) = half * 0.125d0 *  &
                    (15*s(ihi+1,j,khi) - 10*s(ihi+1,j,khi-1) +  &
                    3*s(ihi+1,j,khi-2))
            end do
         else
            do j = slo(2), shi(2)
               s(ihi+1,j,khi+1) = half * half *  &
                    (3*s(ihi+1,j,khi) - s(ihi+1,j,khi-1))
            end do
         end if
         if (ihi-2 .ge. is) then
            do j = slo(2), shi(2)
               s(ihi+1,j,khi+1) = s(ihi+1,j,khi+1) + half * 0.125d0 *  &
                    (15*s(ihi,j,khi+1) - 10*s(ihi-1,j,khi+1) +  &
                    3*s(ihi-2,j,khi+1))
            end do
         else
            do j = slo(2), shi(2)
               s(ihi+1,j,khi+1) = s(ihi+1,j,khi+1) + half * half *  &
                    (3*s(ihi,j,khi+1) - s(ihi-1,j,khi+1))
            end do
         end if
      end if
!
!    Next correct the j-k edges
!
      if ((nbot .gt. 0 .and. bc(2,1) .eq. HOEXTRAP) .and. &
          (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) ) then
         if (klo+2 .le. ke) then
            do i = slo(1), shi(1)
               s(i,jlo-1,klo-1) = half * 0.125d0 * &
                    (15*s(i,jlo-1,klo) - 10*s(i,jlo-1,klo+1) +  &
                    3*s(i,jlo-1,klo+2))
            end do
         else
            do i = slo(1), shi(1)
               s(i,jlo-1,klo-1) = half * half *  &
                    (3*s(i,jlo-1,klo) - s(i,jlo-1,klo+1))
            end do
         end if
         if (jlo+2 .le. je) then
            do i = slo(1), shi(1)
               s(i,jlo-1,klo-1) = s(i,jlo-1,klo-1) + half * 0.125d0 *  &
                    (15*s(i,jlo,klo-1) - 10*s(i,jlo+1,klo-1) +  &
                    3*s(i,jlo+2,klo-1))
            end do
         else
            do i = slo(1), shi(1)
               s(i,jlo-1,klo-1) = s(i,jlo-1,klo-1) + half * half *  &
                    (3*s(i,jlo,klo-1) - s(i,jlo+1,klo-1))
            end do
         end if
      end if
!
! ****************************************************************************
!
      if ((nbot .gt. 0 .and. bc(2,1) .eq. HOEXTRAP) .and. &
          (nup  .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) ) then
         if (khi-2 .ge. ks) then
            do i = slo(1), shi(1)
               s(i,jlo-1,khi+1) = half * 0.125d0 * &
                    (15*s(i,jlo-1,khi) - 10*s(i,jlo-1,khi-1) + &
                    3*s(i,jlo-1,khi-2))
            end do
         else
            do i = slo(1), shi(1)
               s(i,jlo-1,khi+1) = half * half * &
                    (3*s(i,jlo-1,khi) - s(i,jlo-1,khi-1))
            end do
         end if

         if (jlo+2 .le. je) then
            do i = slo(1), shi(1)
               s(i,jlo-1,khi+1) = s(i,jlo-1,khi+1) + half * 0.125d0 * &
                    (15*s(i,jlo,khi+1) - 10*s(i,jlo+1,khi+1) + &
                    3*s(i,jlo+2,khi+1))
            end do
         else
            do i = slo(1), shi(1)
               s(i,jlo-1,khi+1) = s(i,jlo-1,khi+1) + half * half * &
                    (3*s(i,jlo,khi+1) - s(i,jlo+1,khi+1))
            end do
         end if
      end if
!
! ****************************************************************************
!
      if ((ntop .gt. 0 .and. bc(2,2) .eq. HOEXTRAP) .and. &
          (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) ) then
         if (klo+2 .le. ke) then
            do i = slo(1), shi(1)
               s(i,jhi+1,klo-1) = half * 0.125d0 * &
                    (15*s(i,jhi+1,klo) - 10*s(i,jhi+1,klo+1) + &
                    3*s(i,jhi+1,klo+2))
            end do
         else
            do i = slo(1), shi(1)
               s(i,jhi+1,klo-1) = half * half *  &
                    (3*s(i,jhi+1,klo) - s(i,jhi+1,klo+1))
            end do
         end if
         if (jhi-2 .ge. js) then
            do i = slo(1), shi(1)
               s(i,jhi+1,klo-1) = s(i,jhi+1,klo-1) + half * 0.125d0 *  &
                    (15*s(i,jhi,klo-1) - 10*s(i,jhi-1,klo-1) +  &
                    3*s(i,jhi-2,klo-1))
            end do
         else
            do i = slo(1), shi(1)
               s(i,jhi+1,klo-1) = s(i,jhi+1,klo-1) + half * half *  &
                    (3*s(i,jhi,klo-1) - s(i,jhi-1,klo-1))
            end do
         end if
      end if

!
! ****************************************************************************
!
      if ((ntop .gt. 0 .and. bc(2,2) .eq. HOEXTRAP) .and. &
          (nup  .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) ) then
         if (khi-2 .ge. ks) then
            do i = slo(1), shi(1)
               s(i,jhi+1,khi+1) = half * 0.125d0 * &
                    (15*s(i,jhi+1,khi) - 10*s(i,jhi+1,khi-1) + &
                    3*s(i,jhi+1,khi-2))
            end do
         else
            do i = slo(1), shi(1)
               s(i,jhi+1,khi+1) = half * half * &
                    (3*s(i,jhi+1,khi) - s(i,jhi+1,khi-1))
            end do
         end if
         if (jhi-2 .ge. js) then
            do i = slo(1), shi(1)
               s(i,jhi+1,khi+1) = s(i,jhi+1,khi+1) + half * 0.125d0 * &
                    (15*s(i,jhi,khi+1) - 10*s(i,jhi-1,khi+1) + &
                    3*s(i,jhi-2,khi+1))
            end do
         else
            do i = slo(1), shi(1)
               s(i,jhi+1,khi+1) = s(i,jhi+1,khi+1) + half * half * &
                    (3*s(i,jhi,khi+1) - s(i,jhi-1,khi+1))
            end do
         end if
      end if

    end subroutine fill_bc

!
! ****************************************************************************
! ****************************************************************************
! ****************************************************************************
!
    subroutine hoextraptocc(s,slo,shi)

      use iso_c_binding , only: c_int
      use bl_fort_module, only: c_real

      use geometry      , only: domlo, domhi

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)

      real(c_real), intent(inout) ::  s&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), parameter:: half         = 0.5d0

      integer    nlft, nrgt, nbot, ntop, nup, ndwn
      integer    ilo, ihi, jlo, jhi, klo, khi
      integer    is,  ie,  js,  je,  ks,  ke
      integer    i, j, k

      is = max(slo(1),domlo(1))
      ie = min(slo(1),domhi(1))
      js = max(slo(1),domlo(2))
      je = min(slo(1),domhi(2))
      ks = max(slo(1),domlo(3))
      ke = min(shi(3),domhi(3))

      nlft = max(0,domlo(1)-slo(1))
      nbot = max(0,domlo(2)-slo(2))
      ndwn = max(0,domlo(3)-slo(3))

      nrgt = max(0,shi(1)-domhi(1))
      ntop = max(0,shi(2)-domhi(2))
      nup  = max(0,shi(3)-domhi(3))

!
!     First fill sides.
!
      if (nlft .gt. 0) then
         ilo = domlo(1)

         do i = 2, nlft
            do k = slo(3), shi(3)
               do j = slo(2), shi(2)
                  s(ilo-i,j,k) = s(ilo,j,k)
               end do
            end do
         end do
         if (ilo+2 .le. ie) then
            do k = slo(3), shi(3)
               do j = slo(2), shi(2)
                  s(ilo-1,j,k) = 3*s(ilo,j,k) - 3*s(ilo+1,j,k) + &
                                 s(ilo+2,j,k)
               end do
            end do
         else
            do k = slo(3), shi(3)
               do j = slo(2), shi(2)
                  s(ilo-1,j,k) = 2*s(ilo,j,k) - s(ilo+1,j,k)
               end do
            end do
         end if
      end if

      if (nrgt .gt. 0) then
         ihi = domhi(1)

         do i = 2, nrgt
            do k = slo(3), shi(3)
               do j = slo(2), shi(2)
                  s(ihi+i,j,k) = s(ihi,j,k)
               end do
            end do
         end do
         if (ihi-2 .ge. is) then
            do k = slo(3), shi(3)
               do j = slo(2), shi(2)
                  s(ihi+1,j,k) = 3*s(ihi,j,k) - 3*s(ihi-1,j,k) + &
                                 s(ihi-2,j,k)
               end do
            end do
         else
            do k = slo(3), shi(3)
               do j = slo(2), shi(2)
                  s(ihi+1,j,k) = 2*s(ihi,j,k) - s(ihi-1,j,k)
               end do
            end do
         end if
      end if

      if (nbot .gt. 0) then
         jlo = domlo(2)

         do j = 2, nbot
            do k = slo(3), shi(3)
               do i = slo(1), shi(1)
                  s(i,jlo-j,k) = s(i,jlo,k)
               end do
            end do
         end do
         if (jlo+2 .le. je) then
            do k = slo(3), shi(3)
               do i = slo(1), shi(1)
                  s(i,jlo-1,k) = 3*s(i,jlo,k) - 3*s(i,jlo+1,k) + &
                                 s(i,jlo+2,k)
               end do
            end do
         else
            do k = slo(3), shi(3)
               do i = slo(1), shi(1)
                  s(i,jlo-1,k) = 2*s(i,jlo,k) - s(i,jlo+1,k)
               end do
            end do
         end if
      end if

      if (ntop .gt. 0) then
         jhi = domhi(2)

         do j = 2, ntop
            do k = slo(3), shi(3)
               do i = slo(1), shi(1)
                  s(i,jhi+j,k) = s(i,jhi,k)
               end do
            end do
         end do
         if (jhi-2 .ge. js) then
            do k = slo(3), shi(3)
               do i = slo(1), shi(1)
                  s(i,jhi+1,k) = 3*s(i,jhi,k) - 3*s(i,jhi-1,k) + &
                                 s(i,jhi-2,k)
               end do
            end do
         else
            do k = slo(3), shi(3)
               do i = slo(1), shi(1)
                  s(i,jhi+1,k) = 2*s(i,jhi,k) - s(i,jhi-1,k)
               end do
            end do
         end if
      end if

      if (ndwn .gt. 0) then
         klo = domlo(3)

         do k = 2, ndwn
            do j = slo(2), shi(2)
               do i = slo(1), shi(1)
                  s(i,j,klo-k) = s(i,j,klo)
               end do
            end do
         end do
         if (klo+2 .le. ke) then
            do j = slo(2), shi(2)
               do i = slo(1), shi(1)
                  s(i,j,klo-1) = 3*s(i,j,klo) - 3*s(i,j,klo+1) + &
                                 s(i,j,klo+2)
               end do
            end do
         else
            do j = slo(2), shi(2)
               do i = slo(1), shi(1)
                  s(i,j,klo-1) = 2*s(i,j,klo) - s(i,j,klo+1)
               end do
            end do
         end if
      end if

      if (nup .gt. 0) then
         khi = domhi(3)

         do k = 2, nup
            do j = slo(2), shi(2)
               do i = slo(1), shi(1)
                  s(i,j,khi+k) = s(i,j,khi)
               end do
            end do
         end do
         if (khi-2 .ge. ks) then
            do j = slo(2), shi(2)
               do i = slo(1), shi(1)
                  s(i,j,khi+1) = 3*s(i,j,khi) - 3*s(i,j,khi-1) + &
                                 s(i,j,khi-2)
               end do
            end do
         else
            do j = slo(2), shi(2)
               do i = slo(1), shi(1)
                  s(i,j,khi+1) = 2*s(i,j,khi) - s(i,j,khi-1)
               end do
            end do
         end if
      end if
!
!    First correct the i-j edges and all corners
!
      if ((nlft .gt. 0) .and. (nbot .gt. 0)) then
         if (jlo+2 .le. je) then
            do k = slo(3), shi(3)
               s(ilo-1,jlo-1,k) = half * &
                    (3*s(ilo-1,jlo,k) - 3*s(ilo-1,jlo+1,k) + &
                    s(ilo-1,jlo+2,k))
            end do
         else
            do k = slo(3), shi(3)
               s(ilo-1,jlo-1,k) = half * &
                    (2*s(ilo-1,jlo,k) - s(ilo-1,jlo+1,k))
            end do
         end if

         if (ilo+2 .le. ie) then
            do k = slo(3), shi(3)
               s(ilo-1,jlo-1,k) = s(ilo-1,jlo-1,k) + half * &
                    (3*s(ilo,jlo-1,k) - 3*s(ilo+1,jlo-1,k) + &
                    s(ilo+2,jlo-1,k))
            end do
         else
            do k = slo(3), shi(3)
               s(ilo-1,jlo-1,k) = s(ilo-1,jlo-1,k) + half * &
                    (2*s(ilo,jlo-1,k) - s(ilo+1,jlo-1,k))
            end do
         end if

         if (ndwn .gt. 0) then
            if (klo+2 .le. ke) then
               s(ilo-1,jlo-1,klo-1) = &
                    (3*s(ilo-1,jlo-1,klo) - 3*s(ilo-1,jlo-1,klo+1) +&
                    s(ilo-1,jlo-1,klo+2))
            else
               s(ilo-1,jlo-1,klo-1) = &
                    (2*s(ilo-1,jlo-1,klo) - s(ilo-1,jlo-1,klo+1))
            end if
         end if

         if (nup .gt. 0) then
            if (khi-2 .ge. ks) then
               s(ilo-1,jlo-1,khi+1) = &
                    (3*s(ilo-1,jlo-1,khi) - 3*s(ilo-1,jlo-1,khi-1) +&
                    s(ilo-1,jlo-1,khi-2))
            else
               s(ilo-1,jlo-1,khi+1) = &
                    (2*s(ilo-1,jlo-1,khi) - s(ilo-1,jlo-1,khi-1))
            end if
         end if

      end if
!
! ****************************************************************************
!
      if ((nlft .gt. 0) .and. (ntop .gt. 0)) then
         if (jhi-2 .ge. js) then
            do k = slo(3), shi(3)
               s(ilo-1,jhi+1,k) = half * &
                    (3*s(ilo-1,jhi,k) - 3*s(ilo-1,jhi-1,k) +  &
                    s(ilo-1,jhi-2,k))
            end do
         else
            do k = slo(3), shi(3)
               s(ilo-1,jhi+1,k) = half * &
                    (2*s(ilo-1,jhi,k) - s(ilo-1,jhi-1,k))
            end do
         end if

         if (ilo+2 .le. ie) then
            do k = slo(3), shi(3)
               s(ilo-1,jhi+1,k) = s(ilo-1,jhi+1,k) + half * &
                    (3*s(ilo,jhi+1,k) - 3*s(ilo+1,jhi+1,k) +  &
                    s(ilo+2,jhi+1,k))
            end do
         else
            do k = slo(3), shi(3)
               s(ilo-1,jhi+1,k) = s(ilo-1,jhi+1,k) + half * &
                    (2*s(ilo,jhi+1,k) - s(ilo+1,jhi+1,k))
            end do
         end if

         if (ndwn .gt. 0) then
            if (klo+2 .le. ke) then
               s(ilo-1,jhi+1,klo-1) = &
                    (3*s(ilo-1,jhi+1,klo) - 3*s(ilo-1,jhi+1,klo+1) + &
                    s(ilo-1,jhi+1,klo+2))
            else
               s(ilo-1,jhi+1,klo-1) = &
                    (2*s(ilo-1,jhi+1,klo) - s(ilo-1,jhi+1,klo+1))
            end if
         end if

         if (nup .gt. 0) then
            if (khi-2 .ge. ks) then
               s(ilo-1,jhi+1,khi+1) = &
                    (3*s(ilo-1,jhi+1,khi) - 3*s(ilo-1,jhi+1,khi-1) + &
                    s(ilo-1,jhi+1,khi-2))
            else
               s(ilo-1,jhi+1,khi+1) = &
                    (2*s(ilo-1,jhi+1,khi) - s(ilo-1,jhi+1,khi-1))
            end if
         end if

      end if
!
! ****************************************************************************
!
      if ((nrgt .gt. 0) .and. (nbot .gt. 0)) then
         if (jlo+2 .le. je) then
            do k = slo(3), shi(3)
               s(ihi+1,jlo-1,k) = half * &
                    (3*s(ihi+1,jlo,k) - 3*s(ihi+1,jlo+1,k) +  &
                    s(ihi+1,jlo+2,k))
            end do
         else
            do k = slo(3), shi(3)
               s(ihi+1,jlo-1,k) = half * &
                    (2*s(ihi+1,jlo,k) - s(ihi+1,jlo+1,k))
            end do
         end if

         if (ihi-2 .ge. is) then
            do k = slo(3), shi(3)
               s(ihi+1,jlo-1,k) = s(ihi+1,jlo-1,k) + half * &
                    (3*s(ihi,jlo-1,k) - 3*s(ihi-1,jlo-1,k) +  &
                    s(ihi-2,jlo-1,k))
            end do
         else
            do k = slo(3), shi(3)
               s(ihi+1,jlo-1,k) = s(ihi+1,jlo-1,k) + half * &
                    (2*s(ihi,jlo-1,k) - s(ihi-1,jlo-1,k))
            end do
         end if

         if (ndwn .gt. 0) then
            if (klo+2 .le. ke) then
               s(ihi+1,jlo-1,klo-1) = &
                    (3*s(ihi+1,jlo-1,klo) - 3*s(ihi+1,jlo-1,klo+1) + &
                    s(ihi+1,jlo-1,klo+2))
            else
               s(ihi+1,jlo-1,klo-1) = &
                    (2*s(ihi+1,jlo-1,klo) - s(ihi+1,jlo-1,klo+1))
            end if
         end if

         if (nup .gt. 0) then
            if (khi-2 .ge. ks) then
               s(ihi+1,jlo-1,khi+1) = &
                    (3*s(ihi+1,jlo-1,khi) - 3*s(ihi+1,jlo-1,khi-1) + &
                    s(ihi+1,jlo-1,khi-2))
            else
               s(ihi+1,jlo-1,khi+1) = &
                    (2*s(ihi+1,jlo-1,khi) - s(ihi+1,jlo-1,khi-1))
            end if
         end if

      end if
!
! ****************************************************************************
!
      if ((nrgt .gt. 0) .and. (ntop .gt. 0)) then
         if (jhi-2 .ge. js) then
            do k = slo(3), shi(3)
               s(ihi+1,jhi+1,k) = half * &
                    (3*s(ihi+1,jhi,k) - 3*s(ihi+1,jhi-1,k) +  &
                    s(ihi+1,jhi-2,k))
            end do
         else
            do k = slo(3), shi(3)
               s(ihi+1,jhi+1,k) = half * &
                    (2*s(ihi+1,jhi,k) - s(ihi+1,jhi-1,k))
            end do
         end if

         if (ihi-2 .ge. is) then
            do k = slo(3), shi(3)
               s(ihi+1,jhi+1,k) = s(ihi+1,jhi+1,k) + half * &
                    (3*s(ihi,jhi+1,k) - 3*s(ihi-1,jhi+1,k) +  &
                    s(ihi-2,jhi+1,k))
            end do
         else
            do k = slo(3), shi(3)
               s(ihi+1,jhi+1,k) = s(ihi+1,jhi+1,k) + half * &
                    (2*s(ihi,jhi+1,k) - s(ihi-1,jhi+1,k))
            end do
         end if

         if (ndwn .gt. 0) then
            if (klo+2 .le. ke) then
               s(ihi+1,jhi+1,klo-1) = &
                    (3*s(ihi+1,jhi+1,klo) - 3*s(ihi+1,jhi+1,klo+1) + &
                    s(ihi+1,jhi+1,klo+2))
            else
               s(ihi+1,jhi+1,klo-1) = &
                    (2*s(ihi+1,jhi+1,klo) - s(ihi+1,jhi+1,klo+1))
            end if
         end if

         if (nup .gt. 0) then
            if (khi-2 .ge. ks) then
               s(ihi+1,jhi+1,khi+1) = &
                    (3*s(ihi+1,jhi+1,khi) - 3*s(ihi+1,jhi+1,khi-1) + &
                    s(ihi+1,jhi+1,khi-2))
            else
               s(ihi+1,jhi+1,khi+1) = &
                    (2*s(ihi+1,jhi+1,khi) - s(ihi+1,jhi+1,khi-1))
            end if
         end if

      end if
!
!    Next correct the i-k edges
!
      if ((nlft .gt. 0) .and. (ndwn .gt. 0)) then
         if (klo+2 .le. ke) then
            do j = slo(2), shi(2)
               s(ilo-1,j,klo-1) = half * &
                    (3*s(ilo-1,j,klo) - 3*s(ilo-1,j,klo+1) +  &
                    s(ilo-1,j,klo+2))
            end do
         else
            do j = slo(2), shi(2)
               s(ilo-1,j,klo-1) = half * &
                    (2*s(ilo-1,j,klo) - s(ilo-1,j,klo+1))
            end do
         end if

         if (ilo+2 .le. ie) then
            do j = slo(2), shi(2)
               s(ilo-1,j,klo-1) = s(ilo-1,j,klo-1) + half * &
                    (3*s(ilo,j,klo-1) - 3*s(ilo+1,j,klo-1) +  &
                    s(ilo+2,j,klo-1))
            end do
         else
            do j = slo(2), shi(2)
               s(ilo-1,j,klo-1) = s(ilo-1,j,klo-1) + half * &
                    (2*s(ilo,j,klo-1) - s(ilo+1,j,klo-1))
            end do
         end if
      end if
!
! ****************************************************************************
!
      if ((nlft .gt. 0) .and. (nup .gt. 0)) then
         if (khi-2 .ge. ks) then
            do j = slo(2), shi(2)
               s(ilo-1,j,khi+1) = half * &
                    (3*s(ilo-1,j,khi) - 3*s(ilo-1,j,khi-1) +  &
                    s(ilo-1,j,khi-2))
            end do
         else
            do j = slo(2), shi(2)
               s(ilo-1,j,khi+1) = half * &
                    (2*s(ilo-1,j,khi) - s(ilo-1,j,khi-1))
            end do
         end if

         if (ilo+2 .le. ie) then
            do j = slo(2), shi(2)
               s(ilo-1,j,khi+1) = s(ilo-1,j,khi+1) + half * &
                    (3*s(ilo,j,khi+1) - 3*s(ilo+1,j,khi+1) +  &
                    s(ilo+2,j,khi+1))
            end do
         else
            do j = slo(2), shi(2)
               s(ilo-1,j,khi+1) = s(ilo-1,j,khi+1) + half * &
                    (2*s(ilo,j,khi+1) - s(ilo+1,j,khi+1))
            end do
         end if
      end if
!
! ****************************************************************************
!
      if ((nrgt .gt. 0) .and. (ndwn .gt. 0)) then
         if (klo+2 .le. ke) then
            do j = slo(2), shi(2)
               s(ihi+1,j,klo-1) = half * &
                    (3*s(ihi+1,j,klo) - 3*s(ihi+1,j,klo+1) +  &
                    s(ihi+1,j,klo+2))
            end do
         else
            do j = slo(2), shi(2)
               s(ihi+1,j,klo-1) = half * &
                    (2*s(ihi+1,j,klo) - s(ihi+1,j,klo+1))
            end do
         end if

         if (ihi-2 .ge. is) then
            do j = slo(2), shi(2)
               s(ihi+1,j,klo-1) = s(ihi+1,j,klo-1) + half * &
                    (3*s(ihi,j,klo-1) - 3*s(ihi-1,j,klo-1) +  &
                    s(ihi-2,j,klo-1))
            end do
         else
            do j = slo(2), shi(2)
               s(ihi+1,j,klo-1) = s(ihi+1,j,klo-1) + half * &
                    (2*s(ihi,j,klo-1) - s(ihi-1,j,klo-1))
            end do
         end if
      end if
!
! ****************************************************************************
!
      if ((nrgt .gt. 0) .and. (nup .gt. 0)) then
         if (khi-2 .ge. ks) then
            do j = slo(2), shi(2)
               s(ihi+1,j,khi+1) = half * &
                    (3*s(ihi+1,j,khi) - 3*s(ihi+1,j,khi-1) +  &
                    s(ihi+1,j,khi-2))
            end do
         else
            do j = slo(2), shi(2)
               s(ihi+1,j,khi+1) = half * &
                    (2*s(ihi+1,j,khi) - s(ihi+1,j,khi-1))
            end do
         end if
         if (ihi-2 .ge. is) then
            do j = slo(2), shi(2)
               s(ihi+1,j,khi+1) = s(ihi+1,j,khi+1) + half * &
                    (3*s(ihi,j,khi+1) - 3*s(ihi-1,j,khi+1) +  &
                    s(ihi-2,j,khi+1))
            end do
         else
            do j = slo(2), shi(2)
               s(ihi+1,j,khi+1) = s(ihi+1,j,khi+1) + half * &
                    (2*s(ihi,j,khi+1) - s(ihi-1,j,khi+1))
            end do
         end if
      end if
!
!    Next correct the j-k edges
!
      if ((nbot .gt. 0) .and. (ndwn .gt. 0)) then
         if (klo+2 .le. ke) then
            do i = slo(1), shi(1)
               s(i,jlo-1,klo-1) = half * &
                    (3*s(i,jlo-1,klo) - 3*s(i,jlo-1,klo+1) +  &
                    s(i,jlo-1,klo+2))
            end do
         else
            do i = slo(1), shi(1)
               s(i,jlo-1,klo-1) = half * &
                    (2*s(i,jlo-1,klo) - s(i,jlo-1,klo+1))
            end do
         end if
         if (jlo+2 .le. je) then
            do i = slo(1), shi(1)
               s(i,jlo-1,klo-1) = s(i,jlo-1,klo-1) + half * &
                    (3*s(i,jlo,klo-1) - 3*s(i,jlo+1,klo-1) +  &
                    s(i,jlo+2,klo-1))
            end do
         else
            do i = slo(1), shi(1)
               s(i,jlo-1,klo-1) = s(i,jlo-1,klo-1) + half * &
                    (2*s(i,jlo,klo-1) - s(i,jlo+1,klo-1))
            end do
         end if
      end if
!
! ****************************************************************************
!
      if ((nbot .gt. 0) .and. (nup .gt. 0)) then
         if (khi-2 .ge. ks) then
            do i = slo(1), shi(1)
               s(i,jlo-1,khi+1) = half * &
                    (3*s(i,jlo-1,khi) - 3*s(i,jlo-1,khi-1) +  &
                    s(i,jlo-1,khi-2))
            end do
         else
            do i = slo(1), shi(1)
               s(i,jlo-1,khi+1) = half * &
                    (2*s(i,jlo-1,khi) - s(i,jlo-1,khi-1))
            end do
         end if

         if (jlo+2 .le. je) then
            do i = slo(1), shi(1)
               s(i,jlo-1,khi+1) = s(i,jlo-1,khi+1) + half * &
                    (3*s(i,jlo,khi+1) - 3*s(i,jlo+1,khi+1) +  &
                    s(i,jlo+2,khi+1))
            end do
         else
            do i = slo(1), shi(1)
               s(i,jlo-1,khi+1) = s(i,jlo-1,khi+1) + half * &
                    (2*s(i,jlo,khi+1) - s(i,jlo+1,khi+1))
            end do
         end if
      end if
!
! ****************************************************************************
!
      if ((ntop .gt. 0) .and. (ndwn .gt. 0)) then
         if (klo+2 .le. ke) then
            do i = slo(1), shi(1)
               s(i,jhi+1,klo-1) = half * &
                    (3*s(i,jhi+1,klo) - 3*s(i,jhi+1,klo+1) +  &
                    s(i,jhi+1,klo+2))
            end do
         else
            do i = slo(1), shi(1)
               s(i,jhi+1,klo-1) = half * &
                    (2*s(i,jhi+1,klo) - s(i,jhi+1,klo+1))
            end do
         end if
         if (jhi-2 .ge. js) then
            do i = slo(1), shi(1)
               s(i,jhi+1,klo-1) = s(i,jhi+1,klo-1) + half * &
                    (3*s(i,jhi,klo-1) - 3*s(i,jhi-1,klo-1) +  &
                    s(i,jhi-2,klo-1))
            end do
         else
            do i = slo(1), shi(1)
               s(i,jhi+1,klo-1) = s(i,jhi+1,klo-1) + half * &
                    (2*s(i,jhi,klo-1) - s(i,jhi-1,klo-1))
            end do
         end if
      end if
!
! ****************************************************************************
!
      if ((ntop .gt. 0) .and. (nup .gt. 0)) then
         if (khi-2 .ge. ks) then
            do i = slo(1), shi(1)
               s(i,jhi+1,khi+1) = half * &
                    (3*s(i,jhi+1,khi) - 3*s(i,jhi+1,khi-1) +  &
                    s(i,jhi+1,khi-2))
            end do
         else
            do i = slo(1), shi(1)
               s(i,jhi+1,khi+1) = half * &
                    (2*s(i,jhi+1,khi) - s(i,jhi+1,khi-1))
            end do
         end if
         if (jhi-2 .ge. js) then
            do i = slo(1), shi(1)
               s(i,jhi+1,khi+1) = s(i,jhi+1,khi+1) + half * &
                    (3*s(i,jhi,khi+1) - 3*s(i,jhi-1,khi+1) +  &
                    s(i,jhi-2,khi+1))
            end do
         else
            do i = slo(1), shi(1)
               s(i,jhi+1,khi+1) = s(i,jhi+1,khi+1) + half * &
                    (2*s(i,jhi,khi+1) - s(i,jhi-1,khi+1))
            end do
         end if
      end if

      end subroutine hoextraptocc
