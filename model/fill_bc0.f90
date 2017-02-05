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
subroutine fill_bc0(s, slo, shi, flag, vtype)&
   bind(C, name="fill_bc0")


! Global modules
!--------------------------------------------------------------------//
      use geometry      , only: domlo, domhi
      use set_bc0_module, only: set_bc_type
      use iso_c_binding , only: c_int
      use bl_fort_module, only: c_real

      implicit none


! Dummy arguments
!``````````````````````````````````````````````````````````````````````
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: vtype

      real(c_real), intent(inout) ::  s&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)


! Local variables
!``````````````````````````````````````````````````````````````````````
      integer    nlft, nrgt, nbot, ntop, nup, ndwn
      integer    ilo, ihi, jlo, jhi, klo, khi
      integer    i, j, k

      integer :: bc_i_type(2,slo(2):shi(2),slo(3):shi(3))
      integer :: bc_j_type(2,slo(1):shi(1),slo(3):shi(3))
      integer :: bc_k_type(2,slo(1):shi(1),slo(2):shi(2))

      integer :: bc_i_ptr(2,slo(2):shi(2),slo(3):shi(3))
      integer :: bc_j_ptr(2,slo(1):shi(1),slo(3):shi(3))
      integer :: bc_k_ptr(2,slo(1):shi(1),slo(2):shi(2))
!......................................................................

      if(vtype /= 0) return

      nlft = max(0,domlo(1)-slo(1))
      nbot = max(0,domlo(2)-slo(2))
      ndwn = max(0,domlo(3)-slo(3))

      nrgt = max(0,shi(1)-domhi(1))
      ntop = max(0,shi(2)-domhi(2))
      nup  = max(0,shi(3)-domhi(3))

      call set_bc_type(slo, shi, bc_i_type, bc_j_type, bc_k_type, &
         bc_i_ptr, bc_j_ptr, bc_k_ptr, flag)

      if (nlft .gt. 0) then
         ilo = domlo(1)
         do i = 1, nlft
            do k=slo(3),shi(3)
               do j=slo(2),shi(2)
                  if(bc_i_type(1,j,k) >= 100) then
                     s(ilo-i,j,k) = s(ilo,j,k)
                  endif
               end do
            end do
         end do
      endif

      if (nrgt .gt. 0) then
         ihi = domhi(1)
         do i = 1, nrgt
            do k=slo(3),shi(3)
               do j=slo(2),shi(2)
                  if(bc_i_type(2,j,k) >= 100) then
                     s(ihi+i,j,k) = s(ihi,j,k)
                  endif
               end do
            end do
         end do
      endif

      if (nbot .gt. 0) then
         jlo = domlo(2)
         do j = 1, nbot
            do k=slo(3),shi(3)
               do i=slo(1),shi(1)
                  if(bc_j_type(1,i,k) >= 100) then
                     s(i,jlo-j,k) = s(i,jlo,k)
                  endif
               end do
            end do
         end do
      endif

      if (ntop .gt. 0) then
         jhi = domhi(2)
         do j = 1, ntop
            do k=slo(3),shi(3)
               do i=slo(1),shi(1)
                  if(bc_j_type(2,i,k) >= 100) then
                     s(i,jhi+j,k) = s(i,jhi,k)
                  endif
               end do
            end do
         end do
      endif

      if (ndwn .gt. 0) then
         klo = domlo(3)
         do k = 1, ndwn
            do j=slo(2),shi(2)
               do i=slo(1),shi(1)
                  if(bc_k_type(1,i,j) >= 100) then
                     s(i,j,klo-k) = s(i,j,klo)
                  endif
               end do
            end do
         end do
      endif

      if (nup .gt. 0) then
         khi = domhi(3)
         do k = 1, nup
            do j=slo(2),shi(2)
               do i=slo(1),shi(1)
                  if(bc_k_type(1,i,j) >= 100) then
                     s(i,j,khi+k) = s(i,j,khi)
                  endif
               end do
            end do
         end do
      endif

! fill edges --------------------------------------------

      if (nlft .gt. 0) then
         ilo = domlo(1)
         if (nbot .gt. 0) then
            jlo = domlo(2)
            do i = 1, nlft
               do j = 1, nbot
                  do k=slo(3),shi(3)
                     s(ilo-i,jlo-j,k) = s(ilo,jlo,k)
                  end do
               end do
            end do
         endif

         if (ntop .gt. 0) then
            jhi = domhi(2)
            do i = 1, nlft
               do j = 1, ntop
                  do k=slo(3),shi(3)
                     s(ilo-i,jhi+j,k) = s(ilo,jhi,k)
                  end do
               end do
            end do
         endif

         if (ndwn .gt. 0) then
            klo = domlo(3)
            do i = 1, nlft
               do k = 1, ndwn
                  do j=slo(2),shi(2)
                     s(ilo-i,j,klo-k) = s(ilo-i,j,klo)
                  end do
               end do
            end do
         endif

         if (nup .gt. 0) then
            khi = domhi(3)
            do i = 1, nlft
               do k = 1, nup
                  do j=slo(2),shi(2)
                     s(ilo-i,j,khi+k) = s(ilo-i,j,khi)
                  end do
               end do
            end do
         endif
      endif

      if (nrgt .gt. 0) then
         ihi = domhi(1)
         if (nbot .gt. 0) then
            jlo = domlo(2)
            do i = 1, nrgt
               do j = 1, nbot
                  do k=slo(3),shi(3)
                     s(ihi+i,jlo-j,k) = s(ihi,jlo,k)
                  end do
               end do
            end do
         endif

         if (ntop .gt. 0) then
            jhi = domhi(2)
            do i = 1, nrgt
               do j = 1, ntop
                  do k=slo(3),shi(3)
                     s(ihi+i,jhi+j,k) = s(ihi,jhi,k)
                  end do
               end do
            end do
         endif

         if (ndwn .gt. 0) then
            klo = domlo(3)
            do i = 1, nrgt
               do k = 1, ndwn
                  do j=slo(2),shi(2)
                     s(ihi+i,j,klo-k) = s(ihi,j,klo)
                  end do
               end do
            end do
         endif

         if (nup .gt. 0) then
            khi = domhi(3)
            do i = 1, nrgt
               do k = 1, nup
                  do j=slo(2),shi(2)
                     s(ihi+i,j,khi+k) = s(ihi,j,khi)
                  end do
               end do
            end do
         endif
      endif

      return
   end subroutine fill_bc0
