!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: set_ls_inflow                                           !
!  Purpose: This subroutine does the initial setting of all boundary   !
!  conditions. The user specifications of the boundary conditions are  !
!  checked for veracity in various check_data/ routines:               !
!  (e.g., check_boundary_conditions).                                  !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine set_ls_inflow(phi_ls, slo, shi, &
                            bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
                            bc_klo_type, bc_khi_type, domlo, domhi, ng, nref, dx) &
      bind(C, name="set_ls_inflow")

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

      use bc, only: minf_

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)
      integer(c_int), intent(in   ) :: ng, nref
      real(rt)      , intent(inout) :: dx(3)

      ! Recall that phi_ls is nodal
      real(rt), intent(inout) :: phi_ls&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer(c_int), intent(in   ) :: &
           bc_ilo_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_ihi_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_jlo_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_jhi_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_klo_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
           bc_khi_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

! Local variables
!--------------------------------------------------------------------//
! local index for boundary condition
      integer :: i,j,k

      integer    nlft, nrgt, nbot, ntop, nup, ndwn

      real(rt) :: dx_fine(3)
      real(rt) :: offset
!--------------------------------------------------------------------//

      offset = 1.d-8

      ! Here if the level set (slo,shi) is at a finer resolution (by nref) than
      !  the boundary condition routines,
      !  the domain boundaries (domlo,domhi), and dx,
      !     then we make sure to adjust by nref
      dx_fine(:) = dx(:) / dble(nref)

      nlft = max(0,nref*domlo(1)-slo(1))
      nbot = max(0,nref*domlo(2)-slo(2))
      ndwn = max(0,nref*domlo(3)-slo(3))

      nrgt = max(0,shi(1)-(nref*domhi(1)+1))
      ntop = max(0,shi(2)-(nref*domhi(2)+1))
      nup  = max(0,shi(3)-(nref*domhi(3)+1))

      if (nlft .gt. 0) then
         do k=slo(3),shi(3)
            do j=slo(2),shi(2)
               if (any(bc_ilo_type(j/nref:j/nref+1,k/nref:k/nref+1,1) == MINF_)) then

                  do i=slo(1),-1
                      ! Values outside the domain that are positive can be ignored
                      if (phi_ls(i,j,k) .gt. 0.) then
                          phi_ls(i,j,k) = dble(i)*dx_fine(1)-offset

                      ! Values outside the domain that are negative should see the min
                      else
                          phi_ls(i,j,k) = min(phi_ls(i,j,k),dble(i)*dx_fine(1)-offset)
                      end if
                  end do

                  if (phi_ls(0,j,k) .gt. 0.) &
                      phi_ls(0,j,k) = -offset

                  do i=1,shi(1)  ! Only min values that are already interior
                      if (phi_ls(i,j,k) .gt. 0.) &
                         phi_ls(i,j,k) = min(phi_ls(i,j,k),dble(i)*dx_fine(1)-offset)
                  end do

               end if
            end do
         end do
      endif

      if (nrgt .gt. 0) then
         do k=slo(3),shi(3)
            do j=slo(2),shi(2)
               if (any(bc_ihi_type(j/nref:j/nref+1,k/nref:k/nref+1,1) == MINF_)) then

                  do i=slo(1),(domhi(1)+1)*nref-1  ! Only min values that are already interior
                      if (phi_ls(i,j,k) .gt. 0.) &
                          phi_ls(i,j,k) = min(phi_ls(i,j,k),dble(domhi(1)-i+1)*dx_fine(1)-offset)
                  end do

                  if (phi_ls( (domhi(1)+1)*nref,j,k) .gt. 0.) &
                      phi_ls( (domhi(1)+1)*nref,j,k) = -offset

                  do i=(domhi(1)+1)*nref+1,shi(1)
                      ! Values outside the domain that are positive can be ignored
                      if (phi_ls(i,j,k) .gt. 0.) then
                          phi_ls(i,j,k) = dble(domhi(1)-i+1)*dx_fine(1)-offset

                      ! Values outside the domain that are negative should see the min
                      else
                          phi_ls(i,j,k) = min(phi_ls(i,j,k),dble(domhi(1)-i+1)*dx_fine(1)-offset)
                      end if
                  end do
               end if
            end do
         end do
      endif

      if (nbot .gt. 0) then
         do k=slo(3),shi(3)
            do i=slo(1),shi(1)
               if (any(bc_jlo_type(i/nref:i/nref+1,k/nref:k/nref+1,1) == MINF_)) then

                  do j=slo(2),-1
                      ! Values outside the domain that are positive can be ignored
                      if (phi_ls(i,0,k) .gt. 0.) then
                          phi_ls(i,j,k) = dble(j)*dx_fine(2)-offset

                      ! Values outside the domain that are negative should see the min
                      else
                          phi_ls(i,j,k) = min(phi_ls(i,j,k),dble(j)*dx_fine(2)-offset)
                      end if
                  end do

                  if (phi_ls(i,0,k) .gt. 0.) &
                      phi_ls(i,0,k) = -offset

                  do j=1,shi(2)   ! Only min values that are already interior
                      if (phi_ls(i,j,k) .gt. 0.) &
                         phi_ls(i,j,k) = min(phi_ls(i,j,k),dble(j)*dx_fine(2)-offset)
                  end do
               end if
            end do
         end do
      endif

      if (ntop .gt. 0) then
         do k = slo(3),shi(3)
            do i = slo(1),shi(1)
               if (any(bc_jhi_type(i/nref:i/nref+1,k/nref:k/nref+1,1) == MINF_)) then

                  do j=slo(2),(domhi(2)+1)*nref-1  ! Only min values that are already interior
                      if (phi_ls(i,j,k) .gt. 0.) &
                          phi_ls(i,j,k) = min(phi_ls(i,j,k),dble(domhi(2)-j+1)*dx_fine(2)-offset)
                  end do

                  if (phi_ls(i, (domhi(2)+1)*nref,k) .gt. 0.) &
                      phi_ls(i, (domhi(2)+1)*nref,k) = -offset

                  do j=(domhi(2)+1)*nref+1,shi(2)
                      ! Values outside the domain that are positive can be ignored
                      if (phi_ls(i,j,k) .gt. 0.) then
                          phi_ls(i,j,k) = dble(domhi(2)-j+1)*dx_fine(2)-offset

                      ! Values outside the domain that are negative should see the min
                      else
                          phi_ls(i,j,k) = min(phi_ls(i,j,k),dble(domhi(2)-j+1)*dx_fine(2)-offset)
                      end if
                  end do
               end if
            end do
         end do
      endif

      if (ndwn .gt. 0) then
         do j=slo(2),shi(2)
            do i=slo(1),shi(1)
               if (any(bc_klo_type(i/nref:i/nref+1,j/nref:j/nref+1,1) == MINF_)) then

                  do k=slo(3),-1
                      ! Values outside the domain that are positive can be ignored
                      if (phi_ls(i,j,0) .gt. 0.) then
                          phi_ls(i,j,k) = dble(k)*dx_fine(3)-offset

                      ! Values outside the domain that are negative should see the min
                      else
                          phi_ls(i,j,k) = min(phi_ls(i,j,k),dble(k)*dx_fine(3)-offset)
                      end if
                  end do

                  if (phi_ls(i,j,0) .gt. 0.) &
                      phi_ls(i,j,0) = -offset

                  do k=1,shi(3)  ! Only min values that are already interior
                      if (phi_ls(i,j,k) .gt. 0.) &
                          phi_ls(i,j,k) = min(phi_ls(i,j,k),dble(k)*dx_fine(3)-offset)
                  end do
               end if
            end do
         end do
      endif

      if (nup .gt. 0) then
         do j=slo(2),shi(2)
            do i=slo(1),shi(1)
               if (any(bc_khi_type(i/nref:i/nref+1,j/nref:j/nref+1,1) == MINF_)) then

                  do k=slo(3),(domhi(3)+1)*nref-1  ! Only min values that are already interior
                      if (phi_ls(i,j,k) .gt. 0.) &
                          phi_ls(i,j,k) = min(phi_ls(i,j,k),dble(domhi(3)-k+1)*dx_fine(3)-offset)
                  end do

                  if (phi_ls(i,j,(domhi(3)+1)*nref) .gt. 0.) &
                      phi_ls(i,j,(domhi(3)+1)*nref) = -offset

                  do k=(domhi(3)+1)*nref+1,shi(3)
                      ! Values outside the domain that are positive can be ignored
                      if (phi_ls(i,j,k) .gt. 0.) then
                          phi_ls(i,j,k) = dble(domhi(3)-k+1)*dx_fine(3)-offset

                      ! Values outside the domain that are negative should see the min
                      else
                          phi_ls(i,j,k) = min(phi_ls(i,j,k),dble(domhi(3)-k+1)*dx_fine(3)-offset)
                      end if
                  end do

               end if
            end do
         end do
      endif

   end subroutine set_ls_inflow
