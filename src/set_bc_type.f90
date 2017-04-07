module set_bc_type_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc_type                                             C
!                                                                      C
!  Author: J. Musser                                  Date: 05-FEB-17  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   subroutine set_bc_type(slo, shi, &
                          bc_ilo_type, bc_ihi_type, &
                          bc_jlo_type, bc_jhi_type, &
                          bc_klo_type, bc_khi_type, &
                          domlo, domhi, dx, dy, dz, &
                          xlength, ylength, zlength) &
               bind(c,name='set_bc_type')

      use bc, only: bc_defined, bc_type, bc_plane

      use bc, only: BC_X_e, BC_X_w
      use bc, only: BC_Y_n, BC_Y_s
      use bc, only: BC_Z_t, BC_Z_b


      use ic, only: NSW_, FSW_, PSW_
      use ic, only: PINF_, POUT_
      use ic, only: MINF_, MOUT_
      use ic, only: UNDEF_CELL, cycl_

      use geometry, only: cyclic_x, cyclic_y, cyclic_z
      use param, only: dimension_bc
      use toleranc, only: COMPARE
      use param1, only: zero, equal
      use calc_cell_module, only: calc_cell_bc_flow
      use calc_cell_module, only: calc_cell_bc_wall

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)

      real(c_real), intent(in) :: dx, dy, dz
      real(c_real)  , intent(in) :: xlength, ylength, zlength

      integer(c_int), intent(inout) :: bc_ilo_type&
         (domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(inout) :: bc_ihi_type&
         (domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(inout) :: bc_jlo_type&
         (domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(inout) :: bc_jhi_type&
         (domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(inout) :: bc_klo_type&
         (domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2)
      integer(c_int), intent(inout) :: bc_khi_type&
         (domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2)

      ! Local index for boundary condition
      integer :: type, bcv
      integer    nlft, nrgt, nbot, ntop, nup, ndwn
      integer :: i,j,k

      integer :: i_w, j_s, k_b, i_e, j_n, k_t

      nlft = 2!max(0,domlo(1)-slo(1))
      nbot = 2!max(0,domlo(2)-slo(2))
      ndwn = 2!max(0,domlo(3)-slo(3))

      nrgt = 2!max(0,shi(1)-domhi(1))
      ntop = 2!max(0,shi(2)-domhi(2))
      nup  = 2!max(0,shi(3)-domhi(3))


      if(cyclic_x) then
         bc_ilo_type(domlo(2):domhi(2),domlo(3):domhi(3),1) = cycl_
         bc_ihi_type(domlo(2):domhi(2),domlo(3):domhi(3),1) = cycl_
      else
         bc_ilo_type(:,:,1) = undef_cell
         bc_ihi_type(:,:,1) = undef_cell
      endif

      if(cyclic_y)then
         bc_jlo_type(:,:,1) = cycl_
         bc_jhi_type(:,:,1) = cycl_
      else
         bc_jlo_type(:,:,1) = undef_cell
         bc_jhi_type(:,:,1) = undef_cell
      endif

      if(cyclic_z) then
         bc_klo_type(:,:,1) = cycl_
         bc_khi_type(:,:,1) = cycl_
      else
         bc_klo_type(:,:,1) = undef_cell
         bc_khi_type(:,:,1) = undef_cell
      endif

      do bcv = 1, dimension_bc
         if (bc_defined(bcv)) then

            select case (trim(bc_type(bcv)))
               case('FREE_SLIP_WALL','FSW'); type = fsw_
               case('NO_SLIP_WALL'  ,'NSW'); type = nsw_
               case('PAR_SLIP_WALL' ,'PSW'); type = psw_
               case('P_INFLOW'      ,'PI' ); type = pinf_
               case('P_OUTFLOW'     ,'PO' ); type = pout_
               case('MASS_INFLOW'   ,'MI' ); type = minf_
               case('MASS_OUTFLOW'  ,'MO' ); type = mout_
               case default
                  write(6,*) 'unknown bc type'
                  stop 7655
            end select

            select case(type)
            case(nsw_, fsw_, psw_)
               call calc_cell_bc_wall(bcv, domlo, domhi, xlength, ylength, &
                  zlength, dx, dy, dz, i_w, i_e, j_s, j_n, k_b, k_t)
            case(pinf_, pout_, minf_, mout_)
               call calc_cell_bc_flow(bcv, xlength, ylength, zlength, &
                  dx, dy, dz, i_w, i_e, j_s, j_n, k_b, k_t)
            end select

            if (i_w == i_e) then
               if(i_w == domlo(1)-1) then
                  bc_ilo_type(j_s:j_n, k_b:k_t,1) = type
                  bc_ilo_type(j_s:j_n, k_b:k_t,2) = bcv

                  bc_plane(bcv) = 'E' ! Fluid is to East of bc!

               else if(i_w == domhi(1)+1) then
                  bc_ihi_type(j_s:j_n, k_b:k_t,1) = type
                  bc_ihi_type(j_s:j_n, k_b:k_t,2) = bcv

                  bc_plane(bcv) = 'W' ! Fluid is to West of bc!

               endif
            endif

            if (j_s == j_n) then
               if(j_s == domlo(2)-1) then
                  bc_jlo_type(i_w:i_e, k_b:k_t,1) = type
                  bc_jlo_type(i_w:i_e, k_b:k_t,2) = bcv

                  bc_plane(bcv) = 'N' ! Fluid is to North of bc!

               else if(j_s == domhi(2)+1) then
                  bc_jhi_type(i_w:i_e, k_b:k_t,1) = type
                  bc_jhi_type(i_w:i_e, k_b:k_t,2) = bcv

                  bc_plane(bcv) = 'S' ! Fluid is to South of bc!

               endif
            endif

            if (k_b == k_t) then
               if(k_b == domlo(3)-1) then
                  bc_klo_type(i_w:i_e, j_s:j_n,1) = type
                  bc_klo_type(i_w:i_e, j_s:j_n,2) = bcv

                  bc_plane(bcv) = 'T' ! Fluid is to Top of bc!

               elseif(k_b == domhi(3)+1) then
                  bc_khi_type(i_w:i_e, j_s:j_n,1) = type
                  bc_khi_type(i_w:i_e, j_s:j_n,2) = bcv

                  bc_plane(bcv) = 'B' ! Fluid is to Bottom of bc!

               endif
            endif

         endif
      enddo


! Edge cases
! --------------------------------------------------------------------------------------------

      if (cyclic_x)then

         do j=1,nbot
            bc_ilo_type(domlo(2)-j,domlo(3)-1:domhi(3)+1,:) = bc_jlo_type(domlo(1),domlo(3)-1:domhi(3)+1,:)
            bc_ihi_type(domlo(2)-j,domlo(3)-1:domhi(3)+1,:) = bc_jlo_type(domhi(1),domlo(3)-1:domhi(3)+1,:)
         enddo
         do j=1,ntop
            bc_ilo_type(domhi(2)+j,domlo(3)-1:domhi(3)+1,:) = bc_jhi_type(domlo(1),domlo(3)-1:domhi(3)+1,:)
            bc_ihi_type(domhi(2)+j,domlo(3)-1:domhi(3)+1,:) = bc_jhi_type(domhi(1),domlo(3)-1:domhi(3)+1,:)
         enddo
         do k=1,ndwn
            bc_ilo_type(domlo(2)-1:domhi(2)+1,domlo(3)-k,:) = bc_klo_type(domlo(1),domlo(2)-1:domhi(2)+1,:)
            bc_ihi_type(domlo(2)-1:domhi(2)+1,domlo(3)-k,:) = bc_klo_type(domhi(1),domlo(2)-1:domhi(2)+1,:)
         enddo
         do k=1,nup
            bc_ilo_type(domlo(2)-1:domhi(2)+1,domhi(3)+k,:) = bc_khi_type(domlo(1),domlo(2)-1:domhi(2)+1,:)
            bc_ihi_type(domlo(2)-1:domhi(2)+1,domhi(3)+k,:) = bc_khi_type(domhi(1),domlo(2)-1:domhi(2)+1,:)
         enddo

      else

         do j=1,nbot
            bc_ilo_type(domlo(2)-j,domlo(3)-1:domhi(3)+1,:) = bc_ilo_type(domlo(2),domlo(3)-1:domhi(3)+1,:)
            bc_ihi_type(domlo(2)-j,domlo(3)-1:domhi(3)+1,:) = bc_ihi_type(domlo(2),domlo(3)-1:domhi(3)+1,:)
         enddo
         do j=1,ntop
            bc_ilo_type(domhi(2)+j,domlo(3)-1:domhi(3)+1,:) = bc_ilo_type(domhi(2),domlo(3)-1:domhi(3)+1,:)
            bc_ihi_type(domhi(2)+j,domlo(3)-1:domhi(3)+1,:) = bc_ihi_type(domhi(2),domlo(3)-1:domhi(3)+1,:)
         enddo
         do k=1,ndwn
            bc_ilo_type(domlo(2)-1:domhi(2)+1,domlo(3)-k,:) = bc_ilo_type(domlo(2)-1:domhi(2)+1,domlo(3),:)
            bc_ihi_type(domlo(2)-1:domhi(2)+1,domlo(3)-k,:) = bc_ihi_type(domlo(2)-1:domhi(2)+1,domlo(3),:)
         enddo
         do k=1,nup
            bc_ilo_type(domlo(2)-1:domhi(2)+1,domhi(3)+k,:) = bc_ilo_type(domlo(2)-1:domhi(2)+1,domhi(3),:)
            bc_ihi_type(domlo(2)-1:domhi(2)+1,domhi(3)+k,:) = bc_ihi_type(domlo(2)-1:domhi(2)+1,domhi(3),:)
         enddo
      endif

      if(cyclic_y)then
         do i=1,nlft
            bc_jlo_type(domlo(1)-i,domlo(3)-1:domhi(3)+1,:) = bc_ilo_type(domlo(2),domlo(3)-1:domhi(3)+1,:)
            bc_jhi_type(domlo(1)-i,domlo(3)-1:domhi(3)+1,:) = bc_ilo_type(domhi(2),domlo(3)-1:domhi(3)+1,:)
         enddo
         do i=1,nrgt
            bc_jlo_type(domhi(1)+i,domlo(3)-1:domhi(3)+1,:) = bc_ihi_type(domlo(2),domlo(3)-1:domhi(3)+1,:)
            bc_jhi_type(domhi(1)+i,domlo(3)-1:domhi(3)+1,:) = bc_ihi_type(domhi(1),domlo(3)-1:domhi(3)+1,:)
         enddo
         do k=1,ndwn
            bc_jlo_type(domlo(1)-1:domhi(1)+1,domlo(3)-k,:) = bc_klo_type(domlo(1)-1:domhi(1)+1,domlo(2),:)
            bc_jhi_type(domlo(1)-1:domhi(1)+1,domlo(3)-k,:) = bc_klo_type(domlo(1)-1:domhi(1)+1,domhi(2),:)
         enddo
         do k=1,nup
            bc_jlo_type(domlo(1)-1:domhi(1)+1,domhi(3)+k,:) = bc_khi_type(domlo(1)-1:domhi(1)+1,domlo(2),:)
            bc_jhi_type(domlo(1)-1:domhi(1)+1,domhi(3)+k,:) = bc_khi_type(domlo(1)-1:domhi(1)+1,domhi(2),:)
         enddo

      else

         do i=1,nlft
            bc_jlo_type(domlo(1)-i,domlo(3)-1:domhi(3)+1,:) = bc_jlo_type(domlo(1),domlo(3)-1:domhi(3)+1,:)
            bc_jhi_type(domlo(1)-i,domlo(3)-1:domhi(3)+1,:) = bc_jhi_type(domlo(1),domlo(3)-1:domhi(3)+1,:)
         enddo
         do i=1,nrgt
            bc_jlo_type(domhi(1)+i,domlo(3)-1:domhi(3)+1,:) = bc_jlo_type(domhi(1),domlo(3)-1:domhi(3)+1,:)
            bc_jhi_type(domhi(1)+i,domlo(3)-1:domhi(3)+1,:) = bc_jhi_type(domhi(1),domlo(3)-1:domhi(3)+1,:)
         enddo
         do k=1,ndwn
            bc_jlo_type(domlo(1)-1:domhi(1)+1,domlo(3)-k,:) = bc_jlo_type(domlo(1)-1:domhi(1)+1,domlo(3),:)
            bc_jhi_type(domlo(1)-1:domhi(1)+1,domlo(3)-k,:) = bc_jhi_type(domlo(1)-1:domhi(1)+1,domlo(3),:)
         enddo
         do k=1,nup
            bc_jlo_type(domlo(1)-1:domhi(1)+1,domhi(3)+k,:) = bc_jlo_type(domlo(1)-1:domhi(1)+1,domhi(3),:)
            bc_jhi_type(domlo(1)-1:domhi(1)+1,domhi(3)+k,:) = bc_jhi_type(domlo(1)-1:domhi(1)+1,domhi(3),:)
         enddo
      endif


      if(cyclic_z) then
         do i=1,nlft
            bc_klo_type(domlo(1)-i,domlo(2)-1:domhi(2)+1,:) = bc_ilo_type(domlo(2)-1:domhi(2)+1,domlo(3),:)
            bc_khi_type(domlo(1)-i,domlo(2)-1:domhi(2)+1,:) = bc_ilo_type(domlo(2)-1:domhi(2)+1,domhi(3),:)
         enddo
         do i=1,nrgt
            bc_klo_type(domhi(1)+i,domlo(2)-1:domhi(2)+1,:) = bc_ihi_type(domlo(2)-1:domhi(2)+1,domlo(3),:)
            bc_khi_type(domhi(1)+i,domlo(2)-1:domhi(2)+1,:) = bc_ihi_type(domlo(2)-1:domhi(2)+1,domhi(3),:)
         enddo
         do j=1,nbot
            bc_klo_type(domlo(1)-1:domhi(1)+1,domlo(2)-j,:) = bc_jlo_type(domlo(1)-1:domhi(1)+1,domlo(3),:)
            bc_khi_type(domlo(1)-1:domhi(1)+1,domlo(2)-j,:) = bc_jlo_type(domlo(1)-1:domhi(1)+1,domhi(3),:)
         enddo
         do j=1,ntop
            bc_klo_type(domlo(1)-1:domhi(1)+1,domhi(2)+j,:) = bc_jhi_type(domlo(1)-1:domhi(1)+1,domlo(3),:)
            bc_khi_type(domlo(1)-1:domhi(1)+1,domhi(2)+j,:) = bc_jhi_type(domlo(1)-1:domhi(1)+1,domhi(3),:)
         enddo

      else
         do i=1,nlft
            bc_klo_type(domlo(1)-i,domlo(2)-1:domhi(2)+1,:) = bc_klo_type(domlo(1),domlo(2)-1:domhi(2)+1,:)
            bc_khi_type(domlo(1)-i,domlo(2)-1:domhi(2)+1,:) = bc_khi_type(domlo(1),domlo(2)-1:domhi(2)+1,:)
         enddo
         do i=1,nrgt
            bc_klo_type(domhi(1)+i,domlo(2)-1:domhi(2)+1,:) = bc_klo_type(domhi(1),domlo(2)-1:domhi(2)+1,:)
            bc_khi_type(domhi(1)+i,domlo(2)-1:domhi(2)+1,:) = bc_khi_type(domhi(1),domlo(2)-1:domhi(2)+1,:)
         enddo
         do j=1,nbot
            bc_klo_type(domlo(1)-1:domhi(1)+1,domlo(2)-j,:) = bc_klo_type(domlo(1)-1:domhi(1)+1,domlo(2),:)
            bc_khi_type(domlo(1)-1:domhi(1)+1,domlo(2)-j,:) = bc_khi_type(domlo(1)-1:domhi(1)+1,domlo(2),:)
         enddo
         do j=1,ntop
            bc_klo_type(domlo(1)-1:domhi(1)+1,domhi(2)+j,:) = bc_klo_type(domlo(1)-1:domhi(1)+1,domhi(2),:)
            bc_khi_type(domlo(1)-1:domhi(1)+1,domhi(2)+j,:) = bc_khi_type(domlo(1)-1:domhi(1)+1,domhi(2),:)
         enddo
      endif


      if (0.eq.1) then
      i = slo(1)
      write(6,"(2/,'bc_ilo_type for i =',i3)") i
      do k=shi(3),slo(3),-1
         do j=slo(2),shi(2)-1
            write(6,"(i4)",advance='no') bc_ilo_type(j,k,1)
         enddo
         j=shi(2)
         write(6,"(i4)",advance='yes') bc_ilo_type(j,k,1)
      enddo

      i = shi(1)
      write(6,"(2/,'bc_ihi_type for i =',i3)") i
      do k=shi(3),slo(3),-1
         do j=slo(2),shi(2)-1
            write(6,"(i4)",advance='no') bc_ihi_type(j,k,1)
         enddo
         j=shi(2)
         write(6,"(i4)",advance='yes') bc_ihi_type(j,k,1)
      enddo

      j = slo(2)
      write(6,"(2/,'bc_jlo_type for j =',i3)") j
      do k=slo(3),shi(3)
         do i=slo(1),shi(1)-1
            write(6,"(i4)",advance='no') bc_jlo_type(i,k,1)
         enddo
         i=shi(1)
         write(6,"(i4)",advance='yes') bc_jlo_type(i,k,1)
      enddo

      j = shi(2)
      write(6,"(2/,'bc_jhi_type for j =',i3)") j
      do k=slo(3),shi(3)
         do i=slo(1),shi(1)-1
            write(6,"(i4)",advance='no') bc_jhi_type(i,k,1)
         enddo
         i=shi(1)
         write(6,"(i4)",advance='yes') bc_jhi_type(i,k,1)
      enddo

      k = slo(3)
      write(6,"(2/,'bc_klo_type for k =',i3)") j
      do j=slo(2),shi(2)
         do i=slo(1),shi(1)-1
            write(6,"(i4)",advance='no') bc_klo_type(i,j,1)
         end do
         i=shi(1)
         write(6,"(i4)",advance='yes') bc_klo_type(i,j,1)
      enddo


      k = shi(3)
      write(6,"(2/,'bc_khi_type for k =',i3)") j
      do j=slo(2),shi(2)
         do i=slo(1),shi(1)-1
            write(6,"(i4)",advance='no') bc_khi_type(i,j,1)
         enddo
         i=shi(1)
         write(6,"(i4)",advance='yes') bc_khi_type(i,j,1)
      enddo
      end if

   end subroutine set_bc_type

   end module set_bc_type_module
