!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_drag_bcs                                            C
!                                                                      C
!  Purpose: This subroutine does the initial setting of all boundary   C
!  conditions. The user specifications of the boundary conditions are  C
!  checked for veracity in various check_data/ routines:               C
!  (e.g., check_boundary_conditions).                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
subroutine set_drag_bcs( &
                        ulo, uhi, vlo, vhi, wlo, whi, &
                        drag_u, drag_v, drag_w, &
                        bc_ilo_type, bc_ihi_type, &
                        bc_jlo_type, bc_jhi_type, &
                        bc_klo_type, bc_khi_type, domlo, domhi, ng ) &
     bind(C, name="set_drag_bcs")

     use amrex_fort_module, only : c_real => amrex_real
     use iso_c_binding , only: c_int

     use bc, only: nsw_, fsw_, psw_, pinf_, pout_, minf_

     implicit none

     integer(c_int), intent(in   ) :: ulo(3),uhi(3)
     integer(c_int), intent(in   ) :: vlo(3),vhi(3)
     integer(c_int), intent(in   ) :: wlo(3),whi(3)
     integer(c_int), intent(in   ) :: domlo(3),domhi(3), ng

     real(c_real), intent(inout) ::  drag_u&
        (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
     real(c_real), intent(inout) ::  drag_v&
        (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
     real(c_real), intent(inout) ::  drag_w&
        (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

     integer(c_int), intent(in   ) :: &
       bc_ilo_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
       bc_ihi_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
       bc_jlo_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
       bc_jhi_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
       bc_klo_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
       bc_khi_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

! Local variables
!--------------------------------------------------------------------//
  integer(c_int) :: i,j,k
  integer(c_int) :: ilo, ihi, jlo, jhi, klo, khi
  integer(c_int) :: ibc, jbc, kbc

!--------------------------------------------------------------------//
!-  Do tangential corrections first -- we will want to override with normal corrections
!--------------------------------------------------------------------//

  if (ulo(1).le.domlo(1)) then
     do k = vlo(3), vhi(3)
      do j = vlo(2), vhi(2)
        jbc = min(j,domhi(2)+ng)
        if (bc_ilo_type(jbc,k,1) == FSW_) then
            drag_v(domlo(1),j,k) = drag_v(domlo(1),j,k) + drag_v(domlo(1)-1,j,k)
        end if
        drag_v(domlo(1)-1,j,k) = 1.d200
     end do
     end do
     do k = wlo(3), whi(3)
      kbc = min(k,domhi(3)+ng)
      do j = wlo(2), whi(2)
        if (bc_ilo_type(j,kbc,1) == FSW_) then
            drag_w(domlo(1),j,k) = drag_w(domlo(1),j,k) + drag_w(domlo(1)-1,j,k)
        end if
        drag_w(domlo(1)-1,j,k) = 1.d200
     end do
     end do
  endif

  if (uhi(1).ge.domhi(1)+1) then
     do k = vlo(3), vhi(3)
      do j = vlo(2), vhi(2)
        jbc = min(j,domhi(2)+ng)
        if (bc_ihi_type(jbc,k,1) == FSW_) then
            drag_v(domhi(1),j,k) = drag_v(domhi(1),j,k) + drag_v(domhi(1)+1,j,k)
        end if
        drag_v(domhi(1)+1,j,k) = 1.d200
     end do
     end do
     do k = wlo(3), whi(3)
      kbc = min(k,domhi(3)+ng)
      do j = wlo(2), whi(2)
        if (bc_ihi_type(j,kbc,1) == FSW_) then
            drag_w(domhi(1),j,k) = drag_w(domhi(1),j,k) + drag_w(domhi(1)+1,j,k)
        end if
        drag_w(domhi(1)+1,j,k) = 1.d200
     end do
     end do
  endif

  if (vlo(2).le.domlo(2)) then

     do k = ulo(3), uhi(3)
      do i = ulo(1), uhi(1)
      ibc = min(i,domhi(1)+ng)
        if (bc_jlo_type(ibc,k,1) == FSW_) then
            drag_u(i,domlo(2),k) = drag_u(i,domlo(2),k) + drag_u(i, domlo(2)-1,k)
        end if
        drag_u(i,domlo(2)-1,k) = 1.d200
     end do
     end do
     do k = wlo(3), whi(3)
      kbc = min(k,domhi(3)+ng)
      do i = wlo(1), whi(1)
        if (bc_jlo_type(i,kbc,1) == FSW_) then
            drag_w(i,domlo(2),k) = drag_w(i,domlo(2),k) + drag_w(i,domlo(2)-1,k)
        end if
        drag_w(i,domlo(2)-1,k) = 1.d200
     end do
     end do

  endif

  if (vhi(2).ge.domhi(2)+1) then
     do k = ulo(3), uhi(3)
      do i = ulo(1), uhi(1)
      ibc = min(i,domhi(1)+ng)
        if (bc_jhi_type(ibc,k,1) == FSW_) then
            drag_u(i,domhi(2),k) = drag_u(i,domhi(2),k) + drag_u(i,domhi(2)+1,k)
        end if
        drag_u(i,domhi(2)+1,k) = 1.d200
     end do
     end do
     do k = wlo(3), whi(3)
      kbc = min(k,domhi(3)+ng)
      do i = wlo(1), whi(1)
        if (bc_jhi_type(i,kbc,1) == FSW_) then
            drag_w(i,domhi(2),k) = drag_w(i,domhi(2),k) + drag_w(i,domhi(2)+1,k)
        end if
        drag_w(i,domhi(2)+1,k) = 1.d200
     end do
     end do
  endif


  if (wlo(3).le.domlo(3)) then
     do j = ulo(2), uhi(2)
      do i = ulo(1), uhi(1)
      ibc = min(i,domhi(1)+ng)
        if (bc_klo_type(ibc,j,1) == FSW_) then
            drag_u(i,j,domlo(3)) = drag_u(i,j,domlo(3)) + drag_u(i,j,domlo(3)-1)
        end if
        drag_u(i,j,domlo(3)-1) = 1.d200
     end do
     end do
     do j = vlo(2), vhi(2)
      jbc = min(j,domhi(2)+ng)
      do i = vlo(1), vhi(1)
        if (bc_klo_type(i,jbc,1) == FSW_) then
            drag_v(i,j,domlo(3)) = drag_v(i,j,domlo(3)) + drag_v(i,j,domlo(3)-1)
        end if
        drag_v(i,j,domlo(3)-1) = 1.d200
     end do
     end do
  endif


  if (whi(3).ge.domhi(3)+1) then
     do j = ulo(2), uhi(2)
      do i = ulo(1), uhi(1)
      ibc = min(i,domhi(1)+ng)
        if (bc_khi_type(ibc,j,1) == FSW_) then
            drag_u(i,j,domhi(3)) = drag_u(i,j,domhi(3)) + drag_u(i,j,domhi(3)+1)
        end if
        drag_u(i,j,domhi(3)+1) = 1.d200
     end do
     end do
     do j = vlo(2), vhi(2)
      jbc = min(j,domhi(2)+ng)
      do i = vlo(1), vhi(1)
        if (bc_khi_type(i,jbc,1) == FSW_) then
            drag_v(i,j,domhi(3)) = drag_v(i,j,domhi(3)) + drag_v(i,j,domhi(3)+1)
        end if
        drag_v(i,j,domhi(3)+1) = 1.d200
     end do
     end do
  endif 

!--------------------------------------------------------------------//
!-  Now set drag on normal wall faces to zero
!--------------------------------------------------------------------//

  if (ulo(1).le.domlo(1)) then
     ilo = domlo(1)
     do k=ulo(3),uhi(3)
        do j=ulo(2),uhi(2)
           if (bc_ilo_type(j,k,1) == NSW_ .or. &
               bc_ilo_type(j,k,1) == PSW_ .or. &
               bc_ilo_type(j,k,1) == FSW_ .or. &
               bc_ilo_type(j,k,1) == PINF_ .or. &
               bc_ilo_type(j,k,1) == MINF_ .or. &
               bc_ilo_type(j,k,1) == POUT_) then

               drag_u(ilo,j,k) = 0.d0

           end if
        end do
     end do
  endif

  if (uhi(1).ge.domhi(1)+1) then
     ihi = domhi(1)+1
     do k=ulo(3),uhi(3)
        do j=ulo(2),uhi(2)
           if (bc_ihi_type(j,k,1) == NSW_ .or. &
               bc_ihi_type(j,k,1) == PSW_ .or. &
               bc_ihi_type(j,k,1) == FSW_ .or. &
               bc_ihi_type(j,k,1) == PINF_ .or. &
               bc_ihi_type(j,k,1) == MINF_ .or. &
               bc_ihi_type(j,k,1) == POUT_) then

               drag_u(ihi,j,k) = 0.d0

           end if
        end do
     end do
  endif

  if (vlo(2).le.domlo(2)) then
     jlo = domlo(2)
     do k=vlo(3),vhi(3)
        do i=vlo(1),vhi(1)
           if (bc_jlo_type(i,k,1) == NSW_ .or. &
               bc_jlo_type(i,k,1) == PSW_ .or. &
               bc_jlo_type(i,k,1) == FSW_ .or. &
               bc_jlo_type(i,k,1) == PINF_ .or. &
               bc_jlo_type(i,k,1) == MINF_ .or. &
               bc_jlo_type(i,k,1) == POUT_) then

               drag_v(i,jlo,k) = 0.d0

           end if
        end do
     end do
  endif

  if (vhi(2).ge.domhi(2)+1) then
     jhi = domhi(2)+1
     do k=vlo(3),vhi(3)
        do i=vlo(1),vhi(1)
           if (bc_jhi_type(i,k,1) == NSW_ .or. &
               bc_jhi_type(i,k,1) == PSW_ .or. &
               bc_jhi_type(i,k,1) == FSW_ .or. &
               bc_jhi_type(i,k,1) == PINF_ .or. &
               bc_jhi_type(i,k,1) == MINF_ .or. &
               bc_jhi_type(i,k,1) == POUT_) then

               drag_v(i,jhi,k) = 0.d0

           end if
        end do
     end do
  endif

  if (wlo(3).le.domlo(3)) then
     klo = domlo(3)
     do j=wlo(2),whi(2)
        do i=wlo(1),whi(1)
           if (bc_klo_type(i,j,1) == NSW_ .or. &
               bc_klo_type(i,j,1) == PSW_ .or. &
               bc_klo_type(i,j,1) == FSW_ .or. &
               bc_klo_type(i,j,1) == PINF_ .or. &
               bc_klo_type(i,j,1) == MINF_ .or. &
               bc_klo_type(i,j,1) == POUT_) then

               drag_w(i,j,klo) = 0.d0

           end if
        end do
     end do
  endif

  if (whi(3).ge.domhi(3)+1) then
     khi = domhi(3)+1
     do j=wlo(2),whi(2)
        do i=wlo(1),whi(1)
           if (bc_khi_type(i,j,1) == NSW_ .or. &
               bc_khi_type(i,j,1) == PSW_ .or. &
               bc_khi_type(i,j,1) == FSW_ .or. &
               bc_khi_type(i,j,1) == PINF_ .or. &
               bc_khi_type(i,j,1) == MINF_ .or. &
               bc_khi_type(i,j,1) == POUT_) then

               drag_w(i,j,khi) = 0.d0

           end if
        end do
     end do
  endif
end subroutine set_drag_bcs
