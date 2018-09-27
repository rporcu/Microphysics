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
subroutine set_drag_bcs(drag, dlo, dhi, &
                        bc_ilo_type, bc_ihi_type, &
                        bc_jlo_type, bc_jhi_type, &
                        bc_klo_type, bc_khi_type, domlo, domhi, ng ) &
     bind(C, name="set_drag_bcs")

     use amrex_fort_module, only : rt => amrex_real
     use iso_c_binding , only: c_int

     use bc, only: nsw_, fsw_, psw_, pinf_, pout_, minf_

     implicit none

     integer(c_int), intent(in   ) :: dlo(3),dhi(3)
     integer(c_int), intent(in   ) :: domlo(3),domhi(3), ng

     real(rt), intent(inout) ::  drag&
        (dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),3)

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

!--------------------------------------------------------------------//
!-  Do tangential corrections at slip walls only
!--------------------------------------------------------------------//

  if (dlo(1).le.domlo(1)) then
     do k = dlo(3), dhi(3)
      do j = dlo(2), dhi(2)
        if (bc_ilo_type(j,k,1) == FSW_) then
            drag(domlo(1),j,k,2) = drag(domlo(1),j,k,2) + drag(domlo(1)-1,j,k,2)
            drag(domlo(1),j,k,3) = drag(domlo(1),j,k,3) + drag(domlo(1)-1,j,k,3)
        end if
     end do
     end do
  endif

  if (dlo(1).lt.domlo(1)) then
     do k = dlo(3), dhi(3)
        do j = dlo(2), dhi(2)
           drag(domlo(1)-1,j,k,2) = 1.d200
           drag(domlo(1)-1,j,k,3) = 1.d200
        end do
     end do
  endif

  if (dhi(1).ge.domhi(1)+1) then
     do k = dlo(3), dhi(3)
      do j = dlo(2), dhi(2)
        if (bc_ihi_type(j,k,1) == FSW_) then
            drag(domhi(1),j,k,2) = drag(domhi(1),j,k,2) + drag(domhi(1)+1,j,k,2)
            drag(domhi(1),j,k,3) = drag(domhi(1),j,k,3) + drag(domhi(1)+1,j,k,3)
        end if
     end do
     end do
  endif

  if (dhi(1).gt.domhi(1)+1) then
     do k = dlo(3), dhi(3)
        do j = dlo(2), dhi(2)
           drag(domhi(1)+1,j,k,2) = 1.d200
           drag(domhi(1)+1,j,k,3) = 1.d200
        end do
     end do
  endif

  if (dlo(2).le.domlo(2)) then
     do k = dlo(3), dhi(3)
      do i = dlo(1), dhi(1)
        if (bc_jlo_type(i,k,1) == FSW_) then
            drag(i,domlo(2),k,1) = drag(i,domlo(2),k,1) + drag(i,domlo(2)-1,k,1)
            drag(i,domlo(2),k,3) = drag(i,domlo(2),k,3) + drag(i,domlo(2)-1,k,3)
        end if
     end do
     end do
  endif

  if (dlo(2).lt.domlo(2)) then
     do k = dlo(3), dhi(3)
        do i = dlo(1), dhi(1)
           drag(i,domlo(2)-1,k,1) = 1.d200
           drag(i,domlo(2)-1,k,3) = 1.d200
        end do
     end do
  endif

  if (dhi(2).ge.domhi(2)+1) then
     do k = dlo(3), dhi(3)
      do i = dlo(1), dhi(1)
        if (bc_jhi_type(i,k,1) == FSW_) then
            drag(i,domhi(2),k,1) = drag(i,domhi(2),k,1) + drag(i,domhi(2)+1,k,1)
            drag(i,domhi(2),k,3) = drag(i,domhi(2),k,3) + drag(i,domhi(2)+1,k,3)
        end if
     end do
     end do
  endif

  if (dhi(2).gt.domhi(2)+1) then
     do k = dlo(3), dhi(3)
        do i = dlo(1), dhi(1)
           drag(i,domhi(2)+1,k,1) = 1.d200
           drag(i,domhi(2)+1,k,3) = 1.d200
        end do
     end do
  endif

  if (dlo(3).le.domlo(3)) then
     do j = dlo(2), dhi(2)
      do i = dlo(1), dhi(1)
        if (bc_klo_type(i,j,1) == FSW_) then
            drag(i,j,domlo(3),1) = drag(i,j,domlo(3),1) + drag(i,j,domlo(3)-1,1)
            drag(i,j,domlo(3),2) = drag(i,j,domlo(3),2) + drag(i,j,domlo(3)-1,2)
        end if
     end do
     end do
  endif

  if (dlo(3).lt.domlo(3)) then
     do j = dlo(2), dhi(2)
        do i = dlo(1), dhi(1)
           drag(i,j,domlo(3)-1,1) = 1.d200
           drag(i,j,domlo(3)-1,2) = 1.d200
        end do
     end do
  endif

  if (dhi(3).ge.domhi(3)+1) then
     do j = dlo(2), dhi(2)
      do i = dlo(1), dhi(1)
        if (bc_khi_type(i,j,1) == FSW_) then
            drag(i,j,domhi(3),1) = drag(i,j,domhi(3),1) + drag(i,j,domhi(3)+1,1)
            drag(i,j,domhi(3),2) = drag(i,j,domhi(3),2) + drag(i,j,domhi(3)+1,2)
        end if
     end do
     end do
  endif

  if (dhi(3).gt.domhi(3)+1) then
     do j = dlo(2), dhi(2)
        do i = dlo(1), dhi(1)
           drag(i,j,domhi(3)+1,1) = 1.d200
           drag(i,j,domhi(3)+1,2) = 1.d200
        end do
     end do
  endif
end subroutine set_drag_bcs
