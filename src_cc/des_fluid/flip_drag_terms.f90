!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: flip_drag_terms                                         C
!                                                                      C
!  Purpose: This subroutine "flips" any part of "beta" (aka beta)     C
!           or "beta*vpart" (aka drag)                                 C
!           deposited outside the domain at a FSW or NSW               C
!           back into the domain.                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
subroutine flip_drag_terms(beta, blo, bhi, &
                           drag, dlo, dhi, &
                           bct_ilo, bct_ihi, &
                           bct_jlo, bct_jhi, &
                           bct_klo, bct_khi, domlo, domhi, ng ) &
     bind(C, name="flip_drag_terms")

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int

  use bc, only: NSW_, FSW_, cyclic_x, cyclic_y, cyclic_z

  implicit none

  integer(c_int), intent(in   ) :: dlo(3),dhi(3)
  integer(c_int), intent(in   ) :: blo(3),bhi(3)
  integer(c_int), intent(in   ) :: domlo(3),domhi(3), ng

  real(rt), intent(inout) :: &
         beta(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3)  ), &
         drag(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),3)

  integer(c_int), intent(in   ) :: &
       bct_ilo(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
       bct_ihi(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
       bct_jlo(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
       bct_jhi(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
       bct_klo(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
       bct_khi(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

  integer(c_int) :: i,j,k

  ! ***************************************************************************
  ! Note that the algorithm below actually correctly accounts for beta and drag 
  !      deposited outside of edges and corners because of the cumulative effect.
  !      For example, if there was some beta or drag at the (-1,-1,-1) corner, 
  !      it would first get added into (0,-1,-1), then that would get added 
  !      into (0,0,-1), which would finally get added into (0,0,0).
  ! ***************************************************************************

  if (dlo(1).le.domlo(1) .and. .not.cyclic_x) then

     do k = dlo(3), dhi(3)
      do j = dlo(2), dhi(2)
        if (bct_ilo(j,k,1) == FSW_ .or. bct_ilo(j,k,1) == NSW_) then
            beta(domlo(1),j,k  ) = beta(domlo(1),j,k  ) + beta(domlo(1)-1,j,k  )
            drag(domlo(1),j,k,1) = drag(domlo(1),j,k,1) + drag(domlo(1)-1,j,k,1)
            drag(domlo(1),j,k,2) = drag(domlo(1),j,k,2) + drag(domlo(1)-1,j,k,2)
            drag(domlo(1),j,k,3) = drag(domlo(1),j,k,3) + drag(domlo(1)-1,j,k,3)
        end if

           beta(domlo(1)-1,j,k  ) = 1.d200
           drag(domlo(1)-1,j,k,:) = 1.d200

     end do
     end do

  endif

  if (dhi(1).ge.domhi(1)+1 .and. .not.cyclic_x) then

     do k = dlo(3), dhi(3)
      do j = dlo(2), dhi(2)
        if (bct_ihi(j,k,1) == FSW_ .or. bct_ihi(j,k,1) == NSW_) then
            beta(domhi(1),j,k  ) = beta(domhi(1),j,k  ) + beta(domhi(1)+1,j,k  )
            drag(domhi(1),j,k,1) = drag(domhi(1),j,k,1) + drag(domhi(1)+1,j,k,1)
            drag(domhi(1),j,k,2) = drag(domhi(1),j,k,2) + drag(domhi(1)+1,j,k,2)
            drag(domhi(1),j,k,3) = drag(domhi(1),j,k,3) + drag(domhi(1)+1,j,k,3)
        end if

           beta(domhi(1)+1,j,k  ) = 1.d200
           drag(domhi(1)+1,j,k,:) = 1.d200

     end do
     end do

  endif

  if (dlo(2).le.domlo(2) .and. .not.cyclic_y) then

     do k = dlo(3), dhi(3)
      do i = dlo(1), dhi(1)
        if (bct_jlo(i,k,1) == FSW_ .or. bct_jlo(i,k,1) == NSW_) then
            beta(i,domlo(2),k  ) = beta(i,domlo(2),k  ) + beta(i,domlo(2)-1,k  )
            drag(i,domlo(2),k,1) = drag(i,domlo(2),k,1) + drag(i,domlo(2)-1,k,1)
            drag(i,domlo(2),k,2) = drag(i,domlo(2),k,2) + drag(i,domlo(2)-1,k,2)
            drag(i,domlo(2),k,3) = drag(i,domlo(2),k,3) + drag(i,domlo(2)-1,k,3)
        end if

           beta(i,domlo(2)-1,k  ) = 1.d200
           drag(i,domlo(2)-1,k,:) = 1.d200

     end do
     end do

  endif

  if (dhi(2).ge.domhi(2)+1 .and. .not.cyclic_y) then

     do k = dlo(3), dhi(3)
      do i = dlo(1), dhi(1)
        if (bct_jhi(i,k,1) == FSW_ .or. bct_jhi(i,k,1) == NSW_) then
            beta(i,domhi(2),k  ) = beta(i,domhi(2),k  ) + beta(i,domhi(2)+1,k  )
            drag(i,domhi(2),k,1) = drag(i,domhi(2),k,1) + drag(i,domhi(2)+1,k,1)
            drag(i,domhi(2),k,2) = drag(i,domhi(2),k,2) + drag(i,domhi(2)+1,k,2)
            drag(i,domhi(2),k,3) = drag(i,domhi(2),k,3) + drag(i,domhi(2)+1,k,3)
        end if

           beta(i,domhi(2)+1,k  ) = 1.d200
           drag(i,domhi(2)+1,k,:) = 1.d200

     end do
     end do

  endif

  if (dlo(3).le.domlo(3) .and. .not.cyclic_z) then

     do j = dlo(2), dhi(2)
      do i = dlo(1), dhi(1)
        if (bct_klo(i,j,1) == FSW_ .or. bct_klo(i,j,1) == NSW_) then
            beta(i,j,domlo(3)  ) = beta(i,j,domlo(3)  ) + beta(i,j,domlo(3)-1  )
            drag(i,j,domlo(3),1) = drag(i,j,domlo(3),1) + drag(i,j,domlo(3)-1,1)
            drag(i,j,domlo(3),2) = drag(i,j,domlo(3),2) + drag(i,j,domlo(3)-1,2)
            drag(i,j,domlo(3),3) = drag(i,j,domlo(3),3) + drag(i,j,domlo(3)-1,3)
        end if

           beta(i,j,domlo(3)-1  ) = 1.d200
           drag(i,j,domlo(3)-1,:) = 1.d200

     end do
     end do

  endif

  if (dhi(3).ge.domhi(3)+1 .and. .not.cyclic_z) then

     do j = dlo(2), dhi(2)
      do i = dlo(1), dhi(1)

        if (bct_khi(i,j,1) == FSW_ .or. bct_khi(i,j,1) == NSW_) then
            beta(i,j,domhi(3)  ) = beta(i,j,domhi(3)  ) + beta(i,j,domhi(3)+1  )
            drag(i,j,domhi(3),1) = drag(i,j,domhi(3),1) + drag(i,j,domhi(3)+1,1)
            drag(i,j,domhi(3),2) = drag(i,j,domhi(3),2) + drag(i,j,domhi(3)+1,2)
            drag(i,j,domhi(3),3) = drag(i,j,domhi(3),3) + drag(i,j,domhi(3)+1,3)
        end if

           beta(i,j,domhi(3)+1  ) = 1.d200
           drag(i,j,domhi(3)+1,:) = 1.d200
     end do
     end do

  endif
end subroutine flip_drag_terms
