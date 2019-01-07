   !
   ! Average to faces in chosen direction  -- note we only average the "idir"th 
   !    component of cc onto the idir'th face
   ! 
   subroutine average_cc_to_fc ( lo, hi, fx, fxlo, fxhi, fy, fylo, fyhi, fz, fzlo, fzhi,  &
                                 cc, slo, shi) bind(C) 

      use amrex_fort_module,  only: ar => amrex_real
      use iso_c_binding ,     only: c_int
      use param            ,  only: half

      ! Loop bounds (assumed face centered!)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Arrays bounds
      integer(c_int), intent(in   ) :: fxlo(3),fxhi(3)
      integer(c_int), intent(in   ) :: fylo(3),fyhi(3)
      integer(c_int), intent(in   ) :: fzlo(3),fzhi(3)
      integer(c_int), intent(in   ) :: slo(3),shi(3)

      ! Array
      real(ar),       intent(inout) :: &
           fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3)), &
           fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3)), &
           fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3))

      real(ar),       intent(in   ) :: &
           cc(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)
      
      ! Local variables
      integer  :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)+1
               fx(i,j,k) = half * ( cc(i-1,j,k,1) + cc(i,j,k,1) )  
            end do
         end do
      end do

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)+1
            do i = lo(1), hi(1)
               fy(i,j,k) = half * ( cc(i,j-1,k,2) + cc(i,j,k,2) )  
            end do
         end do
      end do

      do k = lo(3), hi(3)+1
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               fz(i,j,k) = half * ( cc(i,j,k-1,3) + cc(i,j,k,3) )  
            end do
         end do
      end do

   end subroutine average_cc_to_fc

  subroutine mfix_sum_mf(lo,hi,rho,r_lo,r_hi,dx, &
                         vol,v_lo,v_hi,mass) bind(c,name='mfix_sum_mf')

    use amrex_fort_module, only: rt => amrex_real, amrex_add

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: r_lo(3), r_hi(3)
    integer,  intent(in   ) :: v_lo(3), v_hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ) :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt), intent(in   ) :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    real(rt), intent(inout) :: mass

    integer  :: i, j, k
    real(rt) :: dm

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             dm = rho(i,j,k) * vol(i,j,k)

             call amrex_add(mass, dm)

          enddo
       enddo
    enddo

  end subroutine mfix_sum_mf
