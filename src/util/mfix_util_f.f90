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
