!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: WRITE_USER                                             !
!  Purpose: Write user-defined output                                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine collect_fluid(slo, shi, lo, hi, domlo, domhi, p_g, ep_g, &
     dx, dy, dz, sums) bind(C, name="mfix_collect_fluid")

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int

  use param, only: is_defined, dim_usr

  use output, only: usr_dt
  use output, only: usr_x_w, usr_y_s, usr_z_b
  use output, only: usr_x_e, usr_y_n, usr_z_t

  implicit none

! Dummy arguments ....................................................//
  integer(c_int), intent(in   ) ::   slo(3),   shi(3)
  integer(c_int), intent(in   ) ::    lo(3),    hi(3)
  integer(c_int), intent(in   ) :: domlo(3), domhi(3)

  real(c_real),   intent(in   ) ::  dx, dy, dz, &
        p_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
       ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

  real(c_real),   intent(inout) :: sums(dim_usr,2)

! Local variables .....................................................//
  integer :: i,j,k, lc, is, ie, js, je, ks, ke

  do lc=1, dim_usr

     if(is_defined(usr_dt(lc))) then

        is = max(lo(1), domlo(1), floor(usr_x_w(lc)/dx + 0.5d0))
        ie = min(hi(1), domhi(1), floor(usr_x_e(lc)/dx + 0.5d0))

        js = max(lo(2), domlo(2), floor(usr_y_s(lc)/dy + 0.5d0))
        je = min(hi(2), domhi(2), floor(usr_y_n(lc)/dy + 0.5d0))

        ks = max(lo(3), domlo(3), floor(usr_z_b(lc)/dz + 0.5d0))
        ke = min(hi(3), domhi(3), floor(usr_z_t(lc)/dz + 0.5d0))

        sums(lc,1) = sums(lc,1) + sum( p_g(is:ie,js:je,ks:ke))
        sums(lc,2) = sums(lc,2) + sum(ep_g(is:ie,js:je,ks:ke))

     endif
  enddo

end subroutine collect_fluid

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: WRITE_USER                                             !
!  Purpose: Write user-defined output                                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine write_fluid(domlo, domhi, dx, dy, dz, time, dt, sums) &
     bind(C, name="mfix_write_fluid")

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int

  use param, only: is_defined, dim_usr

  use output, only: newUnit
  use output, only: usr_dt
  use output, only: usr_x_w, usr_y_s, usr_z_b
  use output, only: usr_x_e, usr_y_n, usr_z_t

  implicit none

! Dummy arguments ....................................................//
  integer(c_int), intent(in   ) :: domlo(3), domhi(3)
  real(c_real),   intent(in   ) :: dx, dy, dz, time, dt, sums(dim_usr,2)


! Local variables .....................................................//
  real(c_real), save :: next_time(dim_usr) = 0.0d0;
  real(c_real) :: cells, avgPg, avgEPg

  integer :: lc, is, ie, js, je, ks, ke, lunit
  character :: clc

  do lc=1, dim_usr

     if(is_defined(usr_dt(lc))) then

        if(time + 0.1d0*dt >= next_time(lc)) then

           next_time(lc) = time + usr_dt(lc)

           is = max(domlo(1), floor(usr_x_w(lc)/dx + 0.5d0))
           ie = min(domhi(1), floor(usr_x_e(lc)/dx + 0.5d0))

           js = max(domlo(2), floor(usr_y_s(lc)/dy + 0.5d0))
           je = min(domhi(2), floor(usr_y_n(lc)/dy + 0.5d0))

           ks = max(domlo(3), floor(usr_z_b(lc)/dz + 0.5d0))
           ke = min(domhi(3), floor(usr_z_t(lc)/dz + 0.5d0))

           cells = real((ie-is+1)*(je-js+1)*(ke-ks+1))

           if(cells > 0.0) then
              avgPg   = sums(lc,1)/cells
              avgEPg  = sums(lc,2)/cells
           else
              avgPg  = 0.0
              avgEPg = 0.0
           endif

           lunit = newUnit()
           write(clc,"(i1)") lc
           open(unit=lunit, file='usr_'//clc//'.csv', &
                status='unknown', position='append')
           write(lunit,"(es15.8,2(',',2x,es15.8))") time, avgPg, avgEPg
           close(lunit)
        endif
     endif
  enddo

end subroutine write_fluid
