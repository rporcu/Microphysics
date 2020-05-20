!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: WRITE_USER                                             !
!  Purpose: Write user-defined output                                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine collect_fluid (lo, hi, domlo, domhi,                    &
     &                    ep_g,    slo, shi,                       &
     &                    p_g,     plo, phi,                       &
     &                    vel_g,   vlo, vhi,                       &
     &                    vratio,  rlo, rhi,                       &
     &                    usr_x_w, usr_x_e,                        &
     &                    usr_y_s, usr_y_n,                        &
     &                    usr_z_b, usr_z_t, dx,                    &
     &                    sum_ep_g, sum_p_g,  sum_vol,             &
     &                    sum_velx, sum_vely, sum_velz)            &
     & bind(C, name="mfix_collect_fluid")

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int

  implicit none

! Dummy arguments ....................................................//
  integer(c_int), intent(in   ) ::    lo(3),    hi(3)
  integer(c_int), intent(in   ) :: domlo(3), domhi(3)
  integer(c_int), intent(in   ) ::   slo(3),   shi(3)
  integer(c_int), intent(in   ) ::   plo(3),   phi(3)
  integer(c_int), intent(in   ) ::   vlo(3),   vhi(3)
  integer(c_int), intent(in   ) ::   rlo(3),   rhi(3)

  real(rt),       intent(in   ) :: usr_x_w, usr_x_e
  real(rt),       intent(in   ) :: usr_y_s, usr_y_n
  real(rt),       intent(in   ) :: usr_z_b, usr_z_t

  real(rt),       intent(in   ) ::  dx(3),                         &
       &    ep_g(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3)),     &
       &     p_g(plo(1):shi(1), plo(2):phi(2), plo(3):phi(3)),     &
       &   vel_g(vlo(1):vhi(1), vlo(2):vhi(2), vlo(3):vhi(3) ,3),  &
       &  vratio(rlo(1):rhi(1), rlo(2):rhi(2), rlo(3):shi(3))

  real(rt),   intent(inout) :: sum_ep_g
  real(rt),   intent(inout) :: sum_p_g
  real(rt),   intent(inout) :: sum_vol
  real(rt),   intent(inout) :: sum_velx, sum_vely, sum_velz

! Local variables .....................................................//
  integer  :: i, j, k, is, ie, js, je, ks, ke

  real(rt) :: regular_cell_volume, this_cell_volume

  regular_cell_volume = dx(1) * dx(2) * dx(3)

  is = max(lo(1), domlo(1), floor(usr_x_w/dx(1) + 0.5d0))
  ie = min(hi(1), domhi(1), floor(usr_x_e/dx(1) + 0.5d0))

  js = max(lo(2), domlo(2), floor(usr_y_s/dx(2) + 0.5d0))
  je = min(hi(2), domhi(2), floor(usr_y_n/dx(2) + 0.5d0))

  ks = max(lo(3), domlo(3), floor(usr_z_b/dx(3) + 0.5d0))
  ke = min(hi(3), domhi(3), floor(usr_z_t/dx(3) + 0.5d0))


  do k=ks,ke
     do j=js,je
        do i=is,ie

           this_cell_volume = regular_cell_volume * vratio(i,j,k)

           sum_vol  = sum_vol  + this_cell_volume

           sum_ep_g = sum_ep_g + this_cell_volume * ep_g(i,j,k)
           sum_p_g  = sum_p_g  + this_cell_volume *  p_g(i,j,k)

           sum_velx = sum_velx + this_cell_volume * vel_g(i,j,k,1)
           sum_vely = sum_vely + this_cell_volume * vel_g(i,j,k,2)
           sum_velz = sum_velz + this_cell_volume * vel_g(i,j,k,3)

        enddo
     enddo
  enddo

end subroutine collect_fluid
