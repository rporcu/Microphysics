!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_solids_volume                                      !
!                                                                      !
!  Purpose: Map particle volume to Eulerian grid. The result is stored !
!  in the gas phase volume fraction variable, ep_g.                    !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
subroutine calc_solids_volume(slo, shi, max_pip, ep_g, &
  pstate, pvol, ppos, dx, dy, dz) &
  bind(C, name="calc_solids_volume")

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int

  use discretelement, only: normal_particle

  implicit none

  integer(c_int), intent(in   ) :: slo(3),shi(3)
  integer(c_int), intent(in   ) :: max_pip

  integer(c_int), intent(in   ) :: pstate(max_pip)
  real(c_real),   intent(in   ) :: ppos(max_pip,3)
  real(c_real),   intent(in   ) :: pvol(max_pip)

  real(c_real),   intent(in   ) :: dx, dy, dz

  real(c_real),   intent(inout) :: ep_g&
    (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))


! Local variables
!-----------------------------------------------
  ! Fluid cell index
  integer :: i, j, k, np
  ! One over cell volume
  real(c_real) :: odx, ody, odz

  odx = 1.0d0/dx
  ody = 1.0d0/dy
  odz = 1.0d0/dz

! Calculate the gas phae forces acting on each particle.
  do np=1,max_pip
     if(pstate(np) == normal_particle) then

! Fluid cell containing the particle
        i = floor(ppos(np,1)*odx)
        j = floor(ppos(np,2)*ody)
        k = floor(ppos(np,3)*odz)

        ep_g(i,j,k) = ep_g(i,j,k) + pvol(np)

     endif
  enddo

end subroutine calc_solids_volume


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_volume_fraction                                    !
!                                                                      !
!  Purpose: Convert the solids volume to the fluid volume fraction.    !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
subroutine calc_volume_fraction(slo, shi, lo, hi, &
  ep_g, rop_g, ro_g, dx, dy, dz) &
  bind(C, name="calc_volume_fraction")

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int

  implicit none

  integer(c_int), intent(in   ) :: slo(3),shi(3)
  integer(c_int), intent(in   ) ::  lo(3), hi(3)

  real(c_real),   intent(in   ) :: dx, dy, dz

  real(c_real),   intent(in   ) :: ro_g&
    (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

  real(c_real),   intent(inout) :: ep_g&
    (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
  real(c_real),   intent(  out) :: rop_g&
    (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))


! Local variables
!-----------------------------------------------
  ! Fluid cell index
  integer :: i, j, k
  ! One over cell volume
  real(c_real) :: oovol

  oovol = 1.0d0/(dx*dy*dz)

  do k = lo(3),hi(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)
           ep_g(i,j,k) = 1.0d0 - ep_g(i,j,k)*oovol
           rop_g(i,j,k) = ep_g(i,j,k) * ro_g(i,j,k)
        enddo
     enddo
  enddo

end subroutine calc_volume_fraction
