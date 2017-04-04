module check_domain_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: check_domain                                              !
!  Purpose: read and verify input data, open files                     !
!                                                                      !
!  Author: P. Nicoletti                               Date: 04-DEC-91  !
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine check_domain(dx,dy,dz,domlo,domhi) &
          bind(C, name="check_domain")

       use get_bc_area_module, only: get_bc_area
       use set_bc_flow_module, only: set_bc_flow
       use check_boundary_conditions_module, only: check_boundary_conditions
       use check_initial_conditions_module, only: check_initial_conditions
       use check_point_sources_module, only: check_point_sources

      integer(c_int), intent(in) :: domlo(3),domhi(3)
      real(c_real)  , intent(in) :: dx, dy ,dz

      call check_initial_conditions(dx,dy,dz,domlo,domhi)
      call check_boundary_conditions(dx,dy,dz,domlo,domhi)
      call check_point_sources(dx,dy,dz)

      ! Compute area of boundary surfaces.
      call get_bc_area(dx,dy,dz)

      ! Convert (mass, volume) flows to velocities.
      call set_bc_flow

      end subroutine check_domain

end module check_domain_module
