!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_FLOW                                             !
!                                                                      !
!  Purpose: Check boundary condition specifications                    !
!     - convert physical locations to i, j, k's                        !
!     - convert mass and volumetric flows to velocities (FLOW_TO_VEL)  !
!     - check specification of physical quantities                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine set_bc_flow(xlength, ylength, zlength, dx, dy, dz) &
     bind(C,name ="set_bc_flow")

     use constant,      only: mmax
     use param1,        only: zero, one, undefined, equal
     use param,         only: dimension_bc, dim_m
     use bc,            only: bc_defined, bc_type, bc_rop_s, bc_ep_s, bc_ep_g

    use flow_to_vel_new_module, only: flow_to_vel
    use amrex_fort_module, only : c_real => amrex_real

    implicit none

    real(c_real)  , intent(in) :: xlength, ylength, zlength
    real(c_real)  , intent(in) :: dx, dy, dz

    integer :: bcv, i, mmax_tot
    logical :: skip(1:dim_m)     ! Flag to skip checks on indexed solid phase.

    ! Total number of solids.
    mmax_tot = mmax

    ! Loop over each defined BC and check the user data.
    do bcv = 1, dimension_bc

       if(.not.bc_defined(bcv)) cycle

       ! Determine which solids phases are present.
       skip = .false.
       do i = 1, dim_m
          if ((equal(bc_rop_s(bcv,i), undefined).or.equal(bc_rop_s(bcv,i), zero)) &
               .and.(equal(bc_ep_s(bcv,i), undefined).or.equal(bc_ep_s(bcv,i), zero))) then
             skip = .true.
          endif
       end do

       if(mmax_tot == 1 .and. .not.equal(bc_ep_g(bcv), one)) skip(1) = .false.

       select case (trim(BC_TYPE(BCV)))

       case ('MASS_INFLOW','MI','MASS_OUTFLOW','MO')
          call flow_to_vel(mmax_tot, skip, bcv, xlength, ylength, zlength, dx, dy, dz)
       end select
    enddo

  end subroutine set_bc_flow
