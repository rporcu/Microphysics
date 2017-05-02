MODULE CFNEWVALUES_MODULE

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CFNEWVALUES                                            !
!                                                                      !
!  Purpose: DES - Calculate the new values of particle velocity,       !
!           position, angular velocity etc                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  SUBROUTINE CFNEWVALUES(max_pip, particle_state, pmass,&
    omoi, drag, des_pos_new, des_vel_new, omega_new, fc, tow, &
    des_acc_old, rot_acc_old)

    use discretelement, only: dtsolid
    use discretelement, only: normal_particle, exiting_particle
    use param, only: zero
    use constant, only: gravity

    implicit none

    integer     , intent(in   ) :: max_pip
    integer     , intent(in   ) :: particle_state(max_pip)
    real(c_real), intent(in   ) :: omoi(max_pip)
    real(c_real), intent(in   ) :: pmass(max_pip)
    real(c_real), intent(in   ) :: drag(max_pip,3)
    real(c_real), intent(inout) :: des_pos_new(max_pip,3)
    real(c_real), intent(inout) :: des_vel_new(max_pip,3)
    real(c_real), intent(inout) :: omega_new(max_pip,3)
    real(c_real), intent(inout) :: fc(max_pip,3)
    real(c_real), intent(inout) :: tow(max_pip,3)
    real(c_real), intent(inout) :: des_acc_old(max_pip,3)
    real(c_real), intent(inout) :: rot_acc_old(max_pip,3)

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
    integer :: L
    real(c_real) :: DD(3)
    logical, save :: first_pass = .true.

    real(c_real) :: lVELo(3), lPOSo(3)

!-----------------------------------------------

! Adams-Bashforth defaults to Euler for the first time step.
    if(first_pass) then
       first_pass = .false.
       do l =1, max_pip
          if(particle_state(l) == normal_particle .or.      &
             particle_state(l) == exiting_particle) then
             des_acc_old(l,:) = (fc(l,:) + drag(l,:))/pmass(l) + gravity(:)
             rot_acc_old(l,:) = tow(l,:)
          endif
       enddo
    endif

    do l = 1, max_pip
       if(particle_state(l) == normal_particle .or.      &
          particle_state(l) == exiting_particle) then

          fc(l,:) = (fc(l,:) + drag(l,:))/pmass(l) + gravity(:)

! Advance particle position, velocity
          lvelo = des_vel_new(l,:)
          lposo = des_pos_new(l,:)

! Second-order Adams-Bashforth/Trapezoidal scheme
          des_vel_new(l,:) = lvelo(:) + 0.5d0*&
            ( 3.d0*fc(l,:)-des_acc_old(l,:) )*dtsolid

          omega_new(l,:)   =  omega_new(l,:) + 0.5d0*&
            ( 3.d0*tow(l,:)*omoi(l)-rot_acc_old(l,:) )*dtsolid

          dd(:) = 0.5d0*( lvelo(:)+des_vel_new(l,:) )*dtsolid

          des_pos_new(l,:) = lposo(:) + dd(:)
          des_acc_old(l,:) = fc(l,:)
          rot_acc_old(l,:) = tow(l,:)*omoi(l)

! Reset total contact force and torque
          fc(l,:) = zero
          tow(l,:) = zero
       endif
    enddo

    return
  end subroutine cfnewvalues
end module cfnewvalues_module
