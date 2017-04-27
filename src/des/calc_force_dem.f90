module calc_force_dem_module

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int

  implicit none
  private

  public calc_force_dem

contains
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  subroutine: calc_force_dem                                          !
  !  author: jay boyalakuntla                           date: 12-jun-04  !
  !                                                                      !
  !  purpose: calculate contact force and torque on particle from        !
  !           particle-particle and particle-wall collisions. treats     !
  !           wall interaction also as a two-particle interaction but    !
  !           accounting for the wall properties                         !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine calc_force_dem(phase, radius, pos, vel, omega, state, fc, tow)

    use cfrelvel_module, only: cfrelvel
    use discretelement, only: des_coll_model_enum, dtsolid
    use discretelement, only: des_etan, des_etat, hert_kt, hert_kn
    use discretelement, only: s_time, des_crossprdct
    use discretelement, only: kn, kt, mew, hertzian

    use drag_gs_des1_module, only: drag_gs_des
    use error_manager, only: init_err_msg, flush_err_msg, err_msg, ival


    integer, intent(in) :: phase(:)

    real(c_real), intent(in) :: radius(:)
    real(c_real), intent(in) :: pos(:,:)
    real(c_real), intent(in) :: vel(:,:)
    real(c_real), intent(in) :: omega(:,:)
    integer, intent(in)      :: state(:)
    ! integer, intent(inout)   :: pairs(:,:)
    ! integer, intent(inout)   :: pair_count
    real(c_real), intent(inout) :: fc(:,:), tow(:,:)

    ! local variables
    !---------------------------------------------------------------------//
    ! percent of particle radius when excess overlap will be flagged
    real(c_real), parameter :: flag_overlap = 0.20d0
    ! particle no. indices
    integer :: i, ll, cc
    ! the overlap occuring between particle-particle or particle-wall
    ! collision in the normal direction
    real(c_real) :: overlap_n, overlap_t(3)
    ! square root of the overlap
    real(c_real) :: sqrt_overlap
    ! distance vector between two particle centers or between a particle
    ! center and wall when the two surfaces are just at contact (i.e. no
    ! overlap)
    real(c_real) :: r_lm,dist_ci,dist_cl
    ! the normal and tangential components of the translational relative
    ! velocity
    real(c_real) :: v_rel_trans_norm, rad
    ! distance vector between two particle centers or between a particle
    ! center and wall at current and previous time steps
    real(c_real) :: dist(3), normal(3), dist_mag, pos_tmp(3)
    ! tangent to the plane of contact at current time step
    real(c_real) :: vrel_t(3)
    ! normal and tangential forces
    real(c_real) :: fn(3), ft(3)
    ! temporary storage of force
    real(c_real) :: fc_tmp(3)
    ! temporary storage of force for torque
    real(c_real) :: tow_force(3)
    ! temporary storage of torque
    real(c_real) :: tow_tmp(3,2)

    ! store solids phase index of particle (i.e. phase(np))
    integer :: phasei, phasell
    ! local values used spring constants and damping coefficients
    real(c_real) :: etan_des, etat_des
    real(c_real) :: kn_des, kt_des

    logical, parameter :: report_excess_overlap = .false.

    real(c_real) :: fnmd, mag_overlap_t, tangent(3)

    integer      :: pairs(6*size(state),3), pair_count
    !-----------------------------------------------

    ! check particle ll neighbor contacts
    !---------------------------------------------------------------------//

    ! calculate pairs of colliding particles
    pairs      = 0
    pair_count = 0
    call calc_collisions( size(state), pairs, pair_count, state, radius, pos)


    do cc = 1, pair_count

       ll = pairs(cc, 1)
       i = pairs(cc, 2)

       pos_tmp = pos(ll,:)
       rad     = radius(ll)

       r_lm = rad + radius(i)
       dist(:) = pos(i,:) - pos_tmp(:)
       dist_mag = dot_product(dist,dist)

       if(abs(dist_mag) < epsilon(dist_mag)) then
          write(*,8550) ll, i
          stop "division by zero"
8550      format('distance between particles is zero:',2(2x,i10))
       endif

       dist_mag = sqrt(dist_mag)
       normal(:)= dist(:)/dist_mag

       ! calcuate the normal overlap
       overlap_n = r_lm-dist_mag
       if(report_excess_overlap) call print_excess_overlap

       ! calculate the components of translational relative velocity for a
       ! contacting particle pair and the tangent to the plane of contact
       call cfrelvel(ll, i, v_rel_trans_norm, vrel_t,            &
            normal(:), dist_mag, vel, radius, omega)

       phasell = phase(ll)
       phasei = phase(i)

       ! hertz spring-dashpot contact model
       if (des_coll_model_enum .eq. hertzian) then
          sqrt_overlap = sqrt(overlap_n)
          kn_des = hert_kn(phasell,phasei)*sqrt_overlap
          kt_des = hert_kt(phasell,phasei)*sqrt_overlap
          sqrt_overlap = sqrt(sqrt_overlap)
          etan_des = des_etan(phasell,phasei)*sqrt_overlap
          etat_des = des_etat(phasell,phasei)*sqrt_overlap

          ! linear spring-dashpot contact model
       else
          kn_des = kn
          kt_des = kt
          etan_des = des_etan(phasell,phasei)
          etat_des = des_etat(phasell,phasei)
       endif

       ! calculate the normal contact force
       fn(:) =  -(kn_des * overlap_n * normal(:) + &
            etan_des * v_rel_trans_norm * normal(:))

       ! calcuate the tangential overlap
       overlap_t(:) = dtsolid*vrel_t(:)
       mag_overlap_t = sqrt(dot_product(overlap_t,overlap_t))

       ! calculate the tangential contact force.
       if(mag_overlap_t > 0.0) then
          ! max force before the on set of frictional slip.
          fnmd = mew*sqrt(dot_product(fn,fn))
          ! direction of tangential force.
          tangent = overlap_t/mag_overlap_t
          ! frictional slip
          ft = -fnmd * tangent
       else
          ft = 0.0
       endif

       ! calculate the distance from the particles' centers to the contact point,
       ! which is taken as the radical line
       ! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2
       dist_cl = dist_mag/2.d0 + (radius(ll)**2 - &
            radius(i)**2)/(2.d0*dist_mag)

       dist_ci = dist_mag - dist_cl

       tow_force(:) = des_crossprdct(normal(:), ft(:))
       tow_tmp(:,1) = dist_cl*tow_force(:)
       tow_tmp(:,2) = dist_ci*tow_force(:)

       ! calculate the total force fc of a collision pair
       ! total contact force
       fc_tmp(:) = fn(:) + ft(:)

       fc(ll,:) = fc(ll,:) + fc_tmp(:)

       fc(i,1) = fc(i,1) - fc_tmp(1)
       fc(i,2) = fc(i,2) - fc_tmp(2)
       fc(i,3) = fc(i,3) - fc_tmp(3)

       ! for each particle the signs of norm and ft both flip, so add the same torque
       tow(ll,:) = tow(ll,:) + tow_tmp(:,1)

       tow(i,1)  = tow(i,1)  + tow_tmp(1,2)
       tow(i,2)  = tow(i,2)  + tow_tmp(2,2)
       tow(i,3)  = tow(i,3)  + tow_tmp(3,2)

    enddo

    return

  contains

    include 'functions.inc'

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !                                                                      !
    !  subroutine: print_excess_overlap                                    !
    !                                                                      !
    !  purpose: print overlap warning messages.                            !
    !                                                                      !
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
    subroutine print_excess_overlap

      if(overlap_n > flag_overlap*radius(ll) .or.                  &
           overlap_n > flag_overlap*radius(i)) then

         write(err_msg,1000) trim(ival(ll)), trim(ival(i)), s_time,    &
              radius(ll), radius(i), overlap_n

         call flush_err_msg(header=.false., footer=.false.)
      endif

1000  format('warning: excessive overplay detected between ',          &
           'particles ',a,' and ',/a,' at time ',g11.4,'.',/             &
           'radii:  ',g11.4,' and ',g11.4,4x,'overlap: ',g11.4)

    end subroutine print_excess_overlap

  end subroutine calc_force_dem



  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  ! Subroutine: calc_collisions                                          !
  ! Purpose: Build a list of collision pairs.                            !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  SUBROUTINE CALC_COLLISIONS(np, pairs, pair_count, &
       state, radius, pos)

    use discretelement, only: nonexistent
    use param, only: small_number

    implicit none

    ! Dummy arguments ......................................................
    integer, intent(in) :: np

    integer     , intent(  out) :: pair_count
    integer     , intent(  out) :: pairs(6*np,2)

    real(c_real), intent(in   ) :: pos(np,3)
    real(c_real), intent(in   ) :: radius(np)
    integer     , intent(in   ) :: state(np)

    ! Local variables ......................................................
    integer :: i, ll
    real(c_real) :: rad
    real(c_real) :: DIST(3), DIST_MAG, pos_tmp(3)

    pair_count = 0

    DO LL = 1, np-1

       IF(NONEXISTENT==STATE(LL)) CYCLE
       pos_tmp = POS(LL,:)
       rad = RADIUS(LL)

       DO I = LL+1, np
          IF(NONEXISTENT==STATE(I)) CYCLE

          DIST(:) = POS(I,:) - POS_tmp(:)
          DIST_MAG = dot_product(DIST,DIST)

          IF(DIST_MAG < (rad + RADIUS(I) - SMALL_NUMBER)**2) THEN
             pair_count = pair_count + 1
             pairs(pair_count, 1) = ll
             pairs(pair_count, 2) = i
          ENDIF
       ENDDO
    ENDDO
  end subroutine calc_collisions


end module calc_force_dem_module
