module calc_force_dem_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   implicit none
   private

   public calc_force_dem, calc_force_dem_nl


contains

   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
   !                                                                      !
   !  subroutine: calc_force_dem                                      !
   !  author: jay boyalakuntla                           date: 12-jun-04  !
   !                                                                      !
   !  purpose: calculate contact force and torque on particle from        !
   !           particle-particle and particle-wall collisions. treats     !
   !           wall interaction also as a two-particle interaction but    !
   !           accounting for the wall properties                         !
   !                                                                      !
   !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine calc_force_dem( particles, fc, tow, dtsolid )

      use particle_mod,   only: particle_t
      use cfrelvel_module, only: cfrelvel
      use discretelement, only: des_coll_model_enum
      use discretelement, only: des_etan, des_etat, hert_kt, hert_kn
      use discretelement, only: des_crossprdct
      use discretelement, only: kn, kt, mew, hertzian
      use param,          only: SMALL_NUMBER
      use error_manager, only: init_err_msg, flush_err_msg, err_msg, ival

      type(particle_t), intent(in   ) :: particles(:)
      real(c_real),     intent(inout) :: fc(:,:), tow(:,:)
      real(c_real),     intent(in   ) :: dtsolid
      logical,      parameter     :: report_excess_overlap = .false.
      real(c_real), parameter     :: flag_overlap = 0.20d0 ! % of particle radius when excess overlap will be flagged
      real(c_real), parameter     :: q2           = 0.5_c_real

      ! particle no. indices
      integer :: ii, ll
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
      integer :: phaseii, phasell
      ! local values used spring constants and damping coefficients
      real(c_real) :: etan_des, etat_des, kn_des, kt_des
      real(c_real) :: fnmd, mag_overlap_t, tangent(3)
      integer      ::  np
      real(c_real) :: radiusii, radiusll

      np = size(particles)

!      print *,'NP IN CALC_FORCE_DEM ',np

      do ll = 1, np-1

         pos_tmp = particles(ll) % pos
         rad     = particles(ll) % radius

         do ii = ll+1, np

            dist     = particles(ii) % pos - pos_tmp(:)
            dist_mag = dot_product( dist, dist )
            r_lm     = rad + particles(ii) % radius

            if ( dist_mag > ( r_lm - SMALL_NUMBER )**2 ) cycle

            if (abs(dist_mag) < epsilon(dist_mag)) then
               write(*,8550) ll, ii
               stop "division by zero"
8550           format('distance between particles is zero:',2(2x,i10))
            endif

            dist_mag  = sqrt( dist_mag )
            normal(:) = dist(:) / dist_mag

            ! calcuate the normal overlap
            overlap_n = r_lm-dist_mag
            if(report_excess_overlap) call print_excess_overlap

            ! calculate the components of translational relative velocity for a
            ! contacting particle pair and the tangent to the plane of contact
            call cfrelvel(ll, ii, v_rel_trans_norm, vrel_t, normal(:), dist_mag, particles )

            radiusll = particles(ll) % radius
            phasell  = particles(ll) % phase
            radiusii = particles(ii) % radius
            phaseii  = particles(ii) % phase

            ! hertz spring-dashpot contact model
            if ( des_coll_model_enum == hertzian ) then
               sqrt_overlap = sqrt(overlap_n)
               kn_des       = hert_kn(phasell,phaseii)*sqrt_overlap
               kt_des       = hert_kt(phasell,phaseii)*sqrt_overlap
               sqrt_overlap = sqrt(sqrt_overlap)
               etan_des     = des_etan(phasell,phaseii)*sqrt_overlap
               etat_des     = des_etat(phasell,phaseii)*sqrt_overlap
               ! linear spring-dashpot contact model
            else
               kn_des   = kn
               kt_des   = kt
               etan_des = des_etan(phasell,phaseii)
               etat_des = des_etat(phasell,phaseii)
            end if

            ! calculate the normal contact force
            fn(:) =  -(kn_des * overlap_n * normal(:) + &
                 etan_des * v_rel_trans_norm * normal(:))

            ! calcuate the tangential overlap
            overlap_t(:) = dtsolid*vrel_t(:)
            mag_overlap_t = sqrt(dot_product(overlap_t,overlap_t))

            ! calculate the tangential contact force.
            if (mag_overlap_t > 0.0) then
               ! max force before the on set of frictional slip.
               fnmd = mew*sqrt(dot_product(fn,fn))
               ! direction of tangential force.
               tangent = overlap_t/mag_overlap_t
               ! frictional slip
               ft = -fnmd * tangent
            else
               ft = 0.0
            end if

            ! calculate the distance from the particles' centers to the contact point,
            ! which is taken as the radical line
            ! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2
            dist_cl = q2 * ( dist_mag  + ( radiusll**2 - radiusii**2) / dist_mag )
            dist_ci = dist_mag - dist_cl

            tow_force(:) = des_crossprdct(normal(:), ft(:))
            tow_tmp(:,1) = dist_cl*tow_force(:)
            tow_tmp(:,2) = dist_ci*tow_force(:)

            ! calculate the total force fc of a collision pair
            ! total contact force
            fc_tmp(:) = fn(:) + ft(:)

            fc(ll,:) = fc(ll,:) + fc_tmp(:)

            fc(ii,1) = fc(ii,1) - fc_tmp(1)
            fc(ii,2) = fc(ii,2) - fc_tmp(2)
            fc(ii,3) = fc(ii,3) - fc_tmp(3)

            ! for each particle the signs of norm and ft both flip, so add the same torque
            tow(ll,:) = tow(ll,:) + tow_tmp(:,1)

            tow(ii,1)  = tow(ii,1)  + tow_tmp(1,2)
            tow(ii,2)  = tow(ii,2)  + tow_tmp(2,2)
            tow(ii,3)  = tow(ii,3)  + tow_tmp(3,2)

         end do
      end do

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

         if(overlap_n > flag_overlap*radiusll .or.                  &
              overlap_n > flag_overlap*radiusii) then

            write(err_msg,1000) trim(ival(ll)), trim(ival(ii)),     &
                 radiusll, radiusii, overlap_n

            call flush_err_msg(header=.false., footer=.false.)
         endif

1000     format('warning: excessive overplay detected between ',          &
              'particles ',a,' and ',/a,'.',/             &
              'radii:  ',g11.4,' and ',g11.4,4x,'overlap: ',g11.4)

      end subroutine print_excess_overlap

   end subroutine calc_force_dem

   subroutine calc_force_dem_nl( particles, nbor_list, size_nl, fc, tow, dtsolid )

      use particle_mod,   only: particle_t
      use cfrelvel_module, only: cfrelvel
      use discretelement, only: des_coll_model_enum
      use discretelement, only: des_etan, des_etat, hert_kt, hert_kn
      use discretelement, only: des_crossprdct
      use discretelement, only: kn, kt, mew, hertzian
      use param,          only: SMALL_NUMBER
      use error_manager, only: init_err_msg, flush_err_msg, err_msg, ival

      integer,          intent(in   ) :: size_nl
      type(particle_t), intent(in   ) :: particles(:)
      integer,          intent(in   ) :: nbor_list(size_nl)
      real(c_real),     intent(inout) :: fc(:,:), tow(:,:)
      real(c_real),     intent(in   ) :: dtsolid
      logical,      parameter     :: report_excess_overlap = .false.
      real(c_real), parameter     :: flag_overlap = 0.20d0 ! % of particle radius when excess overlap will be flagged
      real(c_real), parameter     :: q2           = 0.5_c_real

      ! particle no. indices
      integer :: ii, ll, jj
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
      integer :: phaseii, phasell
      ! local values used spring constants and damping coefficients
      real(c_real) :: etan_des, etat_des, kn_des, kt_des
      real(c_real) :: fnmd, mag_overlap_t, tangent(3)
      integer      ::  np, index, nneighbors
      real(c_real) :: radiusii, radiusll

      np = size(particles)

!      print *,'NP IN CALC_FORCE_DEM ',np

      index = 1
      do ll = 1, np

         pos_tmp = particles(ll) % pos
         rad     = particles(ll) % radius

         nneighbors = nbor_list(index)
         index = index + 1

         do jj = index, index + nneighbors - 1

            ii = nbor_list(jj)

            dist     = particles(ii) % pos - pos_tmp(:)
            dist_mag = dot_product( dist, dist )
            r_lm     = rad + particles(ii) % radius

            if ( dist_mag > ( r_lm - SMALL_NUMBER )**2 ) cycle

            if (abs(dist_mag) < epsilon(dist_mag)) then
               write(*,8550) ll, ii
               stop "division by zero"
8550           format('distance between particles is zero:',2(2x,i10))
            endif

            dist_mag  = sqrt( dist_mag )
            normal(:) = dist(:) / dist_mag

            ! calcuate the normal overlap
            overlap_n = r_lm-dist_mag
            if(report_excess_overlap) call print_excess_overlap

            ! calculate the components of translational relative velocity for a
            ! contacting particle pair and the tangent to the plane of contact
            call cfrelvel(ll, ii, v_rel_trans_norm, vrel_t, normal(:), dist_mag, particles )

            radiusll = particles(ll) % radius
            phasell  = particles(ll) % phase
            radiusii = particles(ii) % radius
            phaseii  = particles(ii) % phase

            ! hertz spring-dashpot contact model
            if ( des_coll_model_enum == hertzian ) then
               sqrt_overlap = sqrt(overlap_n)
               kn_des       = hert_kn(phasell,phaseii)*sqrt_overlap
               kt_des       = hert_kt(phasell,phaseii)*sqrt_overlap
               sqrt_overlap = sqrt(sqrt_overlap)
               etan_des     = des_etan(phasell,phaseii)*sqrt_overlap
               etat_des     = des_etat(phasell,phaseii)*sqrt_overlap
               ! linear spring-dashpot contact model
            else
               kn_des   = kn
               kt_des   = kt
               etan_des = des_etan(phasell,phaseii)
               etat_des = des_etat(phasell,phaseii)
            end if

            ! calculate the normal contact force
            fn(:) =  -(kn_des * overlap_n * normal(:) + &
                 etan_des * v_rel_trans_norm * normal(:))

            ! calcuate the tangential overlap
            overlap_t(:) = dtsolid*vrel_t(:)
            mag_overlap_t = sqrt(dot_product(overlap_t,overlap_t))

            ! calculate the tangential contact force.
            if (mag_overlap_t > 0.0) then
               ! max force before the on set of frictional slip.
               fnmd = mew*sqrt(dot_product(fn,fn))
               ! direction of tangential force.
               tangent = overlap_t/mag_overlap_t
               ! frictional slip
               ft = -fnmd * tangent
            else
               ft = 0.0
            end if

            ! calculate the distance from the particles' centers to the contact point,
            ! which is taken as the radical line
            ! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2
            dist_cl = q2 * ( dist_mag  + ( radiusll**2 - radiusii**2) / dist_mag )
            dist_ci = dist_mag - dist_cl

            tow_force(:) = des_crossprdct(normal(:), ft(:))
            tow_tmp(:,1) = dist_cl*tow_force(:)
            tow_tmp(:,2) = dist_ci*tow_force(:)

            ! calculate the total force fc of a collision pair
            ! total contact force
            fc_tmp(:) = fn(:) + ft(:)

            fc(ll,:) = fc(ll,:) + fc_tmp(:)

            fc(ii,1) = fc(ii,1) - fc_tmp(1)
            fc(ii,2) = fc(ii,2) - fc_tmp(2)
            fc(ii,3) = fc(ii,3) - fc_tmp(3)

            ! for each particle the signs of norm and ft both flip, so add the same torque
            tow(ll,:) = tow(ll,:) + tow_tmp(:,1)

            tow(ii,1)  = tow(ii,1)  + tow_tmp(1,2)
            tow(ii,2)  = tow(ii,2)  + tow_tmp(2,2)
            tow(ii,3)  = tow(ii,3)  + tow_tmp(3,2)

         end do

         index = index + nneighbors

      end do

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

         if(overlap_n > flag_overlap*radiusll .or.                  &
              overlap_n > flag_overlap*radiusii) then

            write(err_msg,1000) trim(ival(ll)), trim(ival(ii)),     &
                 radiusll, radiusii, overlap_n

            call flush_err_msg(header=.false., footer=.false.)
         endif

1000     format('warning: excessive overplay detected between ',          &
              'particles ',a,' and ',/a,'.',/             &
              'radii:  ',g11.4,' and ',g11.4,4x,'overlap: ',g11.4)

      end subroutine print_excess_overlap

   end subroutine calc_force_dem_nl

end module calc_force_dem_module
