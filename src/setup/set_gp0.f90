!                                                                      !
!  Subroutine: set_p0                                                  !
!                                                                      !
!  Purpose: Set the pressure field inside the bed assuming gravity     !
!           is acting in the negative y-direction.                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine set_gp0(domlo, domhi, gp0, &
                   dx, dy, dz, xlength, ylength, zlength, delp_dir_in) &
                   bind(C, name="set_gp0")

   use bc       , only: delp_x, delp_y, delp_z
   use bc       , only: dim_bc, bc_type, bc_p_g, bc_defined
   use bc       , only: pinf_, pout_, minf_
   use constant , only: gravity
   use fld_const, only: ro_g0
   use ic       , only: ic_p_g, ic_defined
   use scales   , only: scale_pressure

   use amrex_fort_module, only : ar => amrex_real
   use iso_c_binding , only: c_int
   use param   , only: zero, undefined
   use param   , only: is_defined, is_undefined
   use param, only: dim_ic

   implicit none

   integer, intent(in) :: domlo(3), domhi(3)

   real(ar), intent(inout) :: gp0(3)

   real(ar), intent(in) :: dx, dy, dz
   real(ar), intent(in) :: xlength, ylength, zlength
   integer , intent(in) :: delp_dir_in

   real(ar) :: offset = - 0.5_ar

   !-----------------------------------------------
   ! Local variables
   !-----------------------------------------------
   ! indices
   integer :: i, j, k, ibc, jbc, kbc
   integer :: icv, bcv, bcv_lo, bcv_hi
   integer :: delp_dir

   ! Gas pressure at the axial location j
   real(ar) :: pj, p_lo, p_hi

   ! Average pressure drop per unit length
   real(ar) :: dpodx, dpody, dpodz

   delp_dir = delp_dir_in

   !  Make sure that ic_p_g is set if using delp pressure conditions
   do icv = 1, dim_ic
      if (ic_defined(icv)) then
         if ( (delp_dir .ge. 0) .and. (delp_dir .eq. delp_dir_in) ) then
            if (.not. is_defined(ic_p_g(icv))) then
               print *,'MUST DEFINE ic_p_g if using the DELP pressure condition'
               stop
            end if

         else if ( (delp_dir .ge. 0) .and. (delp_dir .ne. delp_dir_in) ) then

            if (is_defined(ic_p_g(icv))) then
               print *,'MUST not define ic_p_g if setting p_inflow and p_outflow'
               stop
            end if

         else
            if ( (is_undefined(ic_p_g(icv))) .or. &
                 (gravity(1).ne.0.d0 .or. gravity(2).ne.0.d0 .or. gravity(3).ne.0.d0) ) goto 60
         end if
      end if
   end do

   ! ---------------------------------------------------------------->>>

   ! Here the pressure in each cell is determined from a specified pressure
   ! drop across the domain length. This section requires that the pressure
   ! is already defined in all initial condition regions (otherwise this
   ! section would be skipped)

   !  This hack allows to set the IC pressure  at L-dx/2 for both
   !  nodal and CC pressure -> reference value for pressure, AKA IC_P_G,
   !  is set at the last cell center location.
   if (delp_dir .ne. delp_dir_in) offset = -1.0_ar

   if (abs(delp_x) > epsilon(zero)) &
      gp0(1) = -delp_x/xlength

   if (abs(delp_y) > epsilon(zero)) &
      gp0(2) = -delp_y/ylength

   if (abs(delp_z) > epsilon(zero)) &
      gp0(3) = -delp_z/zlength

   goto 100   ! pressure in all initial condition region cells was defined

   ! ----------------------------------------------------------------<<<

60 continue   ! pressure in an initial condition region cell was undefined

   ! ----------------------------------------------------------------<<<

   ! Set an approximate pressure field assuming that the pressure drop
   ! balances the weight of the bed, if the initial pressure-field is not
   ! specified

   if (abs(gravity(1)) > epsilon(0.0d0)) then

      gp0(1) = ro_g0 * gravity(1)

   else if (abs(gravity(2)) > epsilon(0.0d0)) then

      gp0(2) = ro_g0 * gravity(2)

   else if (abs(gravity(3)) > epsilon(0.0d0)) then

      gp0(3) = ro_g0 * gravity(3)

   endif

   ! ----------------------------------------------------------------<<<

100 continue

end subroutine set_gp0
