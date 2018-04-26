!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: set_p0                                                 !
!  Author: M. Syamlal                                 Date: 21-JAN-92  !
!                                                                      !
!  Purpose: Set the pressure field inside the bed assuming gravity     !
!           is acting in the negative y-direction.                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine set_p0(slo, shi, lo, hi, domlo, domhi, &
                        p0_g, dx, dy, dz, xlength, ylength, zlength, delp_dir) &
                 bind(C, name="set_p0")

      use bc, only: delp_x, delp_y, delp_z
      use bc, only: dim_bc, bc_type, bc_p_g, bc_defined
      use constant , only: gravity
      use eos, ONLY: EOSG
      use fld_const, only: ro_g0
      use ic       , only: ic_p_g, ic_defined
      use scales   , only: scale_pressure

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int
      use param   , only: zero, undefined
      use param   , only: is_defined, is_undefined
      use param, only: dim_ic

      implicit none

      integer, intent(in) :: slo(3), shi(3), lo(3), hi(3)
      integer, intent(in) :: domlo(3), domhi(3)

      real(c_real), intent(inout) :: p0_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: dx, dy, dz
      real(c_real), intent(in   ) :: xlength, ylength, zlength
      integer     , intent(  out) :: delp_dir
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      integer :: i, j, k
! Local loop counter
      integer :: L
! Gas pressure at the axial location j
      real(c_real) :: pj
! Average pressure drop per unit length
      real(c_real) :: dpodx, dpody, dpodz
!-----------------------------------------------

      ! First pass out the direction in which we drop by delp (if any)
      ! so that we can set the correct periodicity flag in the C++
      if (abs(delp_x) > epsilon(zero)) then
         delp_dir = 0
      else if (abs(delp_y) > epsilon(zero)) then
         delp_dir = 1
      else if (abs(delp_z) > epsilon(zero)) then
         delp_dir = 2
      else 
         delp_dir = -1
      end if

!-----------------------------------------------

! If any initial pressures are unspecified skip next section
! calculations.
      do l = 1, dim_ic
         if (ic_defined(l)) then
            if (is_undefined(ic_p_g(l))) goto 60
            pj = ic_p_g(l)
         endif
      enddo

! Here the pressure in each cell is determined from a specified pressure
! drop across the domain length. This section requires that the pressure
! is already defined in all initial condition regions (otherwise this
! section would be skipped)
! ---------------------------------------------------------------->>>
      if (abs(delp_x) > epsilon(zero)) then
         dpodx = delp_x/xlength
         pj = pj - dpodx*dx*(hi(1)-domhi(1)+2)
         do i = hi(1)+1, lo(1)-1, -1
            pj = pj + dpodx*dx
            p0_g(i,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = scale_pressure(pj)
         enddo
      endif

      if (abs(delp_y) > epsilon(zero)) then
         dpody = delp_y/ylength
         pj = pj - dpody*dy*(hi(2)-domhi(2)+2)
         do j = hi(2)+1, lo(2)-1, -1
            pj = pj + dpody*dy
            p0_g(lo(1)-1:hi(1)+1,j,lo(3)-1:hi(3)+1) = scale_pressure(pj)
         enddo
      endif

      if (abs(delp_z) > epsilon(zero)) then
         dpodz = delp_z/zlength
         pj = pj - dpodz*dz*(hi(3)-domhi(3)+2)
         do k = hi(3)+1, lo(3)-1, -1
            pj = pj + dpodz*dz
            p0_g(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,k) = scale_pressure(pj)
         end do
      endif

! ----------------------------------------------------------------<<<

      GOTO 100   ! pressure in all intial condition region cells was defined
   60 CONTINUE   ! pressure in an initial condition region cell was undefined

! ---------------------------------------------------------------->>>
! Search for an outflow boundary condition where pressure is specified
      pj = undefined
      do l = 1, dim_bc
         if (bc_defined(l)) then
            if(bc_type(l)=='P_OUTFLOW' .or. bc_type(l)=='PO') &
               pj = bc_p_g(l)
         endif
      enddo

      ! Either a PO was not specified and/or a PO was specified but not the
      ! pressure at the outlet
      if (is_undefined(pj)) then
         p0_g = zero
         goto 100
      endif

! Set an approximate pressure field assuming that the pressure drop
! balances the weight of the bed, if the initial pressure-field is not
! specified

      if (abs(gravity(1)) > epsilon(0.0d0)) then

         ! Find the average weight per unit area over an x-z slice
         dpodx = -gravity(1)*ro_g0

         if (gravity(1) <= 0.0d0) then
            do i = domhi(1)+1, domlo(1), -1
               if (i <= hi(1)+1 .and. i >= lo(1)-1) &
                  p0_g(i,:,:) = scale_pressure(pj)
               pj = pj + dpodx*dx
            enddo
         else
            do i = domlo(1), domhi(1)+1
               if (i <= hi(1)+1 .and. i >= lo(1)-1) &
                  p0_g(i,:,:) = scale_pressure(pj)
               pj = pj - dpodx*dx
            enddo
         endif

      else if (abs(gravity(2)) > epsilon(0.0d0)) then

         dpody = -gravity(2)*ro_g0

         if (gravity(2) <= 0.0d0) then
            do j = domhi(2)+1, domlo(2), -1
               if (j <= hi(2)+1 .and. j>= lo(2)-1) &
                  p0_g(:,j,:) = scale_pressure(pj)
               pj = pj + dpody*dy
            enddo
         else
            do j = domlo(2),domhi(2)+1
               if (j <= hi(2)+1 .and. j >= lo(2)-1) &
                  p0_g(:,j,:) = scale_pressure(pj)
               pj = pj - dpody*dy
            enddo
         endif

      else if (abs(gravity(3)) > epsilon(0.0d0)) then

         dpodz = -gravity(3)*ro_g0

         if(gravity(3) <= 0.0d0) then
            do k = domhi(3)+1, domlo(3), -1
               if (k <= hi(3)+1 .and. k >= lo(3)-1) &
                  p0_g(:,:,k) = scale_pressure(pj)
               pj = pj + dpodz*dz
            enddo
         else
            do k = domlo(3),domhi(3)+1
               if (k <= hi(3)+1 .and. k >= lo(3)-1) &
                  p0_g(:,:,k) = scale_pressure(pj)
               pj = pj - dpodz*dz
            enddo
         endif
      endif

  100 continue

      return

 1000 FORMAT(/1X,70('*')//' From: SET_FLUIDBED_P'/' Message: Outflow ',&
         'pressure boundary condition (P_OUTFLOW) not found.',/&
         'All the initial pressures (IC_P_g) or at least one P_OUTFLOW',/&
         'condition need to be specified',/1X,70('*')/)

   end subroutine set_p0
