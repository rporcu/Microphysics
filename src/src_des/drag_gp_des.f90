!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: des_drag_gp                                             C
!  Purpose: Calculate the gas-particle drag coefficient using          C
!           the gas velocity interpolated to the particle position     C
!           and the particle velocity.                                 C
!           Invoked from des_drag_gs and calc_des_drag_gs              C
!                                                                      C
!  Comments: The BVK2 drag model and all drag models with the          C
!            polydisperse correction factor (i.e., suffix _PCF)        C
!            require an average particle diameter. This has been       C
!            loosely defined for discrete particles based on their     C
!            solids phase                                              C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine des_drag_gp(p_id, particle_vel, fluid_vel, &
                             ep_g, ro_g, mu_g, f_gp, &
                             i,j,k,des_radius,pvol,ro_s) &
                             bind(C,name="des_drag_gp")

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

      use drag,  only: drag_type_enum, drag_type
      use drag,  only: drag_wen_yu, drag_gidaspow, drag_bvk2
      use drag,  only: user_drag, wen_yu, gidaspow, bvk2, & 
                       wen_yu_pcf, gidaspow_pcf
      use param, only: one

      IMPLICIT NONE

      interface
         subroutine drag_usr(i,j,k,p_id,lDgA,ep_g,mu_g,ro_g, VREL, DPM, &
                             ro_s, lUg, lVg, lWg)
            use amrex_fort_module, only : rt => amrex_real
            use iso_c_binding , only: c_int

            integer(c_int), intent(in   ) :: i,j,k,p_id
            real(rt)      , intent(in   ) :: ep_g, mu_g, ro_g, VREL, DPM, ro_s, lUg, lVg, lWg
            real(rt)      , intent(  out) :: lDgA
         end subroutine drag_usr
      end interface

      ! Grid indices associated with current particle
      integer, intent(in   ) :: i,j,k

      ! Particle ID
      integer, intent(in   ) :: p_id

      real(rt), intent(in) :: des_radius
      real(rt), intent(in) :: pvol, ro_s

!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      ! particle velocity
      real(rt), intent(in) :: particle_vel(3)

      ! fluid velocity interpolated to particle position
      real(rt), intent(in) :: fluid_vel(3)

      ! Gas phase volume fraction.
      real(rt), intent(in) :: ep_g
      real(rt), intent(in) :: ro_g
      real(rt), intent(in) :: mu_g

      real(rt), intent(out) :: f_gp

!-----------------------------------------------
! Local variables
!-----------------------------------------------

      ! Slip velocity and its magnitude
      real(rt) :: VSLP(3), VREL

      ! drag coefficient
      real(rt) :: DgA

      ! total solids volume fraction
      real(rt) :: phi_s

      real(rt) :: rop_g, DPM
!-----------------------------------------------

      ! Gas material and bulk densities
      rop_g = ro_g * ep_g

      ! Slip velocity and its magnitude
      VSLP = fluid_vel - particle_vel
      VREL = SQRT(dot_product(VSLP, VSLP))

      DPM = 2.0d0*des_radius

! Total solids volume fraction.
      phi_s = ONE - ep_g

! determine the drag coefficient
      SELECT CASE(DRAG_TYPE_ENUM)

      CASE (USER_DRAG)
         CALL DRAG_USR(I,J,K,p_id,DgA,ep_g,mu_g,ro_g,VREL,DPM,ro_s, &
                       fluid_vel(1), fluid_vel(2), fluid_vel(3))

      CASE (WEN_YU)
         CALL DRAG_WEN_YU(DgA,ep_g,mu_g,rop_g,VREL,DPM)

      CASE (GIDASPOW)
         CALL DRAG_GIDASPOW(DgA,ep_g,mu_g,ro_g,rop_g,VREL,DPM)

      CASE (BVK2)
         ! calculate the average particle diameter and particle ratio
         ! HACK HACK HACK HACK -- Dependence on rop_s was removed

         CALL DRAG_BVK2(DgA,ep_g,mu_g,rop_g,VREL,DPM,DPM,phi_s)

      CASE DEFAULT
         WRITE (*, '(A,A)') 'Unknown DRAG_TYPE: ', DRAG_TYPE
         stop 20013
      end select   ! end selection of drag_type

! Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
      f_gp = DgA * pvol

      end subroutine des_drag_gp
