!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: set_ps                                                 C
!  Purpose: Check point source specifications.                         C
!  Note: This is called on all MPI processses but we have a flag for   C
!        whether this is the I/O processor                             C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   subroutine set_ps(dx,dy,dz,err,is_ioproc) &
      bind(C, name="set_ps")
 
      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      use ps, only: dimension_ps, point_source, ps_vel_mag_g, ps_massflow_g, ps_vel_mag_g, ps_volume, ps_defined
      use ps, only: ps_i_e, ps_i_w, ps_j_n, ps_j_s, ps_k_t, ps_k_b
      use ps, only: ps_x_w, ps_x_e, ps_y_s, ps_y_n, ps_z_b, ps_z_t
      use ps, only: ps_u_g, ps_v_g, ps_w_g

      use param1  , only: zero, small_number

      implicit none

      real(c_real)  , intent(in   ) :: dx,dy,dz

      ! 1 if on the I/O processor, 0 otherwise
      integer(c_int), intent(in   ) :: is_ioproc

      ! If err = 1  is passed back to the C++ then amrex::Abort is called.
      integer(c_int), intent(  out) :: err

      integer(c_int) :: psv, ps_size
      real(c_real)   :: vol

      CHARACTER(LEN=64) :: eMsg

      logical, parameter :: dbg_ps = .FALSE.

      err = 0

      vol = dx*dy*dz

      if(.not.point_source) return

      ! DETERMINE WHICH BOUNDARY CONDITION INDICES HAVE VALUES
      L50: do PSV = 1, dimension_ps

         IF(.NOT.PS_DEFINED(PSV)) cycle L50

         ! Calculate the velocity magnitude and normalize the axial components.
         call calc_ps_vel_mag(PS_VEL_MAG_g(PSV), PS_U_g(PSV),          &
            PS_V_g(PSV), PS_W_g(PSV))

         ! Calculate the number of cells comprising the point source. 
         ps_size = (PS_I_E(PSV) - PS_I_W(PSV) + 1) * &
                   (PS_J_N(PSV) - PS_J_S(PSV) + 1) * &
                   (PS_K_T(PSV) - PS_K_B(PSV) + 1)

         if(ps_size < 1) then
             eMsg = ''; 
             if (is_ioproc .eq. 1) &
                write(eMsg,"('Invalid PS size: ', I4)")ps_size
             goto 500
         endif

         ! Calculate the volume of the PointSource cells
         ps_volume(PSV) = ps_size * vol

         if(abs(ps_volume(PSV)) < epsilon(ZERO)) then
            eMsg = 'No ps_volume == ZERO'
            if (is_ioproc .eq. 1) &
               call debug_ps(PSV, ps_size)
            goto 501
         endif

         if ( dbg_ps .and. (is_ioproc .eq. 1) ) call debug_ps(PSV, ps_size)

      enddo L50

      return

  500 continue

      if (is_ioproc .eq. 1) &
         write(*,"('PointSource setup Error: ',A)") trim(eMsg)
      err = 1

  501 continue

      if (is_ioproc .eq. 1) &
         write(*,"('PointSource setup Error: ',A)") trim(eMsg)
      err = 1

      contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_PS_VEL_MAG                                        C
!  Purpose: Calculate velocity magnitude and normalize U/V/W.          C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine calc_ps_vel_mag(VEL_MAG, lU, lV, lW)

      real(c_real), intent(inout) :: VEL_MAG
      real(c_real), intent(inout) :: lU, lV, lW

! Normalize velocities:
      VEL_MAG = lU**2 + lV**2 + lW**2

      VEL_MAG = sqrt(VEL_MAG)

      if(VEL_MAG > small_number) then
         lU = lU/VEL_MAG
         lV = lV/VEL_MAG
         lW = lW/VEL_MAG
      else
         VEL_MAG = ZERO
         lU = ZERO
         lV = ZERO
         lW = ZERO
      endif

      end subroutine calc_ps_vel_mag

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: debug_ps                                               C
!                                                                      C
!  Purpose: Write out some information about that point source that    C
!  may be useful in debugging.                                         C
!                                                                      C
!  Author: J. Musser                                  Date: 24-JUN-13  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine debug_ps(lPSV, lps_size)

      use param1  , only: small_number

      implicit none

      ! Index of PS source to debug.
      integer, intent(in) :: lPSV

      ! Number of cells comprising the point source.
      integer, intent(in) :: lps_size

      integer :: lc1

      write(*,"(3/,3x,'Debug Point Source Index: ',I3)") lPSV
      write(*,"(/3x,'Size: ',I4)") lps_size

      lc1 = 0

      ! Write some information to the screen.
      write(*,"(/5x,'Location:')")
      write(*,"( 5x,'X:',2(2x,g12.5),' :: ',2(2x,I4))")&
         PS_X_w(lPSV), PS_X_e(lPSV), PS_I_w(lPSV), PS_I_e(lPSV)
      write(*,"( 5x,'Y:',2(2x,g12.5),' :: ',2(2x,I4))")&
         PS_Y_s(lPSV), PS_Y_n(lPSV), PS_J_s(lPSV), PS_J_n(lPSV)
      write(*,"( 5x,'Z:',2(2x,g12.5),' :: ',2(2x,I4))")&
         PS_Z_b(lPSV), PS_Z_t(lPSV), PS_K_b(lPSV), PS_K_t(lPSV)

      write(*,"(/5x,'Volume: ',g12.5)") ps_volume(lPSV)

      if(PS_MASSFLOW_G(lPSV) > small_number) then
         write(*,"(//5x,'Point Source Gas Phase:')")
         write(*,"(7x,'Mass Flow Rate: ',g12.5)")PS_MASSFLOW_G(lPSV)
         write(*,"(7x,'Velocity Magnitude: ',g12.5)") PS_VEL_MAG_g(lPSV)
         write(*,"(7x,'Normal:')")
         write(*,"(9x,'x-Axis: ',g12.5)")PS_U_g(lPSV)
         write(*,"(9x,'y-Axis: ',g12.5)")PS_V_g(lPSV)
         write(*,"(9x,'z-Axis: ',g12.5)")PS_W_g(lPSV)
      else
         write(*,"(//5x,'No gas phase point source.')")
      endif

      end subroutine debug_ps

      end subroutine set_ps
