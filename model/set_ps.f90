MODULE SET_PS_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   use ps, only: dimension_ps, point_source, ps_vel_mag_g, ps_massflow_g, ps_vel_mag_g, ps_volume, ps_defined
   use ps, only: ps_i_e, ps_i_w, ps_j_n, ps_j_s, ps_k_t, ps_k_b
   use ps, only: ps_x_w, ps_x_e, ps_y_s, ps_y_n, ps_z_b, ps_z_t
   use ps, only: ps_u_g, ps_v_g, ps_w_g

   CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_09                                          C
!  Purpose: Check point source specifications.                         C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine set_ps(slo,shi,lo,hi,flag)

      use compar  , only: myPE, PE_IO, numPEs
      use exit_mod, only: mfix_exit
      use geometry, only: vol
      use param1  , only: zero, small_number

      implicit none

      integer(c_int), intent(in) :: slo(3), shi(3), lo(3), hi(3)
      integer(c_int), intent(in) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      integer :: i,j,k
      integer :: psv, ps_size

      CHARACTER(LEN=64) :: eMsg

      DOUBLE PRECISION, allocatable :: lData_dp(:)
      DOUBLE PRECISION, allocatable :: gData_dp(:)

      logical, parameter :: dbg_PS = .FALSE.

      if(.NOT.POINT_SOURCE) return

! DETERMINE WHICH BOUNDARY CONDITION INDICES HAVE VALUES
      L50: do PSV = 1, DIMENSION_PS

         IF(.NOT.PS_DEFINED(PSV)) cycle L50

! Calculate the velocity magnitude and normalize the axial components.
         CALL CALC_PS_VEL_MAG(PS_VEL_MAG_g(PSV), PS_U_g(PSV),          &
            PS_V_g(PSV), PS_W_g(PSV))


! Calculate the number of cells comprising the point source. This
! information is used to allocate some temp arrays.
!---------------------------------------------------------------------//
         PS_SIZE = (PS_I_E(PSV) - PS_I_W(PSV) + 1) * &
                   (PS_J_N(PSV) - PS_J_S(PSV) + 1) * &
                   (PS_K_T(PSV) - PS_K_B(PSV) + 1)

         if(PS_SIZE < 1) then
             eMsg = ''; write(eMsg,"('Invalid PS size: ', I4)")PS_SIZE
             goto 500
         endif

         allocate(lData_dp( 0:numPEs-1)); lData_dp = ZERO
         allocate(gData_dp( 0:numPEs-1)); gData_dp = ZERO


! Calculate the volume of the PointSource cells
!---------------------------------------------------------------------//
! Initialize the loop counter
         PS_VOLUME(PSV) = ZERO

         do k = PS_K_B(PSV), PS_K_T(PSV)
            do j = PS_J_S(PSV), PS_J_N(PSV)
               do i = PS_I_W(PSV), PS_I_E(PSV)
                  if(1.eq.flag(i,j,k,1)) &
                     lData_dp(myPE) = lData_dp(myPE) + VOL
               enddo
            enddo
         enddo

! Each process (in DMP) only knows about the volume of the point source
! it sees. Invoke send/recv so that all process can calculate the total
! volume.
         ! CALL global_all_sum(lData_dp, gData_dp)

         PS_VOLUME(PSV) = sum(gData_dp)
         if(abs(PS_VOLUME(PSV)) < epsilon(ZERO)) then
            eMsg = 'No PS_VOLUME == ZERO'
            CALL DEBUG_PS(PSV, PS_SIZE, flag)
            goto 501
         endif

         if(allocated(lData_dp)) deallocate(lData_dp)
         if(allocated(gData_dp)) deallocate(gData_dp)


         IF(dbg_PS) CALL DEBUG_PS(PSV, PS_SIZE, flag)


      enddo L50

      return

  500 continue
      if(myPE == PE_IO) &
         write(*,"('PointSource setup Error: ',A)") trim(eMsg)

      call mfix_exit(myPE)


  501 continue
      write(*,"('PointSource setup Error: ',I3,2x,A)") myPE, trim(eMsg)

      call mfix_exit(myPE)


      RETURN

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

      DOUBLE PRECISION, intent(inout) :: VEL_MAG
      DOUBLE PRECISION, intent(inout) :: lU, lV, lW

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

      end subroutine set_ps

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DEBUG_PS                                               C
!                                                                      C
!  Purpose: Write out some information about that point source that    C
!  may be useful in debugging.                                         C
!                                                                      C
!  Author: J. Musser                                  Date: 24-JUN-13  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine debug_ps(lPSV, lPS_SIZE, flag)

      use compar  , only: myPE, PE_IO
      use param1  , only: small_number

      implicit none

! Index of PS source to debug.
      INTEGER, intent(in) :: lPSV
! Number of cells comprising the point source.
      INTEGER, intent(in) :: lPS_SIZE

      INTEGER, DIMENSION(:,:,:,:), INTENT(IN) :: FLAG

      INTEGER :: IJK, I, J, K

      INTEGER :: lc1

      INTEGER, allocatable :: lFlags_i(:,:)
      INTEGER, allocatable :: gFlags_i(:,:)

      if(myPE == PE_IO) then
         write(*,"(3/,3x,'Debug Point Source Index: ',I3)") lPSV
         write(*,"(/3x,'Size: ',I4)") lPS_SIZE
      endif

      allocate(lFlags_i(lPS_SIZE,1:2) );   lFlags_i = 0
      allocate(gFlags_i(lPS_SIZE,1:2) );   gFlags_i = 0

      lc1 = 0

      do k = PS_K_B(lPSV), PS_K_T(lPSV)
         do j = PS_J_S(lPSV), PS_J_N(lPSV)
            do i = PS_I_W(lPSV), PS_I_E(lPSV)

               lc1 = lc1 + 1
               if(1.eq.flag(i,j,k,1)) then
                  lFlags_i(lc1,1) = myPE
                  lFlags_i(lc1,2) = FLAG(i,j,k,1)
               endif
            enddo
         enddo
      enddo

! Collect flag information on root.
      ! CALL global_sum(lFlags_i, gFlags_i)

! Write some information to the screen.
      if(myPE == PE_IO) then
         write(*,"(/5x,'Location:')")
         write(*,"( 5x,'X:',2(2x,g12.5),' :: ',2(2x,I4))")&
            PS_X_w(lPSV), PS_X_e(lPSV), PS_I_w(lPSV), PS_I_e(lPSV)
         write(*,"( 5x,'Y:',2(2x,g12.5),' :: ',2(2x,I4))")&
            PS_Y_s(lPSV), PS_Y_n(lPSV), PS_J_s(lPSV), PS_J_n(lPSV)
         write(*,"( 5x,'Z:',2(2x,g12.5),' :: ',2(2x,I4))")&
            PS_Z_b(lPSV), PS_Z_t(lPSV), PS_K_b(lPSV), PS_K_t(lPSV)

         write(*,"(/5x,'Volume: ',g12.5)") PS_VOLUME(lPSV)


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

         write(*,"(//5x,'Point Source Cells:')")
         write(*,"(9x,'IJK',3(6x,A1),3x,'OWNS',3x,'FLAG')") 'I','J','K'

         lc1 = 0
         do k = PS_K_B(lPSV), PS_K_T(lPSV)
         do j = PS_J_S(lPSV), PS_J_N(lPSV)
         do i = PS_I_W(lPSV), PS_I_E(lPSV)
            lc1 = lc1 + 1
            write(*,"(4x,I8,5(3x,I4))") IJK, I, J, K,  gFlags_i(lc1,:)
         enddo
         enddo
         enddo

      endif

      if(allocated(lFlags_i)) deallocate(lFlags_i)
      if(allocated(gFlags_i)) deallocate(gFlags_i)

      end subroutine debug_ps
end module set_ps_module
