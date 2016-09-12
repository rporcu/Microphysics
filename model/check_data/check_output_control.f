!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CHECK_OUTPUT_CONTROL                                    !
!  Purpose: Check the output control namelist section                  !
!                                                                      !
!  Author: P. Nicoletti                               Date: 27-NOV-91  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_OUTPUT_CONTROL

! Global Variables:
!---------------------------------------------------------------------//
! Time intervalue between updating the RES and SPx files.
      use output, only: RES_DT, SPX_DT
! Time-step intervalue between updating the .LOG file.
      use output, only: NLOG
! VTK
      use vtk
      USE run, only: RUN_NAME
      USE physprop, only: MMAX
      USE mpi_utility, only: XLENGTH,YLENGTH,ZLENGTH
      USE DISCRETELEMENT, only:DISCRETE_ELEMENT
      USE DISCRETELEMENT, only: PARTICLE_ORIENTATION
      USE cutcell, only: USE_STL

! Global Parameters:
!---------------------------------------------------------------------//
! Number aliases
      use param1, only: UNDEFINED, UNDEFINED_I, ZERO, LARGE_NUMBER
! Number of SPx files.
      USE param1, only: N_SPX

! Global Module procedures:
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter
      INTEGER :: LC

      INTEGER :: L,M,N,LV,N_VTK_REGIONS,R

!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_OUTPUT_CONTROL")


! Check the values specified for the RES file.
      IF (RES_DT==UNDEFINED)THEN
         WRITE(ERR_MSG,1000) 'RES_DT'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(RES_DT <= ZERO) THEN
         WRITE(ERR_MSG,1002) 'RES_DT', RES_DT
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF


! Check the SPx Files
      SPx_LP: DO LC = 1, N_SPX

! Disable writing the .SPB file if K-Epsilon is unspecified.
         IF(LC == 11) THEN
            IF (SPX_DT(LC)==UNDEFINED)SPX_DT(LC) = LARGE_NUMBER
            CYCLE SPx_LP

! Verify the remaining SPx files.
         ELSE
            IF(SPX_DT(LC) == UNDEFINED) THEN
               WRITE(ERR_MSG,1000) iVar('SPX_DT',LC)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

            ELSEIF(SPX_DT(LC) <= ZERO) THEN
               WRITE(ERR_MSG,1001) iVar('SPX_DT',LC), SPX_DT(LC)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF
      ENDDO SPx_LP

! Verify that the LOG frequency is valid.
      IF(NLOG <= 0) THEN
         WRITE(ERR_MSG,1003) 'NLOG', NLOG
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Check VTK regions

      IF(FRAME(1)<-1) THEN
         WRITE(ERR_MSG, 2000) trim(iVAL(FRAME(1)))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 2000 FORMAT('Error 2000: Invalid value for FRAME = ',A,'. Acceptable ',&
         'values',/'are integers >= -1. Please correct mfix.dat and',/ &
         'try again.')

      IF(VTK_DT(1)<ZERO) THEN
          WRITE(ERR_MSG,2001) trim(iVal(VTK_DT(1)))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 2001 FORMAT('Error 2001: Invalid value for VTK_DT = ',A,'. Acceptable',&
         ' values',/'are positive numbers (e.g., 0.1).  Please ',      &
         'correct mfix.dat and',/'try again.')

      N_VTK_REGIONS = 0
      DO L = 1, DIMENSION_VTK
         VTK_DEFINED(L) = .FALSE.
         IF (VTK_X_W(L) /= -UNDEFINED)   VTK_DEFINED(L) = .TRUE.
         IF (VTK_X_E(L) /=  UNDEFINED)   VTK_DEFINED(L) = .TRUE.
         IF (VTK_Y_S(L) /= -UNDEFINED)   VTK_DEFINED(L) = .TRUE.
         IF (VTK_Y_N(L) /=  UNDEFINED)   VTK_DEFINED(L) = .TRUE.
         IF (VTK_Z_B(L) /= -UNDEFINED)   VTK_DEFINED(L) = .TRUE.
         IF (VTK_Z_T(L) /=  UNDEFINED)   VTK_DEFINED(L) = .TRUE.

         IF(.NOT.VTK_DEFINED(L)) CYCLE
         N_VTK_REGIONS =  N_VTK_REGIONS + 1
      ENDDO   ! end loop over (l = 1,dimension_vtk)

! There must be at least one VTK region defined
! If this is not the case, define the entire domain as default region
      IF(WRITE_VTK_FILES.AND.N_VTK_REGIONS==0) THEN
         VTK_DEFINED(1) = .TRUE.
         VTK_X_W(1) = ZERO
         VTK_X_E(1) = XLENGTH
         VTK_Y_S(1) = ZERO
         VTK_Y_N(1) = YLENGTH
         VTK_Z_B(1) = ZERO
         VTK_Z_T(1) = ZLENGTH
         VTK_FILEBASE(1) = RUN_NAME
      ENDIF

! If VTK_VAR is defined, fill-up the variable list
! for the vtk subdomains
      DO L = 1, DIM_VTK_VAR
         IF(VTK_VAR(L)/=UNDEFINED_I) VTK_VARLIST(:,L) = VTK_VAR(L)
      ENDDO


      DO L = 1, DIMENSION_VTK

         IF(.NOT.VTK_DEFINED(L)) CYCLE

         DO LV = 1, DIM_VTK_VAR

            SELECT CASE (VTK_VARLIST(L,LV))

               CASE (1)
                  VTK_EP_g(L) = .TRUE.

               CASE (2)
                  VTK_P_g(L)    = .TRUE.

               CASE (3)
                  VTK_VEL_G(L) = .TRUE.

               CASE (4)
                  DO M = 1,MMAX
                     VTK_VEL_S(L,M) = .TRUE.
                  END DO

               CASE (5)
                  DO M = 1,MMAX
                     VTK_ROP_s(L,M) = .TRUE.
                  END DO


               CASE (100)
                  VTK_PARTITION(L) = .TRUE.

               CASE (101)
                  VTK_BC_ID(L) = .TRUE.

               CASE (102)
                  VTK_DWALL(L) = .TRUE.

               CASE (103)
                  IF(DISCRETE_ELEMENT.AND.USE_STL) THEN
                     VTK_FACET_COUNT_DES(L) = .TRUE.
                  ENDIF

               CASE (104)
                  IF(DISCRETE_ELEMENT.AND.USE_STL) THEN
                     VTK_NB_FACET_DES(L) = .TRUE.
                  ENDIF

               CASE(999)
                  VTK_IJK(L) = .TRUE.

               CASE(1000)
                  VTK_NORMAL(L) = .TRUE.

               CASE (1001)
                  VTK_DEBUG(L,1) = .TRUE.

               CASE (1002)
                  VTK_DEBUG(L,2) = .TRUE.

               CASE (1003)
                  VTK_DEBUG(L,3) = .TRUE.

               CASE (1004)
                  VTK_DEBUG(L,4) = .TRUE.

               CASE (1005)
                  VTK_DEBUG(L,5) = .TRUE.

               CASE (1006)
                  VTK_DEBUG(L,6) = .TRUE.

               CASE (1007)
                  VTK_DEBUG(L,7) = .TRUE.

               CASE (1008)
                  VTK_DEBUG(L,8) = .TRUE.

               CASE (1009)
                  VTK_DEBUG(L,9) = .TRUE.

               CASE (1010)
                  VTK_DEBUG(L,10) = .TRUE.

               CASE (1011)
                  VTK_DEBUG(L,11) = .TRUE.

               CASE (1012)
                  VTK_DEBUG(L,12) = .TRUE.

               CASE (1013)
                  VTK_DEBUG(L,13) = .TRUE.

               CASE (1014)
                  VTK_DEBUG(L,14) = .TRUE.

               CASE (1015)
                  VTK_DEBUG(L,15) = .TRUE.


               CASE (0) ! do nothing

               CASE (UNDEFINED_I) ! do nothing

               CASE DEFAULT
                  WRITE(ERR_MSG,2100) trim(iVal(L)),                   &
                     trim(iVal(VTK_VAR(L)))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            END SELECT

 2100 FORMAT(' Error 2100: Unknown VTK variable flag ',A,':',A,       /&
         'Available flags are:',                                      /&
         '  1 : Void fraction (EP_g)',                                /&
         '  3 : Gas velocity (U_g, V_g, W_g)',                        /&
         'Please correct the mfix.dat file.')


         ENDDO

! Activate particle orientation calculation if one vtk region needs it.
         IF(VTK_PART_ORIENTATION(L)) PARTICLE_ORIENTATION = .TRUE.

      ENDDO   ! end loop over (l = 1,dimension_vtk)

! Finalize the error manager.
      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
         'Please correct the mfix.dat file.')

 1002 FORMAT('Error 1002: Illegal or unknown input: ',A,' = ',E14.6,/  &
         'Please correct the mfix.dat file.')

 1003 FORMAT('Error 1003: Illegal or unknown input: ',A,' = ',I4,/     &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_OUTPUT_CONTROL
