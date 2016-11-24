!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: WRITE_RES0                                              C
!  Author: P. Nicoletti                               Date: 13-DEC-91  C
!                                                                      C
!  Purpose: write out the initial restart records (namelist data)      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE WRITE_RES0

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc
      USE compar
      USE constant
      USE funits
      USE geometry
      USE ic
      USE in_binary_512i
      USE leqsol
      USE machine
      USE output
      USE param
      USE param1
      USE physprop
      USE run
      USE scales
      USE sendrecv         ! for filling the boundary information
      USE toleranc
      USE ur_facs
      USE fldvar

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer, allocatable :: arr1(:)
      integer, allocatable :: arr2(:)
! loop counters
      INTEGER :: LC, L
! Pointer to the next record
      INTEGER :: NEXT_RECA
! file version id

! Place holder for deprecated variables:
      LOGICAL, PARAMETER :: lCALL_ISAT = .FALSE.

      CHARACTER(LEN=512) :: VERSION
!-----------------------------------------------

      NEXT_RECA = 5

! note: the value below for version must be the same as the value
!       of version in mfix.f
!       If you change it in this subroutine, you must change it in
!       mfix.f also.

!       The value should also be consistent with the check in
!       read_res0

      VERSION = 'RES = 01.8'

! Add new data entries at the end of the file and identify version no.

      if (myPE.eq.PE_IO) then
         allocate (arr1(ijkmax3))
         allocate (arr2(ijkmax2))
      else
         allocate (arr1(1))
         allocate (arr2(1))
         goto 1100
      end if

!------------------------------------------------------------------------
      WRITE (UNIT_RES, REC=1) VERSION
      WRITE (UNIT_RES, REC=2) RUN_NAME, ID_MONTH, ID_DAY, ID_YEAR, ID_HOUR, &
         ID_MINUTE, ID_SECOND
      WRITE (UNIT_RES, REC=3) NEXT_RECA
      WRITE (UNIT_RES, REC=4) IMIN1, JMIN1, KMIN1, IMAX, JMAX, KMAX, IMAX1, &
         JMAX1, KMAX1, IMAX2, JMAX2, KMAX2, IJMAX2, IJKMAX2, MMAX, DIMENSION_IC&
         , DIMENSION_BC, DIMENSION_C, DT, XLENGTH, YLENGTH&
         , ZLENGTH
      CALL OUT_BIN_512 (UNIT_RES, C, DIMENSION_C, NEXT_RECA)
      NEXT_RECA = 1 + NEXT_RECA                  ! work around for -O3 compiler bug
      NEXT_RECA = NEXT_RECA - 1
      DO LC = 1, DIMENSION_C
         WRITE (UNIT_RES, REC=NEXT_RECA) C_NAME(LC)
         NEXT_RECA = NEXT_RECA + 1
      END DO

!      next_reca = next_reca + 1
!      CALL OUT_BIN_512I (UNIT_RES, work_around, MMAX+1, NEXT_RECA)
!      WRITE (UNIT_RES, REC=NEXT_RECA) (work_around(L),L=1,MMAX+1)



      NEXT_RECA = NEXT_RECA + 1
      write (*,*) ' next_reca = ' , next_reca

      CALL OUT_BIN_512 (UNIT_RES, DX(1), IMAX2, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, DY(1), JMAX2, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, DZ(1), KMAX2, NEXT_RECA)
      WRITE (UNIT_RES, REC=NEXT_RECA) RUN_NAME, DESCRIPTION, UNITS, RUN_TYPE, &
         COORDINATES
      NEXT_RECA = NEXT_RECA + 1
      WRITE (UNIT_RES, REC=NEXT_RECA) (D_P0(L),L=1,MMAX), (RO_S0(L),L=1,MMAX), &
         EP_STAR, RO_G0, MU_G0, MW_AVG
      NEXT_RECA = NEXT_RECA + 1


      CALL OUT_BIN_512 (UNIT_RES, IC_X_W, DIMENSION_IC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, IC_X_E, DIMENSION_IC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, IC_Y_S, DIMENSION_IC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, IC_Y_N, DIMENSION_IC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, IC_Z_B, DIMENSION_IC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, IC_Z_T, DIMENSION_IC, NEXT_RECA)
      CALL OUT_BIN_512I (UNIT_RES, IC_I_W, DIMENSION_IC, NEXT_RECA)
      CALL OUT_BIN_512I (UNIT_RES, IC_I_E, DIMENSION_IC, NEXT_RECA)
      CALL OUT_BIN_512I (UNIT_RES, IC_J_S, DIMENSION_IC, NEXT_RECA)
      CALL OUT_BIN_512I (UNIT_RES, IC_J_N, DIMENSION_IC, NEXT_RECA)
      CALL OUT_BIN_512I (UNIT_RES, IC_K_B, DIMENSION_IC, NEXT_RECA)
      CALL OUT_BIN_512I (UNIT_RES, IC_K_T, DIMENSION_IC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, IC_EP_G, DIMENSION_IC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, IC_P_G, DIMENSION_IC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, IC_T_G, DIMENSION_IC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, IC_U_G, DIMENSION_IC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, IC_V_G, DIMENSION_IC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, IC_W_G, DIMENSION_IC, NEXT_RECA)

      DO LC = 1, MMAX
         CALL OUT_BIN_512 (UNIT_RES, IC_ROP_S(1,LC), DIMENSION_IC, NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES, IC_U_S(1,LC), DIMENSION_IC, NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES, IC_V_S(1,LC), DIMENSION_IC, NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES, IC_W_S(1,LC), DIMENSION_IC, NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES, IC_T_S(1,LC), DIMENSION_IC, NEXT_RECA)
      END DO
      CALL OUT_BIN_512 (UNIT_RES, BC_X_W, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_X_E, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_Y_S, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_Y_N, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_Z_B, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_Z_T, DIMENSION_BC, NEXT_RECA)

      CALL OUT_BIN_512I (UNIT_RES, BC_I_W, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512I (UNIT_RES, BC_I_E, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512I (UNIT_RES, BC_J_S, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512I (UNIT_RES, BC_J_N, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512I (UNIT_RES, BC_K_B, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512I (UNIT_RES, BC_K_T, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_EP_G, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_P_G, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_T_G, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_U_G, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_V_G, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_W_G, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_RO_G, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_ROP_G, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_VOLFLOW_G, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_MASSFLOW_G, DIMENSION_BC, NEXT_RECA)
      DO LC = 1, MMAX
         CALL OUT_BIN_512 (UNIT_RES, BC_ROP_S(1,LC), DIMENSION_BC, NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES, BC_U_S(1,LC), DIMENSION_BC, NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES, BC_V_S(1,LC), DIMENSION_BC, NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES, BC_W_S(1,LC), DIMENSION_BC, NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES, BC_T_S(1,LC), DIMENSION_BC, NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES, BC_VOLFLOW_S(1,LC), DIMENSION_BC, &
            NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES, BC_MASSFLOW_S(1,LC), DIMENSION_BC, &
            NEXT_RECA)
      END DO
      DO LC = 1, DIMENSION_BC
         WRITE (UNIT_RES, REC=NEXT_RECA) BC_TYPE(LC)
         NEXT_RECA = NEXT_RECA + 1
      END DO


 1100 continue

!      call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!      call send_recv (flag,2)
      ! call gather (flag,arr1,root)
! To take care of filling the processor ghost layers with the correct values
      ! call scatter (flag,arr1,root)
!      call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
      if (myPE .ne. PE_IO) goto 1200
      call convert_to_io_i(arr1,arr2,ijkmax2)
      CALL OUT_BIN_512I (UNIT_RES, arr2, IJKMAX2, NEXT_RECA)


      WRITE (UNIT_RES, REC=NEXT_RECA) CYCLIC_X, CYCLIC_Y, CYCLIC_Z, CYCLIC_X_PD&
         , CYCLIC_Y_PD, CYCLIC_Z_PD, DELP_X, DELP_Y, DELP_Z, U_G0, U_S0, V_G0, &
         V_S0, W_G0, W_S0
      NEXT_RECA = NEXT_RECA + 1

! Version 01.09
! ------------------------------------------------------------------------
      WRITE (UNIT_RES, REC=NEXT_RECA) TIME, TSTOP, RES_DT, OUT_DT, &
         NLOG, NO_I, NO_J, NO_K, CALL_USR
      NEXT_RECA = NEXT_RECA + 1


      CALL OUT_BIN_512 (UNIT_RES, USR_DT, DIMENSION_USR, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, USR_X_W, DIMENSION_USR, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, USR_X_E, DIMENSION_USR, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, USR_Y_S, DIMENSION_USR, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, USR_Y_N, DIMENSION_USR, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, USR_Z_B, DIMENSION_USR, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, USR_Z_T, DIMENSION_USR, NEXT_RECA)
      DO LC = 1, DIMENSION_USR
         WRITE (UNIT_RES, REC=NEXT_RECA) USR_FORMAT(LC), USR_EXT(LC), USR_TYPE(&
            LC), USR_VAR(LC)
         NEXT_RECA = NEXT_RECA + 1
      END DO
      DO LC = 1, DIMENSION_IC
         WRITE (UNIT_RES, REC=NEXT_RECA) IC_TYPE(LC)
         NEXT_RECA = NEXT_RECA + 1
      END DO
      CALL OUT_BIN_512 (UNIT_RES, BC_DT_0, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_JET_G0, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_DT_H, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_JET_GH, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_DT_L, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_JET_GL, DIMENSION_BC, NEXT_RECA)



! Version 01.12
! ------------------------------------------------------------------------
      WRITE (UNIT_RES, REC=NEXT_RECA) P_REF, P_SCALE, UR_FAC, &
         TOL_RESID, DT_MAX, DT_MIN, DT_FAC, GRAVITY
      NEXT_RECA = NEXT_RECA + 1
      WRITE (UNIT_RES, REC=NEXT_RECA) LEQ_IT, LEQ_METHOD
      NEXT_RECA = NEXT_RECA + 1
      CALL OUT_BIN_512 (UNIT_RES, BC_HW_G, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_UW_G, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_VW_G, DIMENSION_BC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, BC_WW_G, DIMENSION_BC, NEXT_RECA)
      DO LC = 1, MMAX
         CALL OUT_BIN_512 (UNIT_RES, BC_HW_S(1,LC), DIMENSION_BC, NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES, BC_UW_S(1,LC), DIMENSION_BC, NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES, BC_VW_S(1,LC), DIMENSION_BC, NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES, BC_WW_S(1,LC), DIMENSION_BC, NEXT_RECA)
      END DO
      WRITE (UNIT_RES, REC=NEXT_RECA) MOMENTUM_X_EQ, MOMENTUM_Y_EQ, &
         MOMENTUM_Z_EQ, TOL_DIVERGE, DISCRETIZE, FULL_LOG
      NEXT_RECA = NEXT_RECA + 1

! Version 01.14
! ------------------------------------------------------------------------
      WRITE (UNIT_RES, REC=NEXT_RECA) DETECT_STALL
      NEXT_RECA = NEXT_RECA + 1


! Version 01.2
! ------------------------------------------------------------------------
      WRITE (UNIT_RES, REC=NEXT_RECA) NORM_G
      NEXT_RECA = NEXT_RECA + 1



!  Add new write statements above this line.  Remember to update NEXT_RECA.
!  Remember to change the version number near beginning of this subroutine.
!  Also modify READ_RES0.  The routines such as OUT_BIN_512 etc. writes
!  arrays dimensioned ARRAY(DIM).  So arrays dimensioned ARRAY(DIM1:DIM2)
!  should be passed as ARRAY(DIM1) and array length as DIM2-DIM1+1.
!---------------------------------------------------------------------------
      WRITE (UNIT_RES, REC=3) NEXT_RECA
      FLUSH (UNIT_RES)

 1200 continue


      deallocate (arr1)
      deallocate (arr2)

      RETURN
      END SUBROUTINE WRITE_RES0
