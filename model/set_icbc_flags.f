MODULE set_icbc_flags_module

      use bc, only: bc_x_e, bc_x_w, bc_i_e, bc_i_w
      use bc, only: bc_y_n, bc_y_s, bc_j_s, bc_j_n
      use bc, only: bc_z_t, bc_z_b, bc_k_b, bc_k_t
      use bc, only: dimension_bc, bc_type, bc_x_e, bc_defined
      use compar, only: iend3, jend3, kend3
      use compar, only: istart3, jstart3, kstart3
      use geometry    , only: cyclic_x, cyclic_y, cyclic_z
      use geometry    , only: cyclic_x_pd, cyclic_y_pd, cyclic_z_pd
      use geometry, only: imax2, jmax2, kmax2
      use geometry, only: imax3, jmax3, kmax3
      use geometry, only: imin2, jmin2, kmin2
      use geometry, only: imin3, jmin3, kmin3
      use ic, only: CYCP_, CYCL_, FSW_, UNDEF_CELL, NSW_, FLUID_
      use ic, only: MINF_, MOUT_, POUT_, OUTF_, PINF_, PSW_
      use mod_bc, only: mod_bc_i, mod_bc_j, mod_bc_k
      use param1, only: zero, undefined
      use run, only: run_type

   CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: SET_ICBC_FLAG                                            !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE SET_ICBC_FLAG(flag)

      integer, intent(inout) ::  flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,4)

      CALL INIT_ICBC_FLAG(flag)

      CALL SET_IC_FLAGS(flag)

      CALL SET_BC_FLAGS_WALL(flag)

      CALL SET_BC_FLAGS_FLOW(flag)

! Verify that ICBC flags are set for all fluid cells.
      CALL CHECK_ICBC_FLAG(flag)

      END SUBROUTINE SET_ICBC_FLAG



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: INIT_ICBC_FLAG                                           !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE INIT_ICBC_FLAG(flag)

      implicit none

      integer, intent(inout) ::  flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,4)

      INTEGER :: I, J, K

! Initialize the icbc_flag array.
      DO k = kstart3, kend3
      DO j = jstart3, jend3
      DO i = istart3, iend3

! Initialize the ICBC Flag
         if (run_type == 'NEW') then
            flag(i,j,k,1) = UNDEF_CELL
         else
            flag(i,j,k,1) = FLUID_
         endif

! If at domain boundaries then set default values (wall or, if
! specified, cyclic)
         IF(K==KMIN3 .OR. K==KMIN2 .OR. K==KMAX2 .OR. K==KMAX3)THEN
            IF (CYCLIC_Z_PD) THEN
               FLAG(i,j,k,1) = CYCP_
            ELSEIF (CYCLIC_Z) THEN
               FLAG(i,j,k,1) = CYCL_
            ELSE
               FLAG(i,j,k,1) = NSW_
            ENDIF
         ENDIF

         IF(J==JMIN3 .OR. J==JMIN2 .OR. J==JMAX2 .OR. J==JMAX3)THEN
            IF (CYCLIC_Y_PD) THEN
               FLAG(i,j,k,1) = CYCP_
            ELSEIF (CYCLIC_Y) THEN
               FLAG(i,j,k,1) = CYCL_
            ELSE
               FLAG(i,j,k,1) = NSW_
            ENDIF
         ENDIF

         IF(I==IMIN3 .OR. I==IMIN2 .OR. I==IMAX2 .OR. I==IMAX3)THEN
            IF (CYCLIC_X_PD) THEN
               FLAG(i,j,k,1) = CYCP_
            ELSEIF (CYCLIC_X) THEN
               FLAG(i,j,k,1) = CYCL_
            ELSE
               FLAG(i,j,k,1) = NSW_
            ENDIF
         ENDIF
! corner cells are wall cells
         IF ((I==IMIN3 .OR. I==IMIN2 .OR. I==IMAX2 .OR. I==IMAX3) .AND. &
             (J==JMIN3 .OR. J==JMIN2 .OR. J==JMAX2 .OR. J==JMIN3) .AND. &
             (K==KMIN3 .OR. K==KMIN2 .OR. K==KMAX2 .OR. K==KMAX3)) THEN
            IF (FLAG(i,j,k,1) /= FSW_) FLAG(i,j,k,1) = NSW_
         ENDIF

      ENDDO ! end do loop (i=istart3, iend3)
      ENDDO ! end do loop (j=jstart3, jend3)
      ENDDO ! end do loop (k=kstart3, kend3)

      RETURN

      END SUBROUTINE INIT_ICBC_FLAG



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_ICBC_FLAG                                         !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Verify that data was not given for undefined BC regions.   !
!  Note that the error message may be incomplete
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_ICBC_FLAG(flag)

      use compar, only: iend2, jend2, kend2
      use compar, only: istart2, jstart2, kstart2
      use compar, only: mype
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ival
      use ic, only: UNDEF_CELL
      use open_files_mod, only: open_pe_log
      use run, only: RUN_TYPE

      IMPLICIT NONE

      integer, intent(inout) ::  flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,4)

      INTEGER :: I, J ,K, IER

      LOGICAL :: ERROR = .FALSE.

      IF(RUN_TYPE(1:3) /= 'NEW') RETURN

      CALL INIT_ERR_MSG("CHECK_ICBC_FLAG")

      ! First check for any errors.
      DO K = kStart2, kEnd2
         DO J = jStart2, jEnd2
            DO I = iStart2, iEnd2
               IF (FLAG(i,j,k,1) == UNDEF_CELL) ERROR = .TRUE.
            ENDDO
         ENDDO
      ENDDO

! Sync up the error flag across all processes.
      ! CALL GLOBAL_ALL_OR(ERROR)

! If an error is detected, have each rank open a log file and write
! it's own message. Otherwise, we need to send all the data back to
! PE_IO and that's too much work!
      IF(ERROR) THEN

         CALL OPEN_PE_LOG(IER)

         WRITE(ERR_MSG, 1100) trim(iVal(myPE))
         CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)

         DO K = kStart2, kEnd2
            DO J = jStart2, jEnd2
               DO I = iStart2, iEnd2
                  IF (FLAG(i,j,k,1) == UNDEF_CELL) then
                     WRITE(ERR_MSG,1101) I, J, K
                     CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
                  ENDIF

               ENDDO
            ENDDO
         ENDDO

         WRITE(ERR_MSG, 1102)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)

      ELSE
! If no erros, sync up the ghost cell layers.
         ! CALL SEND_RECV(FLAG,1)
      ENDIF

! Clean up and return.
      CALL FINL_ERR_MSG

      RETURN

 1100 FORMAT('Error 1100 (PE ',A,') : No initial or boundary ',        &
         'condtions specified in','the following cells:',/             &
         '    I       J       K')

 1101 FORMAT(I5,3X,I5,3X,I5)

 1102 FORMAT('Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_ICBC_FLAG



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_IC_FLAGS                                            !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Set the IC portions of the ICBC_Flag array.                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_IC_FLAGS(flag)

      use ic, only: IC_DEFINED
      use ic, only: IC_TYPE

      use ic, only: IC_I_W, IC_I_E
      use ic, only: IC_J_S, IC_J_N
      use ic, only: IC_K_B, IC_K_T
      use ic, only: FLUID_

      use param, only: dimension_ic

      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar

      IMPLICIT NONE

      integer, intent(inout) ::  flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,4)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: ICV
      INTEGER :: I, J, K

      CALL INIT_ERR_MSG("SET_IC_FLAGS")


      IC_LP: DO ICV=1, DIMENSION_IC

         IF(.NOT.IC_DEFINED(ICV)) CYCLE IC_LP

! Skip checks for PATCH restarts.
         IF (IC_TYPE(ICV) == 'PATCH') CYCLE IC_LP

!  Set ICBC flag
         DO K = IC_K_B(ICV), IC_K_T(ICV)
            DO J = IC_J_S(ICV), IC_J_N(ICV)
               DO I = IC_I_W(ICV), IC_I_E(ICV)
                  flag(i,j,k,1) = FLUID_
               ENDDO
            ENDDO
         ENDDO


      ENDDO IC_LP

! Update the ICBC flag on ghost cells.
      ! CALL SEND_RECV(FLAG, 1)


! Clean up and return.
      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE SET_IC_FLAGS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_FLAGS_WALL                                       !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Find and validate i, j, k locations for walls BC's         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_BC_FLAGS_WALL(flag)

      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar

      IMPLICIT NONE

      integer, intent(inout) ::  flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,4)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! loop/variable indices
      INTEGER :: I , J , K
! loop index
      INTEGER :: BCV

!-----------------------------------------------

      CALL INIT_ERR_MSG("SET_BC_FLAGS_WALL")

! Set the wall flags.
      DO BCV=1, DIMENSION_BC
         IF(.NOT.BC_DEFINED(BCV)) CYCLE

         IF(BC_TYPE(BCV)=='FREE_SLIP_WALL' .OR. &
            BC_TYPE(BCV)=='NO_SLIP_WALL'   .OR. &
            BC_TYPE(BCV)=='PAR_SLIP_WALL') THEN

            DO K = BC_K_B(BCV), BC_K_T(BCV)
            DO J = BC_J_S(BCV), BC_J_N(BCV)
            DO I = BC_I_W(BCV), BC_I_E(BCV)

               SELECT CASE (TRIM(BC_TYPE(BCV)))
                  CASE('FREE_SLIP_WALL'); FLAG(i,j,k,1) = FSW_
                  CASE('NO_SLIP_WALL');   FLAG(i,j,k,1) = NSW_
                  CASE('PAR_SLIP_WALL');  FLAG(i,j,k,1) = PSW_
               END SELECT
            ENDDO
            ENDDO
            ENDDO

         ENDIF
      ENDDO

      ! CALL SEND_RECV(FLAG,1)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_BC_FLAGS_WALL

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_FLAGS_FLOW                                       !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Find and validate i, j, k locations for flow BC's. Also    !
!           set value of bc_plane for flow BC's.                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_BC_FLAGS_FLOW(flag)

      use compar       , only: nodesi, nodesj, nodesk
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar
      use geometry     , only: imax2, jmax2, kmax2
      use geometry    , only: cyclic_x, cyclic_y, cyclic_z

      USE open_files_mod, only: open_pe_log

      IMPLICIT NONE

      integer, intent(inout) ::  flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,4)

! loop/variable indices
      INTEGER :: BCV, I, J, K

      INTEGER :: IER

! error indicator
      LOGICAL :: ERROR
! surface indictors
      LOGICAL :: X_CONSTANT, Y_CONSTANT, Z_CONSTANT

      CALL INIT_ERR_MSG("SET_BC_FLAGS_FLOW")


! FIND THE FLOW SURFACES
      ERROR = .FALSE.

      DO BCV = 1, DIMENSION_BC

         IF(.NOT.BC_DEFINED(BCV)) CYCLE

         IF(BC_TYPE(BCV)=='MASS_INFLOW'  .OR. &
            BC_TYPE(BCV)=='MASS_OUTFLOW' .OR. &
            BC_TYPE(BCV)=='P_INFLOW'     .OR. &
            BC_TYPE(BCV)=='P_OUTFLOW'    .OR. &
            BC_TYPE(BCV)=='OUTFLOW') THEN

            X_CONSTANT = (BC_X_W(BCV) == BC_X_E(BCV))
            Y_CONSTANT = (BC_Y_S(BCV) == BC_Y_N(BCV))
            Z_CONSTANT = (BC_Z_B(BCV) == BC_Z_T(BCV))


            write(*,*) 'call mod_bc_i',bcv
            IF(X_CONSTANT .AND. BC_X_W(BCV)/=UNDEFINED)                &
               CALL MOD_BC_I(BCV,flag)

            write(*,*) 'call mod_bc_j',bcv
            IF(Y_CONSTANT .AND. BC_Y_S(BCV)/=UNDEFINED)                &
               CALL MOD_BC_J(BCV,flag)

            write(*,*) 'call mod_bc_k',bcv
            IF(Z_CONSTANT .AND. BC_Z_B(BCV)/=UNDEFINED)                &
               CALL MOD_BC_K(BCV,flag)

! Extend the boundaries for cyclic implementation
            IF(BC_I_W(BCV) == 2 .AND. BC_I_E(BCV) == (IMAX2 - 1) .AND. &
               CYCLIC_X .AND. NODESI > 1) THEN
                   BC_I_W(BCV) = 1
                   BC_I_E(BCV) = IMAX2
            ENDIF
            IF(BC_J_S(BCV) == 2 .AND. BC_J_N(BCV) == (JMAX2 - 1) .AND. &
               CYCLIC_Y .AND. NODESJ > 1) THEN
               BC_J_S(BCV) = 1
               BC_J_N(BCV) = JMAX2
            ENDIF
            IF(BC_K_B(BCV) == 2 .AND. BC_K_T(BCV) == (KMAX2 - 1) .AND. &
               CYCLIC_Z .AND. NODESK > 1) THEN
               BC_K_B(BCV) = 1
               BC_K_T(BCV) = KMAX2
            ENDIF

! Set add the BC to the FLAG. If a "non-wall" BC is found, then flag
! this as an error. The next triple-loop will take care of reporting the
! error.
            ERROR = .FALSE.
            DO K = BC_K_B(BCV), BC_K_T(BCV)
            DO J = BC_J_S(BCV), BC_J_N(BCV)
            DO I = BC_I_W(BCV), BC_I_E(BCV)

! Verify that the FLOW BC is overwriting a wall.
               IF(flag(i,j,k,1) == NSW_ .or. &
                  flag(i,j,k,1) == FSW_ .or. &
                  flag(i,j,k,1) == PSW_) then

                  SELECT CASE (TRIM(BC_TYPE(BCV)))
                     CASE ('P_OUTFLOW');    FLAG(i,j,k,1) = POUT_
                     CASE ('MASS_INFLOW');  FLAG(i,j,k,1) = MINF_
                     CASE ('MASS_OUTFLOW'); FLAG(i,j,k,1) = MOUT_
                     CASE ('OUTFLOW');      FLAG(i,j,k,1) = OUTF_
                     CASE ('P_INFLOW');     FLAG(i,j,k,1) = PINF_
                  END SELECT

               ELSE
                  ERROR = .TRUE.
               ENDIF

            ENDDO
            ENDDO
            ENDDO

! Sync the error flag over all ranks.
            ! CALL GLOBAL_ALL_OR(ERROR)

! Report errors and exit.
            IF(ERROR)THEN

               CALL OPEN_PE_LOG(IER)

               WRITE(ERR_MSG, 1200) BCV
               CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)

 1200 FORMAT('Error 1200: Boundary condition ',I3,' overlaps with ',&
         'another BC.',2/7x,'I',7x,'J',7x,'K',3x,'ICBC')

               DO K = BC_K_B(BCV), BC_K_T(BCV)
               DO J = BC_J_S(BCV), BC_J_N(BCV)
               DO I = BC_I_W(BCV), BC_I_E(BCV)

! Verify that the FLOW BC is overwriting a wall.
                  IF(flag(i,j,k,1) /= NSW_ .and. &
                     flag(i,j,k,1) /= FSW_ .and. &
                     flag(i,j,k,1) /= PSW_) then
                     WRITE(ERR_MSG, 1201) I,J,K, FLAG(i,j,k,1)
                     CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
                  ENDIF

 1201 FORMAT(1x,3(2x,I6),3x,A3)

               ENDDO
               ENDDO
               ENDDO

               WRITE(ERR_MSG,"('Please correct the mfix.dat file.')")
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)

            ENDIF ! IF(ERROR)
         ENDIF ! IF(not a wall BC)
      ENDDO ! BC Loop

! Sync the ICBC flag across ghost layers
      ! CALL SEND_RECV(FLAG,1)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_BC_FLAGS_FLOW
END MODULE set_icbc_flags_module
