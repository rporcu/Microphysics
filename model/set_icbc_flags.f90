MODULE set_icbc_flags_module

      use bc, only: bc_x_e, bc_x_w, bc_i_e, bc_i_w
      use bc, only: bc_y_n, bc_y_s, bc_j_s, bc_j_n
      use bc, only: bc_z_t, bc_z_b, bc_k_b, bc_k_t
      use bc, only: dimension_bc, bc_type, bc_x_e, bc_defined
      use geometry, only: domlo, domhi
      use geometry    , only: cyclic_x, cyclic_y, cyclic_z
      use geometry    , only: cyclic_x_pd, cyclic_y_pd, cyclic_z_pd
      use ic, only: CYCP_, CYCL_, FSW_, UNDEF_CELL, NSW_, FLUID_
      use ic, only: MINF_, MOUT_, POUT_, OUTF_, PINF_, PSW_
      use mod_bc, only: mod_bc_i, mod_bc_j, mod_bc_k
      use param1, only: zero, is_defined, equal
      use run, only: run_type, IFILE_NAME

   CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: SET_ICBC_FLAG                                            !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE SET_ICBC_FLAG(slo,shi,flag)

      integer     , intent(in   ) :: slo(3),shi(3)

      integer, intent(inout) ::  flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      CALL INIT_ICBC_FLAG(slo,shi,flag)

      CALL SET_IC_FLAGS(slo,shi,flag)

      CALL SET_BC_FLAGS_WALL(slo,shi,flag)

      CALL SET_BC_FLAGS_FLOW(slo,shi,flag)

      END SUBROUTINE SET_ICBC_FLAG

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: INIT_ICBC_FLAG                                           !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE INIT_ICBC_FLAG(slo,shi,flag)

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3)

      integer, intent(inout) ::  flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      INTEGER :: I, J, K

! Initialize the icbc_flag array.
      DO k = slo(3),shi(3)
      DO j = slo(2),shi(2)
      DO i = slo(1),shi(1)

! Initialize the ICBC Flag
         if (run_type == 'NEW') then
            flag(i,j,k,1) = UNDEF_CELL
         else
            flag(i,j,k,1) = FLUID_
         endif

! If at domain boundaries then set default values (wall or, if
! specified, cyclic)
         IF(K==DOMLO(3)-1 .OR. K==DOMHI(3)+1)THEN
            IF (CYCLIC_Z_PD) THEN
               FLAG(i,j,k,1) = CYCP_
            ELSEIF (CYCLIC_Z) THEN
               FLAG(i,j,k,1) = CYCL_
            ELSE
               FLAG(i,j,k,1) = NSW_
            ENDIF
         ENDIF

         IF(J==DOMLO(2)-1 .OR. J==DOMHI(2)+1)THEN
            IF (CYCLIC_Y_PD) THEN
               FLAG(i,j,k,1) = CYCP_
            ELSEIF (CYCLIC_Y) THEN
               FLAG(i,j,k,1) = CYCL_
            ELSE
               FLAG(i,j,k,1) = NSW_
            ENDIF
         ENDIF

         IF(I==DOMLO(1)-1 .OR. I==DOMHI(1)+1)THEN
            IF (CYCLIC_X_PD) THEN
               FLAG(i,j,k,1) = CYCP_
            ELSEIF (CYCLIC_X) THEN
               FLAG(i,j,k,1) = CYCL_
            ELSE
               FLAG(i,j,k,1) = NSW_
            ENDIF
         ENDIF
! corner cells are wall cells
         if ((i==domlo(1)-1 .or. i==domhi(1)+1) .and. &
             (j==domlo(2)-1 .or. j==domhi(2)+1) .and. &
             (k==domlo(3)-1 .or. k==domhi(3)+1)) then
            IF (FLAG(i,j,k,1) /= FSW_) FLAG(i,j,k,1) = NSW_
         ENDIF

      ENDDO ! end do loop
      ENDDO ! end do loop
      ENDDO ! end do loop

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
      SUBROUTINE CHECK_ICBC_FLAG(slo,shi,lo,hi,flag)

      use compar, only: mype
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ival
      use ic, only: UNDEF_CELL
      use open_files_mod, only: open_pe_log
      use run, only: RUN_TYPE

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      integer, intent(inout) ::  flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      INTEGER :: I, J ,K, IER

      LOGICAL :: ERROR = .FALSE.

      IF(RUN_TYPE(1:3) /= 'NEW') RETURN

      CALL INIT_ERR_MSG("CHECK_ICBC_FLAG")

      ! First check for any errors.
      do K = lo(3)-1,hi(3)+1
         do J = lo(2)-1,hi(2)+1
            do I = lo(1)-1,hi(1)+1
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

         do K = lo(3)-1,hi(3)+1
            do J = lo(2)-1,hi(2)+1
               do I = lo(1)-1,hi(1)+1
                  IF (FLAG(i,j,k,1) == UNDEF_CELL) then
                     WRITE(ERR_MSG,1101) I, J, K
                     CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
                  ENDIF

               ENDDO
            ENDDO
         ENDDO

         WRITE(ERR_MSG, 1102) trim(IFILE_NAME)
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

 1102 format('Please correct the ',A,' file.')

      END SUBROUTINE CHECK_ICBC_FLAG

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_IC_FLAGS                                            !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Set the IC portions of the ICBC_Flag array.                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_IC_FLAGS(slo,shi,flag)

      use ic, only: IC_DEFINED
      use ic, only: IC_TYPE

      use ic, only: IC_I_W, IC_I_E
      use ic, only: IC_J_S, IC_J_N
      use ic, only: IC_K_B, IC_K_T
      use ic, only: FLUID_

      use param, only: dimension_ic

      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3)

      integer, intent(inout) ::  flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: ICV
      INTEGER :: I, J, K

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

      END SUBROUTINE SET_IC_FLAGS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_FLAGS_WALL                                       !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Find and validate i, j, k locations for walls BC's         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_BC_FLAGS_WALL(slo,shi,flag)

      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3)

      integer, intent(inout) ::  flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer :: istart, iend
      integer :: jstart, jend
      integer :: kstart, kend

      ! loop index
      integer :: BCV

!-----------------------------------------------

      ! Set the wall flags.
      DO BCV=1, DIMENSION_BC

         IF(.NOT.BC_DEFINED(BCV)) CYCLE

         IF(BC_TYPE(BCV)=='FREE_SLIP_WALL' .OR. &
            BC_TYPE(BCV)=='NO_SLIP_WALL'   .OR. &
            BC_TYPE(BCV)=='PAR_SLIP_WALL') THEN

            ! We don't want to write outside of the current grid
            istart = max(slo(1), bc_i_w(bcv))
            jstart = max(slo(2), bc_j_s(bcv))
            kstart = max(slo(3), bc_k_b(bcv))
            iend   = min(shi(1), bc_i_e(bcv))
            jend   = min(shi(2), bc_j_n(bcv))
            kend   = min(shi(3), bc_k_t(bcv))

            SELECT CASE (TRIM(BC_TYPE(BCV)))
               CASE('FREE_SLIP_WALL'); FLAG(istart:iend,jstart:jend,kstart:kend,1) = FSW_
               CASE('NO_SLIP_WALL'  ); FLAG(istart:iend,jstart:jend,kstart:kend,1) = NSW_
               CASE('PAR_SLIP_WALL' ); FLAG(istart:iend,jstart:jend,kstart:kend,1) = PSW_
            END SELECT

         ENDIF
      ENDDO

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
      SUBROUTINE SET_BC_FLAGS_FLOW(slo,shi,flag)

      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar
      use geometry    , only: cyclic_x, cyclic_y, cyclic_z
      use geometry, only: domhi

      USE open_files_mod, only: open_pe_log

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3)

      integer, intent(inout) ::  flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      ! loop/variable indices
      INTEGER :: bcv, i, j, k
      integer :: istart, iend
      integer :: jstart, jend
      integer :: kstart, kend

      INTEGER :: IER

      ! error indicator
      LOGICAL :: ERROR

      ! surface indictors
      LOGICAL :: X_CONSTANT, Y_CONSTANT, Z_CONSTANT

      CALL INIT_ERR_MSG("SET_BC_FLAGS_FLOW")

      ! Find the flow surfaces
      ERROR = .FALSE.

      DO BCV = 1, DIMENSION_BC

         IF(.NOT.BC_DEFINED(BCV)) CYCLE

         IF(BC_TYPE(BCV)=='MASS_INFLOW'  .OR. &
            BC_TYPE(BCV)=='MASS_OUTFLOW' .OR. &
            BC_TYPE(BCV)=='P_INFLOW'     .OR. &
            BC_TYPE(BCV)=='P_OUTFLOW'    .OR. &
            BC_TYPE(BCV)=='OUTFLOW') THEN

            X_CONSTANT = EQUAL(BC_X_W(BCV), BC_X_E(BCV))
            Y_CONSTANT = EQUAL(BC_Y_S(BCV), BC_Y_N(BCV))
            Z_CONSTANT = EQUAL(BC_Z_B(BCV), BC_Z_T(BCV))

            IF(X_CONSTANT .AND. IS_DEFINED(BC_X_W(BCV)))                &
               CALL MOD_BC_I(BCV,flag,slo,shi)

            IF(Y_CONSTANT .AND. IS_DEFINED(BC_Y_S(BCV)))                &
               CALL MOD_BC_J(BCV,flag,slo,shi)

            IF(Z_CONSTANT .AND. IS_DEFINED(BC_Z_B(BCV)))                &
               CALL MOD_BC_K(BCV,flag,slo,shi)

            ! Extend the boundaries for cyclic implementation
            IF (BC_I_W(BCV) == domlo(1) .and. &
                BC_I_E(BCV) == DOMHI(1) .and. &
                CYCLIC_X) then
                   BC_I_W(BCV) = 1
                   BC_I_E(BCV) = DOMHI(1)
            ENDIF
            IF(BC_J_S(BCV) == domlo(2) .and. &
               BC_J_N(BCV) == DOMHI(2) .and. &
               CYCLIC_Y) then
               BC_J_S(BCV) = 1
               BC_J_N(BCV) = DOMHI(2)
            ENDIF
            IF(BC_K_B(BCV) == domlo(3) .and. &
               BC_K_T(BCV) == DOMHI(3) .and. &
               CYCLIC_Z) then
               BC_K_B(BCV) = 1
               BC_K_T(BCV) = DOMHI(3)
            ENDIF

            ! Set add the BC to the FLAG. If a "non-wall" BC is found, then flag
            ! this as an error. The next triple-loop will take care of reporting the
            ! error.
            ERROR = .FALSE.

            ! We don't want to write outside of the current grid

            istart = max(slo(1), bc_i_w(bcv))
            jstart = max(slo(2), bc_j_s(bcv))
            kstart = max(slo(3), bc_k_b(bcv))
            iend   = min(shi(1), bc_i_e(bcv))
            jend   = min(shi(2), bc_j_n(bcv))
            kend   = min(shi(3), bc_k_t(bcv))

            DO K = kstart, kend
            DO J = jstart, jend
            DO I = istart, iend

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

               write(ERR_MSG,"('Please correct the ',A,' file.')") trim(IFILE_NAME)
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
