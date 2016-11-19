!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc0                                                 C
!  Purpose: This subroutine does the initial setting of all boundary   C
!  conditions. The user specifications of the boundary conditions are  C
!  checked for veracity in various check_data/ routines:               C
!  (e.g., check_boundary_conditions).                                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_BC0

! Modules
!--------------------------------------------------------------------//
      use bc, only: bc_type, bc_defined
      use fldvar, only: p_g

      use param, only: dimension_bc
      use param1, only: undefined_i

      use mpi_utility
      use sendrecv
      IMPLICIT NONE

! Local variables
!--------------------------------------------------------------------//
! Local index for boundary condition
      INTEGER ::  L
!--------------------------------------------------------------------//

! Incompressible cases require that Ppg specified for one cell.
! The following attempts to pick an appropriate cell.
      CALL SET_IJK_P_G


      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN

            SELECT CASE (TRIM(BC_TYPE(L)))

            CASE ('FREE_SLIP_WALL')
            CASE ('NO_SLIP_WALL')
            CASE ('PAR_SLIP_WALL')
            CASE ('P_OUTFLOW')
               CALL SET_BC0_INIT_BCDT_CALCS(L)
               CALL SET_BC0_OUTFLOW(L)
            CASE ('MASS_OUTFLOW')
               CALL SET_BC0_INIT_BCDT_CALCS(L)
               CALL SET_BC0_INFLOW(L)
            CASE ('OUTFLOW')
               CALL SET_BC0_INIT_BCDT_CALCS(L)
               CALL SET_BC0_OUTFLOW(L)
            CASE ('MASS_INFLOW')
               CALL SET_BC0_INIT_JET(L)
               CALL SET_BC0_INFLOW(L)
            CASE ('P_INFLOW')
               CALL SET_BC0_INFLOW(L)
            END SELECT
         ENDIF
      ENDDO

! Make T_g nonzero in k=0,1 ghost layers when k-decomposition employed
      call send_recv(P_G,2)

      RETURN
      END SUBROUTINE SET_BC0



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc0_outflow                                         C
!  Purpose: Set the initial settings of the boundary conditions (if    C
!  they are defined) for pressure outflow (PO) or outflow (O)          C
!  boundary types.                                                     C
!                                                                      C
!  Comments: For a new run the field variables are undefined in the    C
!  boundary cell locations, while for a restart run the field variable C
!  may have an existing value based on the preceding simulation.       C
!  Regardless, a user defined BC value will supercede any existing     C
!  value.                                                              C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_BC0_OUTFLOW(BCV)

! Modules
!--------------------------------------------------------------------//
      use bc, only: bc_k_b, bc_k_t
      use bc, only: bc_j_s, bc_j_n
      use bc, only: bc_i_w, bc_i_e
      use bc, only: bc_p_g
      use bc, only: bc_ep_g

      use constant, only: pi

      use fldvar, only: p_g
      use fldvar, only: ep_g

      use param1, only: undefined, zero
      use scales, only: scale_pressure
      use toleranc, only: tmin

      use functions, only: is_on_mype_plus2layers
      use compar, only: dead_cell_at
      use boundfunijk, only: bound_funijk
      IMPLICIT NONE

! Dummy arguments
!--------------------------------------------------------------------//
! index of boundary condition
      INTEGER, INTENT(IN) :: BCV

! Local variables
!--------------------------------------------------------------------//
! indices
      INTEGER :: I, J, K, IJK
!--------------------------------------------------------------------//


      DO K = BC_K_B(BCV), BC_K_T(BCV)
      DO J = BC_J_S(BCV), BC_J_N(BCV)
      DO I = BC_I_W(BCV), BC_I_E(BCV)
         IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
         IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
         IJK = BOUND_FUNIJK(I,J,K)

         P_G(IJK) = SCALE_PRESSURE(BC_P_G(BCV))
         IF (BC_EP_G(BCV) /= UNDEFINED) EP_G(IJK) = BC_EP_G(BCV)

      ENDDO   ! do i
      ENDDO   ! do j
      ENDDO   ! do k

      RETURN
      END SUBROUTINE SET_BC0_OUTFLOW

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc0_init_jet                                        C
!  Purpose: initializing time dependent jet conditions                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_BC0_INIT_JET(BCV)

! Modules
!--------------------------------------------------------------------//
      use bc, only: bc_plane
      use bc, only: bc_jet_g, bc_jet_g0
      use bc, only: bc_dt_0, bc_time
      use bc, only: bc_u_g, bc_v_g, bc_w_g
      use param1, only: undefined
      use run, only: time
      IMPLICIT NONE

! Dummy arguments
!--------------------------------------------------------------------//
! index of boundary
      INTEGER, INTENT(IN) :: BCV
!--------------------------------------------------------------------//

      BC_JET_G(BCV) = UNDEFINED
      IF (BC_DT_0(BCV) /= UNDEFINED) THEN
         BC_TIME(BCV) = TIME + BC_DT_0(BCV)
         BC_JET_G(BCV) = BC_JET_G0(BCV)
         IF (BC_JET_G(BCV) /= UNDEFINED) THEN
            SELECT CASE (TRIM(BC_PLANE(BCV)))
            CASE ('W', 'E')
               BC_U_G(BCV) = BC_JET_G(BCV)
            CASE ('S', 'N')
               BC_V_G(BCV) = BC_JET_G(BCV)
            CASE ('B', 'T')
               BC_W_G(BCV) = BC_JET_G(BCV)
            END SELECT
         ENDIF
      ELSE
         BC_TIME(BCV) = UNDEFINED
      ENDIF
      RETURN
      END SUBROUTINE SET_BC0_INIT_JET



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc0_inflow                                          C
!  Purpose: Set the initial settings of the boundary conditions        C
!  pressure inflow (PI) and. mass inflow (MI) boundary types. Also do  C
!  mass outflow (MO) boundary types... due to velocity...              C
!                                                                      C
!  Comments: Unlike the treament of PO or O boundary types no checks   C
!  are made for these boundary types to determine whether the given    C
!  BC value is defined before it is assigned to the field variable.    C
!  However, the corresponding check routines generally ensure such BC  C
!  quantities are defined for MI or PI boundaries if they are needed   C
!  for the simulation.                                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_BC0_INFLOW(BCV)

! Modules
!--------------------------------------------------------------------//
      use bc, only: bc_plane
      use bc, only: bc_k_b, bc_k_t
      use bc, only: bc_j_s, bc_j_n
      use bc, only: bc_i_w, bc_i_e
      use bc, only: bc_u_g, bc_v_g, bc_w_g
      use bc, only: bc_p_g
      use bc, only: bc_ep_g

      use constant, only: pi

      use fldvar, only: u_g, v_g, w_g
      use fldvar, only: p_g
      use fldvar, only: ep_g

      use param1, only: zero
      use scales, only: scale_pressure

      use indices, only: im1, jm1, km1
      use functions, only: is_on_mype_plus2layers
      use boundfunijk, only: bound_funijk
      use compar, only: dead_cell_at
      IMPLICIT NONE

! Dummy arguments
!--------------------------------------------------------------------//
! index for boundary condition
      INTEGER, INTENT(IN) :: BCV

! Local variables
!--------------------------------------------------------------------//
! indices
      INTEGER :: I, J, K, IJK
! ijk index for setting normal component of velocity
      INTEGER :: FIJK
!--------------------------------------------------------------------//

      DO K = BC_K_B(BCV), BC_K_T(BCV)
      DO J = BC_J_S(BCV), BC_J_N(BCV)
      DO I = BC_I_W(BCV), BC_I_E(BCV)
         IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
         IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
         IJK = BOUND_FUNIJK(I,J,K)

         P_G(IJK) = SCALE_PRESSURE(BC_P_G(BCV))
         EP_G(IJK) = BC_EP_G(BCV)

         U_G(IJK) = BC_U_G(BCV)
         V_G(IJK) = BC_V_G(BCV)
         W_G(IJK) = BC_W_G(BCV)

! When the boundary plane is located on the E, N, T side of the domain
! (fluid cell is located w, s, b), set the component of velocity normal
! to the boundary plane of the adjacent fluid cell
         SELECT CASE (TRIM(BC_PLANE(BCV)))
            CASE ('W')
               FIJK = BOUND_FUNIJK(IM1(I),J,K)
               U_G(FIJK) = BC_U_G(BCV)
            CASE ('S')
               FIJK = BOUND_FUNIJK(I,JM1(J),K)
               V_G(FIJK) = BC_V_G(BCV)
            CASE ('B')
               FIJK = BOUND_FUNIJK(I,J,KM1(K))
               W_G(FIJK) = BC_W_G(BCV)
         END SELECT

      ENDDO   ! do i
      ENDDO   ! do j
      ENDDO   ! do k

      RETURN
      END SUBROUTINE SET_BC0_INFLOW




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc0_init_bcdt_calcs                                 C
!  Purpose: initializing time dependent outflow calculations for       C
!  modifying outflow conditions (MO type) or simple reporting          C
!  outflow conditions (PO or O types)                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_BC0_INIT_BCDT_CALCS(BCV)

! Modules
!--------------------------------------------------------------------//
      use bc, only: bc_dt_0, bc_time
      use bc, only: bc_mout_g
      use bc, only: bc_out_n
      use bc, only: bc_vout_g
      use param1, only: undefined, zero
      use run, only: time
      IMPLICIT NONE
! Dummy arguments
!--------------------------------------------------------------------//
! index of boundary
      INTEGER, INTENT(IN) :: BCV
!--------------------------------------------------------------------//

! initializing for time dependent outflow reporting calculation
      IF (BC_DT_0(BCV) /= UNDEFINED) THEN
         BC_TIME(BCV) = TIME + BC_DT_0(BCV)
         BC_OUT_N(BCV) = 0
         BC_MOUT_G(BCV) = ZERO
         BC_VOUT_G(BCV) = ZERO
      ELSE
         BC_TIME(BCV) = UNDEFINED
      ENDIF
      RETURN
      END SUBROUTINE SET_BC0_INIT_BCDT_CALCS




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_IJK_P_G                                             !
!  Purpose: Pick an appropriate control volume to specify Ppg.         !
!                                                                      !
!  Author: J. Musser                                  Date: 07-Nov-13  !
!  Reviewer:                                          Date:            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_IJK_P_G

! IJK location where Ppg is fixed.
      use bc, only: IJK_P_G
! Specified constant gas density.
      use fldvar, only: RO_G0

      use geometry, only: CYCLIC_X, CYCLIC_X_PD, CYCLIC_X_MF
      use geometry, only: CYCLIC_Y, CYCLIC_Y_PD, CYCLIC_Y_MF
      use geometry, only: CYCLIC_Z, CYCLIC_Z_PD, CYCLIC_Z_MF
      use geometry, only: iMAX1, iMin1
      use geometry, only: jMAX1, jMin1
      use geometry, only: kMAX1, kMin1

      use geometry, only: do_K

      use funits, only: DMP_LOG

      use bc, only: BC_DEFINED
      use bc, only: BC_TYPE

! MFIX Runtime parameters:
      use param, only: DIMENSION_BC
      use param1, only: UNDEFINED
      use param1, only: UNDEFINED_I

      use mpi_utility

      implicit none
!--------------------------------------------------------------------//
      INTEGER :: BCV

      CHARACTER(len=7) :: Map
      CHARACTER(len=128) :: lMsg

      INTEGER :: l3
      INTEGER :: l2, u2
      INTEGER :: l1, u1

      LOGICAL, parameter :: setDBG = .FALSE.
      LOGICAL :: dFlag
      INTEGER :: iErr

!--------------------------------------------------------------------//

      dFlag = (DMP_LOG .AND. setDBG)

! Initialize.
      iErr = 0
      IJK_P_G = UNDEFINED_I

! This is not needed for compressible cases.
      IF(RO_G0 == UNDEFINED) THEN
         IF(dFlag) write(*,"(3x,A)")                                   &
            'Compressible: IJK_P_g remaining undefined.'
         return
      ELSEIF(RO_G0 == 0.0d0) THEN
         IF(dFlag) write(*,"(3x,A)")                                   &
            'No gas phase: IJK_P_g remaining undefined.'
         return
      ENDIF

! If there are no cyclic boundaries, look for a pressure outflow.
      lpBCV: DO BCV = 1, DIMENSION_BC
         IF (.NOT.BC_DEFINED(BCV)) cycle lpBCV
         IF (BC_TYPE(BCV) == 'P_OUTFLOW' .OR. &
             BC_TYPE(BCV) == 'P_INFLOW') THEN
            IF(dFlag) write(*,"(3x,A)")                                &
               'Outflow PC defined: IJK_P_g remaining undefined.'
            RETURN
         ENDIF
      ENDDO lpBCV

! Initialize.
         l3 = UNDEFINED_I
         l2 = UNDEFINED_I; u2=l2
         l1 = UNDEFINED_I; u1=l1

! If there are cyclic boundaries, flag a cell along the positive
! domain extreme in the cyclic direction (e.g., JMAX1).
      IF(CYCLIC_Y .OR. CYCLIC_Y_PD .OR. CYCLIC_Y_MF) THEN

         Map = 'JKI_MAP'
         l3 = JMAX1
         l2 = KMIN1;  u2 = KMAX1
         l1 = IMIN1;  u1 = IMAX1
         lMsg='Cyclic in Y'

      ELSEIF(CYCLIC_X .OR. CYCLIC_X_PD .OR. CYCLIC_X_MF) THEN

         Map = 'IKJ_MAP'
         l3 = IMAX1
         l2 = KMIN1;  u2 = KMAX1
         l1 = JMIN1;  u1 = JMAX1
         lMsg='Cyclic in X'

      ELSEIF(CYCLIC_Z .OR. CYCLIC_Z_PD .OR. CYCLIC_Z_MF) THEN

         Map = 'KIJ_MAP'
         l3 = KMAX1
         l2 = IMIN1;  u2 = IMAX1
         l1 = JMIN1;  u1 = JMAX1
         lMsg='Cyclic in Z'

      ENDIF

! No cyclic boundaries or pressure outflows. The IJ plane is used in
! this case to maximize search region for 2D problems.
      IF(l3 == UNDEFINED_I) THEN
         Map = 'KIJ_MAP'
         l3 = merge(max((KMAX1-KMIN1)/2+1,2), KMIN1, do_K)
         l2 = IMIN1;  u2 = IMAX1
         l1 = JMIN1;  u1 = JMAX1
         lMsg='Center of domain'
      ENDIF

! Debugging messages.
      IF(dFlag) THEN
         write(*,"(3/,3x,'Map: ',A)") Map
         write(*,"(/5x,'l3:',2x,I4)") l3
         write(*,"( 5x,'l2:',2(2x,I4))") l2, u2
         write(*,"( 5x,'l1:',2(2x,I4))") l1, u1
         write(*,"( 5x,'Msg: ',A)") trim(lMsg)
      ENDIF

! Invoke the search routine.
      CALL IJK_Pg_SEARCH(l3, l2, u2, l1, u1, MAP, dFlag, iErr)

      IF(iErr == 0) RETURN

! Error management.
      IF(DMP_LOG) THEN
         SELECT CASE (iErr)
         CASE ( 1001);  WRITE(UNIT_LOG, 1001); WRITE(*,1001)
         CASE ( 2000);  WRITE(UNIT_LOG, 2000); WRITE(*,2000)
         CASE ( 2001);  WRITE(UNIT_LOG, 2001); WRITE(*,2001)
         CASE ( 2002);  WRITE(UNIT_LOG, 2002); WRITE(*,2002)
         CASE DEFAULT
            WRITE(UNIT_LOG, 1000) iErr
            WRITE(*,1000) iErr
         END SELECT

         WRITE(UNIT_LOG, 9000) MAP(1:1), l3, MAP(2:2),                 &
            l2, u2, MAP(3:3), l1, u1
         WRITE(*, 9000) MAP(1:1), l3, MAP(2:2),                        &
            l2, u2, MAP(3:3), l1, u1

         WRITE(*, 9999)
         WRITE(UNIT_LOG, 9999)

      ENDIF


      CALL MFIX_EXIT(myPE)


 1000 FORMAT(//1X,70('*')/' From: SET_IJK_Pg',/,                       &
         ' Error 1000: Unknown error reported. x', I4.4)

 1001 FORMAT(//1X,70('*')/' From: SET_IJK_Pg',/,                       &
         ' Error 1001: Invalid mapping function.')

 2000 FORMAT(//1X,70('*')/' From: SET_IJK_Pg > IJK_Pg_SEARCH',/,       &
         ' Error 2000: Unknown error reported from IJK_Pg_SEARCH.')

 2001 FORMAT(//1X,70('*')/' From: SET_IJK_Pg > IJK_Pg_SEARCH',/,       &
         ' Error 2001: Unable to locate fluid cell in search region.')

 2002 FORMAT(//1X,70('*')/' From: SET_IJK_Pg > IJK_Pg_SEARCH',/,       &
         ' Error 2002: Unable to locate fluid cell owner.')

 9000 FORMAT(/' Search plane information:',/,3x,A1,': ',I8,            &
          2(/3x,A1,': ',I8,' x ',I8))

 9999 FORMAT(/' Fatal Error --> Invoking MFIX_EXIT',/1x,70('*'),2/)

      END SUBROUTINE SET_IJK_P_G

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Author: J. Musser                                  Date: 07-Nov-13  !
!  Reviewer:                                          Date:            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE IJK_Pg_SEARCH(ll3, ll2, lu2, ll1, lu1, lMAP,          &
         ldFlag, iErr)

! Modules
!--------------------------------------------------------------------//
! IJK location where Ppg is fixed.
      use bc, only: IJK_P_g
      use indices
      use param1, only: undefined_i
      use mpi_utility
      use functions
      implicit none

! Dummy arguments
!--------------------------------------------------------------------//
      INTEGER, INTENT(IN)  :: ll3
      INTEGER, INTENT(IN)  :: ll2, lu2
      INTEGER, INTENT(IN)  :: ll1, lu1
      LOGICAL, INTENT(IN) :: ldFlag
      INTEGER, INTENT(OUT)  :: iErr
      CHARACTER(len=*), INTENT(IN) :: lMAP

! Local variables
!--------------------------------------------------------------------//
      INTEGER :: lc2, lS2, lE2
      INTEGER :: lc1, lS1, lE1
      INTEGER :: I, J, K, IJK
      LOGICAL :: recheck
      INTEGER :: IJK_Pg_Owner, proc
      INTEGER :: gIJK(0:numPEs-1)
      INTEGER :: I_J_K_Pg(3)
      INTEGER :: lpCnt

      CHARACTER(len=32) :: cInt

!--------------------------------------------------------------------//

! Initialize Error Flag
      iErr = 2000

! Initialize the Owner ID
      IJK_Pg_Owner = UNDEFINED_I

! Set the initial search region, a single cell.
      lS1 = ll1 + (lu1-ll1)/2 + 1; lE1 = lS1
      lS2 = ll2 + (lu2-ll2)/2 + 1; lE2 = lS2

      lpCnt = 1
      recheck = .TRUE.
      do while(recheck)

! Initialize the global IJK array to zero. Resetting this array inside
! this do-loop is most likely overkill. This loop should only cycle
! if gIJK is zero.
         gIJK = 0

! Report debugging information for the search region.
         if(ldFlag) then
            write(*,"(/5x,'Pass: ',I4)") lpCnt
            write(*,"( 5x,'lp2 bounds:',2(2x,I4))")lS2, lE2
            write(*,"( 5x,'lp1 bounds:',2(2x,I4))")lS1, lE1
         endif

         lp2: do lc2 = lS2, lE2
         lp1: do lc1 = lS1, lE1
! Map the loop counters to I/J/K indices.
            SELECT CASE (lMap)
            CASE ('JKI_MAP')
               I=lc1; J=ll3; K=lc2
            CASE ('IKJ_MAP')
               I=ll3; J=lc1; K=lc2
            CASE ('KIJ_MAP')
               I=lc2; J=lc1; K=ll3
            CASE DEFAULT
               iErr = 1001
            END SELECT

! Only the rank that owns this I/J/K proceeds.
            if(.NOT.IS_ON_myPE_owns(I,J,K)) cycle
! Calculate the triple loop index.
            IJK = funijk(I,J,K)
! If there is fluid at this location, store the IJK and exit loops.
            if(fluid_at(IJK)) then
               gIJK(myPE) = IJK
               exit lp2
            endif
         enddo lp1
         enddo lp2

! Sync gIJK across all processes. Select the lowest ranked process that
! has gIJK defined. This choice is arbitray and doesn't really matter.
! It just needs to be consistent.
         CALL global_all_sum(gIJK)
         proc_lp: do proc=0, numPEs-1
            if(gIJK(proc) /= 0) then
               IJK_P_g = gIJK(proc)
               IJK_Pg_Owner = proc
               recheck = .FALSE.
               exit proc_lp
            endif
         enddo proc_lp

! If the proceeding section did not find a fluid cell, expand the search
! area and try again.
         if(recheck) then
            if(lS1 > ll1 .OR. lE1 < lu1 .OR.                           &
               lS2 > ll2 .OR. lE2 < lu2) then
! Expand the 1-axis
               lS1 = max((lS1-1), ll1)
               lE1 = min((lE1+1), lu1)
! Expand the 2-axis
               lS2 = max((lS2-1), ll2)
               lE2 = min((lE2+1), lu2)
! The entire seach plane was checked with no fluid cell identified.
! Force IJK_P_g to undefined for later error checking.
            else
               recheck = .FALSE.
               IJK_P_g = UNDEFINED_I
            endif
         endif
      enddo

! Verify that one fluid cell was detected. Otherwise flag the possible
! errors and return.
      if(IJK_P_G == UNDEFINED_I) then
         iErr = 2001
         return
      elseif(IJK_Pg_Owner == UNDEFINED_I) then
         iErr = 2002
         return
      endif

! The owner if the IJK_Pg gets the global I/J/K values and sends
! them to all ranks.
      if(myPE == IJK_Pg_Owner) then
         I_J_K_Pg(1) = I_OF(IJK_P_G)
         I_J_K_Pg(2) = J_OF(IJK_P_G)
         I_J_K_Pg(3) = K_OF(IJK_P_G)
      endif
      CALL BCAST(I_J_K_Pg, IJK_Pg_Owner)

      I = I_J_K_Pg(1)
      J = I_J_K_Pg(2)
      K = I_J_K_Pg(3)

! If debugging, have PE_IO report some information before the
! data is overwritten.
      if(ldFlag) then
         write(*,"(/3x,'IJK_P_g successfully identified!')")
         cInt=''; write(cInt,*) IJK_Pg_Owner
         write(*,"( 5x,'Owner Rank: ',A)")trim(adjustl(cInt))
         cInt=''; write(cInt,*) IJK_P_G
         write(*,"(5x, 'IJK: ',A)", advance='no') trim(adjustl(cInt))
         write(*,"(' :: ')", advance='no')
         cInt=''; write(cInt,*) I
         write(*,"('(',A)",advance='no') trim(adjustl(cInt))
         cInt=''; write(cInt,*) J
         write(*,"(',',A)",advance='no') trim(adjustl(cInt))
         cInt=''; write(cInt,*) K
         write(*,"(',',A,')',2/)") trim(adjustl(cInt))
      endif

! Ranks that 'see' IJK_P_g store their local IJK value. Everyone else
! resets IJK_P_g to UNDEFINED_I. This removes the need for getting
! I/J/K values later on in source_PPg.
!      IJK_P_g = merge(funijk(I,J,K), UNDEFINED_I,                      &
!         IS_ON_myPE_plus2layers(I,J,K))

      IF(IS_ON_myPE_plus2layers(I,J,K)) THEN
         IJK_P_g = funijk(I,J,K)
      ELSE
         IJK_P_g = UNDEFINED_I
      ENDIF

      IERR = 0
      RETURN
      END SUBROUTINE IJK_Pg_SEARCH
