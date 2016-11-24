!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_DEM_MO                                           !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 23-Nov-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_BC_DEM_MO

      use bc, only: BC_PLANE
      use bc, only: BC_X_w, BC_X_e, BC_Y_s, BC_Y_n, BC_Z_b, BC_Z_t

      use des_bc, only: DEM_BCMO, DEM_BCMO_MAP, DEM_BCMO_IJK
      use des_bc, only: DEM_BCMO_IJKSTART, DEM_BCMO_IJKEND

      use funits, only: DMP_LOG

      use error_manager
      use functions

      use desgrid, only: DG_FUNIJK
      use desgrid, only: IofPOS, JofPOS, KofPOS
      use desgrid, only: dg_is_ON_myPE_plus1layers

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: BCV, BCV_I      ! BC loop counter

      INTEGER :: LC

      LOGICAL, parameter :: setDBG = .FALSE.
      LOGICAL :: dFlag

      INTEGER :: MAX_CELLS, BND1, BND2

      INTEGER, ALLOCATABLE :: LOC_DEM_BCMO_IJK(:)

      INTEGER :: I,J,K,IJK
      INTEGER :: I_w, I_e, J_s, J_n, K_b, K_t

      CALL INIT_ERR_MSG("SET_BC_DEM_MO")


! Initialize the data structures:
      allocate( DEM_BCMO_IJKSTART(DEM_BCMO) )
      allocate( DEM_BCMO_IJKEND(DEM_BCMO) )

      DEM_BCMO_IJKSTART = -1
      DEM_BCMO_IJKEND   = -1

      dFlag = (DMP_LOG .AND. setDBG)
      if(dFlag) write(*,"(2/,2x,'DEM outlet count: ',I4)") DEM_BCMO

! Loop over the outflow BCs to get an approximate count of the number
! of fluid cells that are adjacent to the outlet.
      MAX_CELLS = 0
      DO BCV_I = 1, DEM_BCMO
         BCV = DEM_BCMO_MAP(BCV_I)

! Set the search area to the dimensions of the inlet.
         if(dFlag) WRITE(*,"(/2x,'Adding cells for BCV: ',I3)") BCV
         SELECT CASE (BC_PLANE(BCV))
         CASE('N','S')
            BND1 = IofPOS(BC_X_e(BCV)) - IofPOS(BC_X_w(BCV))
            BND2 = KofPOS(BC_Z_t(BCV)) - KofPOS(BC_Z_b(BCV))

         CASE('E','W')
            BND1 = JofPOS(BC_Y_n(BCV)) - JofPOS(BC_Y_s(BCV))
            BND2 = KofPOS(BC_Z_t(BCV)) - KofPOS(BC_Z_b(BCV))

         CASE('T','B')
            BND1 = IofPOS(BC_X_e(BCV)) - IofPOS(BC_X_w(BCV))
            BND2 = JofPOS(BC_Y_n(BCV)) - JofPOS(BC_Y_s(BCV))
         CASE DEFAULT
            STOP __LINE__
         END SELECT

         MAX_CELLS = MAX_CELLS +                                      &
            2*(BND1+1)*(BND2+1) + 2*(BND1+2) + 2*(BND2+2)

         if(dFlag) WRITE(*,"(4x,'Plane:   ',A)") BC_PLANE(BCV)
         if(dFlag) WRITE(*,"(4x,'Cells: ', I8)") (BND1+1)*(BND2+1)
      ENDDO

      if(dFlag) write(*,"(2x,'Max Cells: ',I8)") MAX_CELLS

! Allocate an array to hold the IJK values. This array should be
! more than enough to store all the IJKs.
      allocate( LOC_DEM_BCMO_IJK(MAX_CELLS) )

! Loop over the IJKs for each BC and store only the IJKs that you
! own as well as the start/end array positions for each BC.
      LC = 1
      DO BCV_I = 1, DEM_BCMO

         DEM_BCMO_IJKSTART(BCV_I) = LC
         BCV = DEM_BCMO_MAP(BCV_I)

         if(dFlag) write(*,"(/2x,'Searching for fluid cells:',I3)") BCV

         I_w = IofPOS(BC_X_w(BCV)); I_e = IofPOS(BC_X_e(BCV))
         J_s = JofPOS(BC_Y_s(BCV)); J_n = JofPOS(BC_Y_n(BCV))
         IF(DO_K) THEN
            K_b = KofPOS(BC_Z_b(BCV)); K_t = KofPOS(BC_Z_t(BCV))
         ELSE
            K_b = 1; K_t = 1
         ENDIF

! Depending on the flow plane, the 'common' index needs shifted to
! reference the fluid cell.
         SELECT CASE (BC_PLANE(BCV))
         CASE('N'); J_s = J_s+1;  J_n = J_s
         CASE('S'); J_s = J_s-1;  J_n = J_s
         CASE('E'); I_w = I_w+1;  I_e = I_w
         CASE('W'); I_w = I_w-1;  I_e = I_w
         CASE('T'); K_b = K_b+1;  K_t = K_b
         CASE('B'); K_b = K_b-1;  K_t = K_b
        END SELECT

         if(dFlag) then
            write(*,"(4x,'Search bounds: ')")
            write(*,"(6x,'I_w/I_e:',2(2x,I6))") I_w, I_e
            write(*,"(6x,'J_s/J_n:',2(2x,I6))") J_s, J_n
            write(*,"(6x,'K_b/K_t:',2(2x,I6))") K_b, K_t
         endif

! Store the IJKs.
         DO K = K_b, K_t
         DO J = J_s, J_n
         DO I = I_w, I_e
! Skip cells that this rank does not own or are considered dead.
            IF(.NOT.dg_is_ON_myPE_plus1layers(I,J,K))CYCLE

            IJK = DG_FUNIJK(I,J,K)
            LOC_DEM_BCMO_IJK(LC) = IJK
            LC = LC+1
         ENDDO
         ENDDO
         ENDDO

         if(dFlag) write(*,"(/2x,'Adding boundary cells:',I3)") BCV

         I_w = IofPOS(BC_X_w(BCV))-1; I_e = IofPOS(BC_X_e(BCV))+1
         J_s = JofPOS(BC_Y_s(BCV))-1; J_n = JofPOS(BC_Y_n(BCV))+1

         IF(DO_K) THEN
            K_b = KofPOS(BC_Z_b(BCV))-1; K_t = KofPOS(BC_Z_t(BCV))+1
         ELSE
            K_b = 1;   K_t = 1
         ENDIF

! Depending on the flow plane, the 'common' index needs shifted to
! reference the fluid cell.
         SELECT CASE (BC_PLANE(BCV))
         CASE('N','S'); J_s = J_s+1;  J_n = J_n-1
         CASE('E','W'); I_w = I_w+1;  I_e = I_e-1
         CASE('T','B'); K_b = K_b+1;  K_t = K_t-1
         END SELECT

         if(dFlag) then
            write(*,"(4x,'Search bounds: ')")
            write(*,"(6x,'I_w/I_e:',2(2x,I6))") I_w, I_e
            write(*,"(6x,'J_s/J_n:',2(2x,I6))") J_s, J_n
            write(*,"(6x,'K_b/K_t:',2(2x,I6))") K_b, K_t
         endif

! Store the IJKs.
         DO K = K_b, K_t
         DO J = J_s, J_n
         DO I = I_w, I_e
! Skip cells that this rank does not own or are considered dead.
            IF(.NOT.dg_is_ON_myPE_plus1layers(I,J,K))CYCLE

            IJK = DG_FUNIJK(I,J,K)
            LOC_DEM_BCMO_IJK(LC) = IJK
            LC = LC+1
         ENDDO
         ENDDO
         ENDDO

         DEM_BCMO_IJKEND(BCV_I) = LC-1

         if(dFLAG) write(*,1111) BCV, BCV_I,                           &
            DEM_BCMO_IJKSTART(BCV_I),DEM_BCMO_IJKEND(BCV_I)

      ENDDO

 1111 FORMAT(/2x,'DEM Mass Outflow:',/4x,'BC:',I4,3x,'MAP:',I4,&
         /4x,'IJKSTART:',I6,/4x,'IJKEND:  ',I6)

! Allocate the global store arrary array. This changes across MPI ranks.
      IF(LC > 1) THEN
         allocate( DEM_BCMO_IJK(LC-1) )
         DEM_BCMO_IJK(1:LC-1) = LOC_DEM_BCMO_IJK(1:LC-1)
      ELSE
         allocate( DEM_BCMO_IJK(1) )
         DEM_BCMO_IJK(1) = LOC_DEM_BCMO_IJK(1)
      ENDIF

      deallocate(LOC_DEM_BCMO_IJK)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_BC_DEM_MO
