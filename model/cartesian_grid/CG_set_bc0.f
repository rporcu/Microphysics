!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_SET_BC0                                             C
!  Purpose: This module does the initial setting of boundary           C
!           conditions for cut cells only                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CG_SET_BC0
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE bc
      USE boundfunijk
      USE compar
      USE cutcell
      USE eos, ONLY: EOSG
      USE fldvar
      USE functions
      USE funits
      USE geometry
      USE indices
      USE mpi_utility
      USE param
      USE param1
      USE physprop
      USE quadric
      USE run
      USE scales
      USE sendrecv
      USE toleranc

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

!
!                      Local index for boundary condition
      INTEGER          L
!
!                      indices
      INTEGER          IJK, M, N ,IJKW,IJKS,IJKB
!
!----------------------------------------------

      INTEGER, DIMENSION(8) :: ACCEPTABLE_DEFAULT_WALL=-1
      LOGICAL :: GLOBAL_CORNER

!
!  Define global corners as acceptable default walls
!  These cells should never be used
!

      IF(.NOT.RE_INDEXING.AND.NumPEs==1) THEN

      ACCEPTABLE_DEFAULT_WALL(1) = FUNIJK(IMIN3,JMIN3,KMIN3)
      ACCEPTABLE_DEFAULT_WALL(2) = FUNIJK(IMAX3,JMIN3,KMIN3)
      ACCEPTABLE_DEFAULT_WALL(3) = FUNIJK(IMIN3,JMAX3,KMIN3)
      ACCEPTABLE_DEFAULT_WALL(4) = FUNIJK(IMAX3,JMAX3,KMIN3)
      ACCEPTABLE_DEFAULT_WALL(5) = FUNIJK(IMIN3,JMIN3,KMAX3)
      ACCEPTABLE_DEFAULT_WALL(6) = FUNIJK(IMAX3,JMIN3,KMAX3)
      ACCEPTABLE_DEFAULT_WALL(7) = FUNIJK(IMIN3,JMAX3,KMAX3)
      ACCEPTABLE_DEFAULT_WALL(8) = FUNIJK(IMAX3,JMAX3,KMAX3)

      ENDIF

!      DO N = 1,8
!         print*,'acceptable default=',ACCEPTABLE_DEFAULT_WALL(N)
!      ENDDO


      DO IJK = ijkstart3, ijkend3

         L = BC_ID(IJK)

         IF(L>0) THEN
            IF(BC_TYPE(L)=='CG_PO') THEN

               P_G(IJK) = SCALE_PRESSURE(BC_P_G(L))
               IF (BC_EP_G(L) /= UNDEFINED) EP_G(IJK) = BC_EP_G(L)

            ELSEIF(BC_TYPE(L)=='CG_MI') THEN

               EP_G(IJK) = BC_EP_G(L)
               P_G(IJK) = SCALE_PRESSURE(BC_P_G(L))

               IF(BC_U_g(L)/=UNDEFINED) THEN
                  U_G(IJK) =  BC_U_g(L)
               ELSE
                  U_G(IJK) =  BC_VELMAG_g(L)*NORMAL_S(IJK,1)
               ENDIF

               IF(BC_V_g(L)/=UNDEFINED) THEN
                  V_G(IJK) =  BC_V_g(L)
               ELSE
                  V_G(IJK) =  BC_VELMAG_g(L)*NORMAL_S(IJK,2)
               ENDIF

               IF(BC_W_g(L)/=UNDEFINED) THEN
                  W_G(IJK) =  BC_W_g(L)
               ELSE
                  W_G(IJK) =  BC_VELMAG_g(L)*NORMAL_S(IJK,3)
               ENDIF

               IJKW = WEST_OF(IJK)
               IJKS = SOUTH_OF(IJK)
               IJKB = BOTTOM_OF(IJK)

               IF(FLUID_AT(IJKW)) THEN
                  IF(BC_U_g(L)/=UNDEFINED) THEN
                     U_G(IJKW) =  BC_U_g(L)
                  ELSE
                     U_G(IJKW) =  BC_VELMAG_g(L)*NORMAL_S(IJK,1)
                  ENDIF
               ENDIF

               IF(FLUID_AT(IJKS)) THEN
                  IF(BC_V_g(L)/=UNDEFINED) THEN
                     V_G(IJKS) =  BC_V_g(L)
                  ELSE
                     V_G(IJKS) =  BC_VELMAG_g(L)*NORMAL_S(IJK,2)
                  ENDIF
               ENDIF

               IF(FLUID_AT(IJKB)) THEN
                  IF(BC_W_g(L)/=UNDEFINED) THEN
                     W_G(IJKB) =  BC_W_g(L)
                  ELSE
                     W_G(IJKB) =  BC_VELMAG_g(L)*NORMAL_S(IJK,3)
                  ENDIF
               ENDIF

               IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_AVG,&
                  P_G(IJK),295.15d0)
               ROP_G(IJK) = EP_G(IJK)*RO_G(IJK)

            ENDIF
         ENDIF

         IF(DEFAULT_WALL_AT(IJK)) THEN

!            print*,'Default_wall_at IJK=',IJK,I_OF(IJK),J_OF(IJK),K_OF(IJK)

            GLOBAL_CORNER = .FALSE.
            DO N = 1,8
               IF(IJK==ACCEPTABLE_DEFAULT_WALL(N)) GLOBAL_CORNER = .TRUE.
            ENDDO

            IF(.NOT.GLOBAL_CORNER.AND..NOT.BLOCKED_CELL_AT(IJK)) THEN

               ICBC_FLAG(IJK)(2:3) = 'CG'

               IF((MyPE == PE_IO).AND.PRINT_WARNINGS) THEN
                  WRITE(*,*) 'WARNING: DEFAULT WALL DETECTED AT I,J,K = ',I_OF(IJK),J_OF(IJK),K_OF(IJK) ,BLOCKED_CELL_AT(IJK)
                  WRITE(*,*) '         WHEN USING CARTESIAN GRID CUT-CELL FEATURE.'
                  WRITE(*,*) '         DEFAULT WALLS ARE NOT ALLOWED WITH CUT-CELLS.'
                  WRITE(*,*) '         THE DEFAULT WALL WAS REMOVED ALONG THIS CELL.'
                  WRITE(*,*) ''
               ENDIF
!               CALL MFIX_EXIT(MYPE)

            ENDIF

         ENDIF

      ENDDO

      RETURN
      END SUBROUTINE CG_SET_BC0

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 020 New local variables for parallelization: FLAG_G , FLUID_AT_G
!// 360 Check if i,j,k resides on current processor
