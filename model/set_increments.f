!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: SET_INCREMENTS                                         !
!  Author: M. Syamlal, W. Rogers                      Date: 10-DEC-91  !
!                                                                      !
!  Purpose: The purpose of this module is to create increments to be   !
!           stored in the array STORE_INCREMENT which will be added    !
!           to cell index ijk to find the effective indices of its     !
!           neighbors. These increments are found using the 'class'    !
!           of cell ijk. The class is determined based on the          !
!           neighboring cell type, i.e. wall or fluid.                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_INCREMENTS

      USE compar
      USE fldvar
      USE functions
      USE funits
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop

! Module procedures
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM
      use error_manager

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IJK
      INTEGER :: IMJK, IPJK, IJKW, IJKE  ! I+, I-, east/west
      INTEGER :: IJMK, IJPK, IJKS, IJKN  ! J+, J-, north/south
      INTEGER :: IJKM, IJKP, IJKB, IJKT  ! K+, K-, top/bottom
! DO-loop index, ranges from 1 to ICLASS
      INTEGER :: IC
! Index for the solids phase.
      INTEGER :: M
! Local DO-loop index
      INTEGER :: L
! Index denoting cell class
      INTEGER :: ICLASS
! Array of sum of increments to make the class determination faster.
      INTEGER :: DENOTE_CLASS(MAX_CLASS)
! Flags for using the 'real' I/J/K value (not cyclic.)
      LOGICAL :: SHIFT
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("SET_INCREMENTS")

! Allocate increment arrays and report an allocation errors.
      CALL ALLOCATE_ARRAYS_INCREMENTS

! Initialize the default values to Undefined_I
      IP1(:) = UNDEFINED_I
      IM1(:) = UNDEFINED_I
      JP1(:) = UNDEFINED_I
      JM1(:) = UNDEFINED_I
      KP1(:) = UNDEFINED_I
      KM1(:) = UNDEFINED_I

      DO I = ISTART3, IEND3
         SHIFT = .NOT.(I==IMIN3 .OR. I==IMIN2 .OR. &
                       I==IMAX3 .OR. I==IMAX2)

         IF(CYCLIC_X .AND. NODESI.EQ.1 .AND. DO_I .AND. SHIFT) THEN
            IP1(I) = IMAP_C(IMAP_C(I)+1)
            IM1(I) = IMAP_C(IMAP_C(I)-1)
         ELSE
            IM1(I) = MAX(ISTART3, I - 1)
            IP1(I) = MIN(IEND3,   I + 1)
         ENDIF
      ENDDO

      DO J = JSTART3, JEND3

         SHIFT = .NOT.(J==JMIN3 .OR. J==JMIN2 .OR. &
                       J==JMAX3 .OR. J==JMAX2)

         IF (CYCLIC_Y .AND. NODESJ.EQ.1 .AND. DO_J .AND. SHIFT) THEN
            JP1(J) = JMAP_C(JMAP_C(J)+1)
            JM1(J) = JMAP_C(JMAP_C(J)-1)
         ELSE
            JM1(J) = MAX(JSTART3,J - 1)
            JP1(J) = MIN(JEND3,  J + 1)
         ENDIF
      ENDDO


      DO K = KSTART3, KEND3

         SHIFT = .NOT.(K==KMIN3 .OR. K==KMIN2 .OR. &
                       K==KMAX3 .OR. K==KMAX2)

         IF(CYCLIC_Z .AND. NODESK.EQ.1 .AND. DO_K .AND. SHIFT) THEN
            KP1(K) = KMAP_C(KMAP_C(K)+1)
            KM1(K) = KMAP_C(KMAP_C(K)-1)
         ELSE
            KM1(K) = MAX(KSTART3,K - 1)
            KP1(K) = MIN(KEND3,K + 1)
         ENDIF
      ENDDO

      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE SET_INCREMENTS
