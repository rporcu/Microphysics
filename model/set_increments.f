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
      USE cutcell, ONLY: CARTESIAN_GRID
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
! Used for checking iteration over core cells
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: ALREADY_VISITED
      INTEGER :: interval, j_start(2), j_end(2)
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

! Loop over all cells
      DO K = KSTART3, KEND3
      DO J = JSTART3, JEND3
      DO I = ISTART3, IEND3

         IJK = FUNIJK(I,J,K)  ! Find value of IJK

         I_OF(IJK) = I
         J_OF(IJK) = J
         K_OF(IJK) = K

      ENDDO
      ENDDO
      ENDDO

      ICLASS = 0

! Loop over all cells (minus the ghost layers)
      DO K = KSTART3, KEND3
      DO J = JSTART3, JEND3
      L100: DO I = ISTART3, IEND3

         IJK = FUNIJK(I,J,K)

! Find the the effective cell-center indices for all neighbor cells
         CALL SET_INDEX1A (I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, &
            IJKP, IJKW, IJKE, IJKS, IJKN, IJKB, IJKT)

         ICLASS = ICLASS + 1               !Increment the ICLASS counter
         IF(ICLASS > MAX_CLASS) THEN
            WRITE(ERR_MSG, 1200) trim(iVal(MAX_CLASS))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1200 FORMAT('Error 1200: The number of classes has exceeded the ',    &
         'maximum: ',A,/'Increase the MAX_CLASS parameter in param1',  &
         '_mod.f and recompile.')

         INCREMENT_FOR_N(ICLASS)  = IJKN - IJK
         INCREMENT_FOR_S(ICLASS)  = IJKS - IJK
         INCREMENT_FOR_E(ICLASS)  = IJKE - IJK
         INCREMENT_FOR_W(ICLASS)  = IJKW - IJK
         INCREMENT_FOR_T(ICLASS)  = IJKT - IJK
         INCREMENT_FOR_B(ICLASS)  = IJKB - IJK

         INCREMENT_FOR_IM(ICLASS) = IMJK - IJK
         INCREMENT_FOR_IP(ICLASS) = IPJK - IJK
         INCREMENT_FOR_JM(ICLASS) = IJMK - IJK
         INCREMENT_FOR_JP(ICLASS) = IJPK - IJK
         INCREMENT_FOR_KM(ICLASS) = IJKM - IJK
         INCREMENT_FOR_KP(ICLASS) = IJKP - IJK

         INCREMENT_FOR_NB(1,ICLASS) = INCREMENT_FOR_E(ICLASS)
         INCREMENT_FOR_NB(2,ICLASS) = INCREMENT_FOR_W(ICLASS)
         INCREMENT_FOR_NB(3,ICLASS) = INCREMENT_FOR_S(ICLASS)
         INCREMENT_FOR_NB(4,ICLASS) = INCREMENT_FOR_N(ICLASS)
         INCREMENT_FOR_NB(5,ICLASS) = INCREMENT_FOR_B(ICLASS)
         INCREMENT_FOR_NB(6,ICLASS) = INCREMENT_FOR_T(ICLASS)

         INCREMENT_FOR_MP(1,ICLASS) = INCREMENT_FOR_IM(ICLASS)
         INCREMENT_FOR_MP(2,ICLASS) = INCREMENT_FOR_IP(ICLASS)
         INCREMENT_FOR_MP(3,ICLASS) = INCREMENT_FOR_JM(ICLASS)
         INCREMENT_FOR_MP(4,ICLASS) = INCREMENT_FOR_JP(ICLASS)
         INCREMENT_FOR_MP(5,ICLASS) = INCREMENT_FOR_KM(ICLASS)
         INCREMENT_FOR_MP(6,ICLASS) = INCREMENT_FOR_KP(ICLASS)


         DENOTE_CLASS(ICLASS) = INCREMENT_FOR_N(ICLASS) + INCREMENT_FOR_S&
            (ICLASS) + INCREMENT_FOR_E(ICLASS) + INCREMENT_FOR_W(ICLASS)&
             + INCREMENT_FOR_T(ICLASS) + INCREMENT_FOR_B(ICLASS) + &
            INCREMENT_FOR_IM(ICLASS) + INCREMENT_FOR_IP(ICLASS) + &
            INCREMENT_FOR_JM(ICLASS) + INCREMENT_FOR_JP(ICLASS) + &
            INCREMENT_FOR_KM(ICLASS) + INCREMENT_FOR_KP(ICLASS)

         CELL_CLASS(IJK) = ICLASS

! Place the cell in a class based on its DENOTE_CLASS(ICLASS) value
         DO IC = 1, ICLASS - 1             !Loop over previous and present classes
!                                                !IF a possible match in cell types
            IF(DENOTE_CLASS(ICLASS) == DENOTE_CLASS(IC)) THEN
!                                                !is found, compare all increments
               IF(INCREMENT_FOR_N(ICLASS) /= INCREMENT_FOR_N(IC)) CYCLE
               IF(INCREMENT_FOR_S(ICLASS) /= INCREMENT_FOR_S(IC)) CYCLE
               IF(INCREMENT_FOR_E(ICLASS) /= INCREMENT_FOR_E(IC)) CYCLE
               IF(INCREMENT_FOR_W(ICLASS) /= INCREMENT_FOR_W(IC)) CYCLE
               IF(INCREMENT_FOR_T(ICLASS) /= INCREMENT_FOR_T(IC)) CYCLE
               IF(INCREMENT_FOR_B(ICLASS) /= INCREMENT_FOR_B(IC)) CYCLE
               IF(INCREMENT_FOR_IM(ICLASS) /= INCREMENT_FOR_IM(IC)) CYCLE
               IF(INCREMENT_FOR_IP(ICLASS) /= INCREMENT_FOR_IP(IC)) CYCLE
               IF(INCREMENT_FOR_JM(ICLASS) /= INCREMENT_FOR_JM(IC)) CYCLE
               IF(INCREMENT_FOR_JP(ICLASS) /= INCREMENT_FOR_JP(IC)) CYCLE
               IF(INCREMENT_FOR_KM(ICLASS) /= INCREMENT_FOR_KM(IC)) CYCLE
               IF(INCREMENT_FOR_KP(ICLASS) /= INCREMENT_FOR_KP(IC)) CYCLE
               CELL_CLASS(IJK) = IC        !Assign cell to a class
               ICLASS = ICLASS - 1
               CYCLE  L100                 !Go to next cell
            ENDIF
         END DO

      ENDDO L100
      ENDDO
      ENDDO

      DO M = 1, MMAX
      DO L = M, MMAX
         IF(L == M) THEN
            STORE_LM(L,M) = 0
         ELSE
            STORE_LM(L,M) = M + (L - 2)*(L - 1)/2
            STORE_LM(M,L) = M + (L - 2)*(L - 1)/2
         ENDIF
      ENDDO
      ENDDO

      USE_CORECELL_LOOP = .not.CARTESIAN_GRID

      if (USE_CORECELL_LOOP) then
         Allocate( already_visited(DIMENSION_3))
         already_visited(:) = .false.

         core_istart = istart+2
         core_iend = iend-2

         core_jstart = jstart+2
         core_jend = jend-2

         if (do_k) then
            core_kstart = kstart+2
            core_kend = kend-2
         else
            core_kstart = 1
            core_kend = 1
            kstart = 1
            kend = 1
         endif

         iclass = cell_class(funijk(core_istart,core_jstart,core_kstart))

         outer: do k = core_kstart,core_kend
            do i = core_istart,core_iend
               do j = core_jstart,core_jend
                  IJK = funijk(i,j,k)
                  ! this shouldn't happen, but we might as well check
                  if (ijk.ne. (j + c0 + i*c1 + k*c2)) then
                     USE_CORECELL_LOOP = .false.
                     exit outer
                  endif

                  ijk = (j + c0 + i*c1 + k*c2)

                  if (already_visited(ijk)) then
                     USE_CORECELL_LOOP = .false.
                     exit outer
                  endif
                  already_visited(ijk) = .true.

                  if (iclass.ne.cell_class(ijk)) then
                     USE_CORECELL_LOOP = .false.
                     exit outer
                  endif
               enddo
            enddo
         enddo outer

         j_start(1) = jstart
         j_end(1) = jend
         j_start(2) = 0 ! no iterations
         j_end(2) = -1  ! no iterations

         outer2: do k = kstart,kend
            do i = istart,iend

               if  (USE_CORECELL_LOOP) then
                  if (core_istart<= i .and. i <= core_iend .and. core_kstart <= k .and. k<=core_kend) then
                     j_start(1) = jstart
                     j_end(1) = core_jstart-1
                     j_start(2) = core_jend+1
                     j_end(2) = jend
                  else
                     j_start(1) = jstart
                     j_end(1) = jend
                     j_start(2) = 0 ! no iterations
                     j_end(2) = -1  ! no iterations
                  endif
               endif

               do interval=1,2
                  do j = j_start(interval),j_end(interval)
                     if (already_visited(funijk(i,j,k))) then
                        USE_CORECELL_LOOP = .false.
                        exit outer2
                     endif
                     already_visited(funijk(i,j,k)) = .true.
                  enddo
               enddo
            enddo
         enddo outer2

         outer3: do k = kstart,kend
            do i = istart,iend
               do j = jstart,jend
                  if (.not.already_visited(funijk(i,j,k))) then
                     USE_CORECELL_LOOP = .false.
                     exit outer3
                  endif
               enddo
            enddo
         enddo outer3

         deallocate(already_visited)

      endif

      IF(.NOT.INCREMENT_ARRAYS_ALLOCATED) THEN
         allocate(WEST_ARRAY_OF(ijkstart3:ijkend3))
         allocate(EAST_ARRAY_OF(ijkstart3:ijkend3))
         allocate(SOUTH_ARRAY_OF(ijkstart3:ijkend3))
         allocate(NORTH_ARRAY_OF(ijkstart3:ijkend3))
         allocate(BOTTOM_ARRAY_OF(ijkstart3:ijkend3))
         allocate(TOP_ARRAY_OF(ijkstart3:ijkend3))

         allocate(IM_ARRAY_OF(ijkstart3:ijkend3))
         allocate(IP_ARRAY_OF(ijkstart3:ijkend3))
         allocate(JM_ARRAY_OF(ijkstart3:ijkend3))
         allocate(JP_ARRAY_OF(ijkstart3:ijkend3))
         allocate(KM_ARRAY_OF(ijkstart3:ijkend3))
         allocate(KP_ARRAY_OF(ijkstart3:ijkend3))
      ENDIF

      INCREMENT_ARRAYS_ALLOCATED = .TRUE.

      DO IJK = ijkstart3,ijkend3
         WEST_ARRAY_OF(ijk)   = WEST_OF_0(IJK)
         EAST_ARRAY_OF(ijk)   = EAST_OF_0(IJK)
         SOUTH_ARRAY_OF(ijk)  = SOUTH_OF_0(IJK)
         NORTH_ARRAY_OF(ijk)  = NORTH_OF_0(IJK)
         BOTTOM_ARRAY_OF(ijk) = BOTTOM_OF_0(IJK)
         TOP_ARRAY_OF(ijk)    = TOP_OF_0(IJK)

         IM_ARRAY_OF(ijk) = IM_OF_0(IJK)
         IP_ARRAY_OF(ijk) = IP_OF_0(IJK)
         JM_ARRAY_OF(ijk) = JM_OF_0(IJK)
         JP_ARRAY_OF(ijk) = JP_OF_0(IJK)
         KM_ARRAY_OF(ijk) = KM_OF_0(IJK)
         KP_ARRAY_OF(ijk) = KP_OF_0(IJK)
      ENDDO

      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE SET_INCREMENTS
