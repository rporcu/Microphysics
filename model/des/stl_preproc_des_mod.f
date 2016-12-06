!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: stl_preproc_des                                        !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose: This module containd routines for geometric interaction    !
!  required for STL files.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE STL_PREPROC_DES

      IMPLICIT NONE

! Use this module only to define functions and subroutines.
      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DES_STL_PREPROCESSING                                   !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_STL_PREPROCESSING

! Number of facets from STL files, (plus DES generated)
      use stl, only: N_FACETS_DES
! Start/End position of different STLs
      use stl, only: STL_START, STL_END
! All STLS
      use stl, only: ALL_STL
! STLs read from geometry files
      use stl, only: BASE_STL
! STLs for user specified walls (NSW, PSW, FSW)
      use stl, only: BCWALLS_STL
! STLs for impermeable surfaces
      use stl, only: IMPRMBL_STL
! STLs for default walls
      use stl, only: DEFAULT_STL

      use stl_dbg_des, only: stl_dbg_write_facets

      use error_manager, only: err_msg, flush_err_msg

      IMPLICIT NONE

! Pre-procssing for the des in order to assign facets to grid cells.
      WRITE(ERR_MSG,"('Pre-Processing geometry for DES.')")
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)


! Process the STL files
      N_FACETS_DES = 0
! Store the Start/End of the base STLs from geometry files
      STL_START(BASE_STL)=1;   STL_END(BASE_STL)=N_FACETS_DES

! Process stair-step geometries
      CALL CONVERT_BC_WALLS_TO_STL

! Bin the STL to the DES grid.
      CALL BIN_FACETS_TO_DG

! Some functions for debugging.
!      CALL STL_DBG_WRITE_FACETS(BASE_STL)
!      CALL STL_DBG_WRITE_FACETS(BCWALLS_STL)
!      CALL STL_DBG_WRITE_FACETS(IMPRMBL_STL)
!      CALL STL_DBG_WRITE_FACETS(DEFAULT_STL)
      CALL STL_DBG_WRITE_FACETS(ALL_STL)
!      CALL STL_DBG_WRITE_STL_FROM_DG(STL_TYPE=BASE_STL)

! Pre-procssing for the des in order to assign facets to grid cells.
      WRITE(ERR_MSG,"('DES geometry pre-processing complete.')")
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      RETURN
      END SUBROUTINE DES_STL_PREPROCESSING


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: BIN_FACETS_TO_DG                                        !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine BIN_FACETS_TO_DG

      use desgrid, only: DG_IJKSIZE2
      use desgrid, only: DG_IEND2, DG_ISTART2
      use desgrid, only: DG_JEND2, DG_JSTART2
      use desgrid, only: DG_KEND2, DG_KSTART2

      use stl, only: FACETS_AT_DG

      use geometry, only: XLENGTH, YLENGTH, ZLENGTH
      use stl, only: N_FACETS_DES
      use stl, only: VERTEX

      use stl, only: TOL_STL
      use param1, only: ZERO, ONE

      use desgrid, only: DG_FUNIJK
      use desgrid, only: IofPOS, JofPOS, KofPOS
      use desgrid, only: dg_is_ON_myPE_plus1layers

      IMPLICIT NONE

! DES Grid cell index.
      INTEGER :: IJK
! Loop counters:
      INTEGER :: I1, I2, II  ! X-axis
      INTEGER :: J1, J2, JJ  ! Y-axis
      INTEGER :: K1, K2, KK  ! Z-axis
      INTEGER :: NN          ! STLs

! Maximum and minimum extents of the indexed STL
      DOUBLE PRECISION:: X1,Y1,Z1
      DOUBLE PRECISION:: X2,Y2,Z2

! Allocate the data storage array.
      IF(.not.allocated(FACETS_AT_DG)) &
         allocate(FACETS_AT_DG(DG_IJKSIZE2))

      FACETS_AT_DG(:)%COUNT = 0

      DO NN = 1,N_FACETS_DES

         X1 = minval(VERTEX(1:3,1,NN))
         X2 = maxval(VERTEX(1:3,1,NN))
         Y1 = minval(VERTEX(1:3,2,NN))
         Y2 = maxval(VERTEX(1:3,2,NN))
         Z1 = minval(VERTEX(1:3,3,NN))
         Z2 = maxval(VERTEX(1:3,3,NN))

         I1 = DG_IEND2
         I2 = DG_ISTART2
         IF(X2>=-TOL_STL .AND. X1<=XLENGTH+TOL_STL) THEN
            I1 = max(iofpos(X1)-1, dg_istart2)
            I2 = min(iofpos(X2)+1, dg_iend2)
         ENDIF

         J1 = DG_JEND2
         J2 = DG_JSTART2
         IF(Y2>=-TOL_STL .AND. Y1<=YLENGTH+TOL_STL) THEN
            J1 = max(jofpos(Y1)-1, dg_jstart2)
            J2 = min(jofpos(Y2)+1, dg_jend2)
         ENDIF

         K1 = DG_KEND2
         K2 = DG_KSTART2
         IF(Z2>=-TOL_STL .AND. Z1<=ZLENGTH+TOL_STL) THEN
            K1 = max(kofpos(Z1)-1, dg_kstart2)
            K2 = min(kofpos(Z2)+1, dg_kend2)
         ENDIF

         DO KK=K1,K2
         DO JJ=J1,J2
         DO II=I1,I2
            IF(dg_is_ON_myPE_plus1layers(II,JJ,KK)) THEN
               IJK = DG_FUNIJK(II,JJ,KK)
               CALL ADD_FACET_FOR_DES(II,JJ,KK,IJK,NN)
            ENDIF
         ENDDO
         ENDDO
         ENDDO

      ENDDO

      RETURN
      END SUBROUTINE BIN_FACETS_TO_DG




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ADD_FACET_FOR_DES                                       !
!  Author: Rahul Garg                                  Date: 24-Oct-13 !
!                                                                      !
!  Purpose: Add facets to DES grid cells.                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ADD_FACET_FOR_DES(I,J,K,IJK,N)

      use desgrid, only: dg_dxinv, dg_xstart, dg_istart1
      use desgrid, only: dg_dyinv, dg_ystart, dg_jstart1
      use desgrid, only: dg_dzinv, dg_zstart, dg_kstart1

      use discretelement, only: MAX_RADIUS

      use stl, only: VERTEX

      use stl_functions_des, only: TRI_BOX_OVERLAP

      use param1, only: ZERO, HALF, ONE
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ival

      IMPLICIT NONE

! DES grid index and facet index
      INTEGER, INTENT(IN) :: I,J,K,IJK, N

! Center of DES grid cell and half size. Note that a buffer is added to
! the half size to make the cell appear a little larger. This ensures
! that paricles near the edge 'see' STLs that are nearby but do not
! directly intersect the DES grid cell contain the particle center.
      DOUBLE PRECISION :: CENTER(3), HALFSIZE(3)
! Flag: STL intersects the DES grid cell
      LOGICAL :: OVERLAP
! DES grid cell dimensions
      DOUBLE PRECISION :: lDX, lDY, lDZ
! Buffer to ensure all particle-STL collisions are captured.
      DOUBLE PRECISION :: BUFFER

      BUFFER = 1.1d0*MAX_RADIUS

      lDX = ONE/DG_DXINV
      lDY = ONE/DG_DYINV
      lDZ = ONE/DG_DZINV

      CENTER(1) = dg_xstart + (dble(I-dg_istart1)+HALF)*lDX
      HALFSIZE(1) = HALF*lDX + BUFFER

      CENTER(2) = dg_ystart + (dble(J-dg_jstart1)+HALF)*lDY
      HALFSIZE(2) = HALF*lDY + BUFFER
      CENTER(3) = dg_zstart + (dble(K-dg_kstart1)+HALF)*lDZ
      HALFSIZE(3) = HALF*lDZ + BUFFER

      CALL TRI_BOX_OVERLAP(CENTER, HALFSIZE, VERTEX(:,:,N), OVERLAP)

      IF(OVERLAP) CALL ADD_FACET(IJK, N)

      RETURN
      END SUBROUTINE ADD_FACET_FOR_DES



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ADD_FACET                                               !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ADD_FACET(IJK, FACET_ID)

      use stl, only: VERTEX
      use stl, only: FACETS_AT_DG
      use param1, only: ZERO

      implicit none

      INTEGER, INTENT(IN) :: IJK, facet_id

      INTEGER, ALLOCATABLE :: int_tmp(:)
      DOUBLE PRECISION, ALLOCATABLE :: real_tmp(:)

      INTEGER :: lSIZE, II,FC
      DOUBLE PRECISION :: smallest_extent, min_temp, max_temp


      FC = FACETS_AT_DG(IJK)%COUNT
      IF(FC > 0) THEN
!      IF(FACETS_AT_DG(IJK)%COUNT > 0) THEN

         DO II=1, FACETS_AT_DG(IJK)%COUNT
            IF(FACET_ID == FACETS_AT_DG(IJK)%ID(II)) RETURN
         ENDDO

         FACETS_AT_DG(IJK)%COUNT = FACETS_AT_DG(IJK)%COUNT+1

         lSIZE = size(FACETS_AT_DG(IJK)%ID)
         IF(FACETS_AT_DG(IJK)%COUNT +1> lSIZE) THEN
            allocate(int_tmp(2*lSIZE)); int_tmp=0
            int_tmp(1:lSIZE) = FACETS_AT_DG(IJK)%ID(1:lSIZE)
            call move_alloc(int_tmp,FACETS_AT_DG(IJK)%ID)

            allocate(int_tmp(2*lSIZE)); int_tmp=0
            int_tmp(1:lSIZE) = FACETS_AT_DG(IJK)%DIR(1:lSIZE)
            call move_alloc(int_tmp, FACETS_AT_DG(IJK)%DIR)

            allocate(real_tmp(2*lSIZE)); real_tmp=ZERO
            real_tmp(1:lSIZE) = FACETS_AT_DG(IJK)%MIN(1:lSIZE)
            call move_alloc(real_tmp, FACETS_AT_DG(IJK)%MIN)

            allocate(real_tmp(2*lSIZE)); real_tmp=ZERO
            real_tmp(1:lSIZE) = FACETS_AT_DG(IJK)%MAX(1:lSIZE)
            call move_alloc(real_tmp, FACETS_AT_DG(IJK)%MAX)
         ENDIF

      ELSE
         FACETS_AT_DG(IJK)%COUNT = 1
         IF(allocated(FACETS_AT_DG(IJK)%ID)) &
            deallocate(FACETS_AT_DG(IJK)%ID)
         allocate(FACETS_AT_DG(IJK)%ID(4))
         IF(allocated(FACETS_AT_DG(IJK)%DIR)) &
            deallocate(FACETS_AT_DG(IJK)%DIR)
         allocate(FACETS_AT_DG(IJK)%DIR(4))
         IF(allocated(FACETS_AT_DG(IJK)%MIN)) &
            deallocate(FACETS_AT_DG(IJK)%MIN)
         allocate(FACETS_AT_DG(IJK)%MIN(4))
         IF(allocated(FACETS_AT_DG(IJK)%MAX)) &
            deallocate(FACETS_AT_DG(IJK)%MAX)
         allocate(FACETS_AT_DG(IJK)%MAX(4))
      ENDIF

      FACETS_AT_DG(IJK)%ID(FACETS_AT_DG(IJK)%COUNT) = FACET_ID

      SMALLEST_EXTENT = HUGE(0.0)

      DO II=1,3
         MIN_TEMP = MINVAL(VERTEX(:,II,FACET_ID))
         MAX_TEMP = MAXVAL(VERTEX(:,II,FACET_ID))
         IF(ABS(MAX_TEMP - MIN_TEMP) < SMALLEST_EXTENT ) THEN
            FACETS_AT_DG(IJK)%DIR(FACETS_AT_DG(IJK)%COUNT) = II
            FACETS_AT_DG(IJK)%MIN(FACETS_AT_DG(IJK)%COUNT) = MIN_TEMP
            FACETS_AT_DG(IJK)%MAX(FACETS_AT_DG(IJK)%COUNT) = MAX_TEMP
            SMALLEST_EXTENT = ABS(MAX_TEMP - MIN_TEMP)
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE ADD_FACET



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CONVERT_BC_WALLS_TO_STL                                 !
!  Author: J.Musser                                   Date: 03-Nov-15  !
!                                                                      !
!  Purpose: Convert user specified walls to STLs for particle-wall     !
!  collision detection.                                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine CONVERT_BC_WALLS_TO_STL

      use bc, only: BC_DEFINED, BC_TYPE
      use bc, only: BC_I_w, BC_I_e
      use bc, only: BC_J_s, BC_J_n
      use bc, only: BC_K_b, BC_K_t
      use bc, only: BC_X_w, BC_X_e
      use bc, only: BC_Y_s, BC_Y_n
      use bc, only: BC_Z_b, BC_Z_t

      use geometry, only: xLength, yLength, zLength
      use stl, only: N_FACETS_DES
      use stl, only: STL_START, STL_END, BCWALLS_STL

      use discretelement, only: XE, YN, ZT

      use param, only: DIMENSION_BC
      USE param1, only: ZERO

      IMPLICIT NONE

! Loop counter.
      INTEGER :: BCV
! Extents of the BC region with respect to the fluid grid.
      DOUBLE PRECISION :: lXw, lXe, lYs, lYn, lZb, lZt

      STL_START(BCWALLS_STL)=N_FACETS_DES+1

      DO BCV=1, DIMENSION_BC
         IF(.NOT.BC_DEFINED(BCV)) CYCLE

         IF(BC_TYPE(BCV) == 'FREE_SLIP_WALL' .OR.   &
            BC_TYPE(BCV) == 'NO_SLIP_WALL'   .OR.   &
            BC_TYPE(BCV) == 'PAR_SLIP_WALL') THEN

            lXw = XE(BC_I_w(BCV)-1)
            lXe = XE(BC_I_e(BCV))
            IF(BC_X_w(BCV) == BC_X_e(BCV)) then
               if(BC_X_w(BCV) == 0.0d0) then
                  lXw = 0.0d0
                  lXe = 0.0d0
               elseif(BC_X_e(BCV) == xLength) then
                  lXw = xLength
                  lXe = xLength
               endif
            endif

            lYs = YN(BC_J_s(BCV)-1)
            lYn =YN(BC_J_n(BCV))
            if(BC_Y_s(bcv) == BC_Y_n(bcv)) then
               if(BC_Y_s(BCV) == 0.0d0) then
                  lYs = 0.0d0
                  lYn = 0.0d0
               elseif(BC_Y_s(BCV) == yLength) then
                  lYs = yLength
                  lYn = yLength
               endif
            endif

            lZt = ZT(BC_K_t(BCV))
            lZb = ZT(BC_K_b(BCV)-1)
            if(BC_Z_b(bcv) == BC_Z_t(bcv)) then
               if(BC_Z_b(bcv) == 0.0d0) then
                  lZb = 0.0d0
                  lZt = 0.0d0
               elseif(BC_Z_b(bcv) == zLength) then
                  lZb = zLength
                  lZt = zLength
               endif
            endif

            CALL GENERATE_STL_BOX(lXw, lXe, lYs, lYn, lZb, lZt)
         ENDIF
      ENDDO
      STL_END(BCWALLS_STL)=N_FACETS_DES

      RETURN
      END SUBROUTINE CONVERT_BC_WALLS_TO_STL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GENERATE_STL_BOX                                        !
!  Author: J.Musser                                   Date: 03-Nov-15  !
!                                                                      !
!  Purpose: Given the six corners of a box, create the 12 STLs needed  !
!  to define the geometry.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GENERATE_STL_BOX(pXw, pXe, pYs, pYn, pZb, pZt)

      use stl, only: VERTEX, NORM_FACE
      use stl, only: N_FACETS_DES
      use geometry, only: xLength, yLength, zLength

      use param1, only: ZERO, ONE

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: pXw, pXe, pYs, pYn, pZb, pZt

! West Face
      if(pYs /= pYn .and. pZb /= pZt) then
         vertex(1,:,n_facets_des + 1) = (/pXw, pYs, pZb/)
         vertex(2,:,n_facets_des + 1) = (/pXw, pYn, pZb/)
         vertex(3,:,n_facets_des + 1) = (/pXw, pYn, pZt/)

         vertex(1,:,n_facets_des + 2) = (/pXw, pYs, pZt/)
         vertex(2,:,n_facets_des + 2) = (/pXw, pYs, pZb/)
         vertex(3,:,n_facets_des + 2) = (/pXw, pYn, pZt/)

         if(pXw == pXe) then
            if(pXw == 0.0d0) then
               norm_face(:,n_facets_des + 1) = (/ one, zero, zero/)
               norm_face(:,n_facets_des + 2) = (/ one, zero, zero/)
            elseif(pXw == xLength) then
               norm_face(:,n_facets_des + 1) = (/-one, zero, zero/)
               norm_face(:,n_facets_des + 2) = (/-one, zero, zero/)
            else
               stop 9800
            endif
         else
            norm_face(:,n_facets_des + 1) = (/-one, zero, zero/)
            norm_face(:,n_facets_des + 2) = (/-one, zero, zero/)
         endif
         n_facets_des = n_facets_des+2

! East Face
         if(pXw /= pXe) then
            vertex(3,:,n_facets_des + 1) = (/pXe, pYs, pZb/)
            vertex(2,:,n_facets_des + 1) = (/pXe, pYn, pZb/)
            vertex(1,:,n_facets_des + 1) = (/pXe, pYn, pZt/)

            vertex(3,:,n_facets_des + 2) = (/pXe, pYs, pZt/)
            vertex(2,:,n_facets_des + 2) = (/pXe, pYs, pZb/)
            vertex(1,:,n_facets_des + 2) = (/pXe, pYn, pZt/)

            norm_face(:,n_facets_des + 1) = (/ one, zero, zero/)
            norm_face(:,n_facets_des + 2) = (/ one, zero, zero/)
            n_facets_des = n_facets_des + 2
         endif
      endif


! South Face
      if(pXw /= pXe .and. pZb /= pZt) then
         vertex(1,:,n_facets_des + 1) = (/pXw, pYs, pZb/)
         vertex(2,:,n_facets_des + 1) = (/pXe, pYs, pZb/)
         vertex(3,:,n_facets_des + 1) = (/pXe, pYs, pZt/)

         vertex(1,:,n_facets_des + 2) = (/pXe, pYs, pZt/)
         vertex(2,:,n_facets_des + 2) = (/pXw, pYs, pZt/)
         vertex(3,:,n_facets_des + 2) = (/pXw, pYs, pZb/)

         if(pYs == pYn) then
            if(pYs == 0.0d0) then
               norm_face(:,n_facets_des + 1) = (/zero, one, zero/)
               norm_face(:,n_facets_des + 2) = (/zero, one, zero/)
            elseif(pYs == yLength) then
               norm_face(:,n_facets_des + 1) = (/zero,-one, zero/)
               norm_face(:,n_facets_des + 2) = (/zero,-one, zero/)
            else
               stop 9801
            endif
            norm_face(:,n_facets_des + 1) = (/zero,-one, zero/)
            norm_face(:,n_facets_des + 2) = (/zero,-one, zero/)
         endif
         n_facets_des = n_facets_des + 2

! North Face
         if(pYs /= pYn) then
            vertex(3,:,n_facets_des + 1) = (/pXw, pYn, pZb/)
            vertex(2,:,n_facets_des + 1) = (/pXe, pYn, pZb/)
            vertex(1,:,n_facets_des + 1) = (/pXe, pYn, pZt/)

            vertex(3,:,n_facets_des + 2) = (/pXe, pYn, pZt/)
            vertex(2,:,n_facets_des + 2) = (/pXw, pYn, pZt/)
            vertex(1,:,n_facets_des + 2) = (/pXw, pYn, pZb/)
            norm_face(:,n_facets_des + 1) = (/zero, one, zero/)
            norm_face(:,n_facets_des + 2) = (/zero, one, zero/)
            n_facets_des = n_facets_des + 2
         endif
      endif

! Bottom Face
      if(pXw /= pXe .and. pYs /= pYn) then
         vertex(1,:,n_facets_des + 1) = (/pXw, pYs, pZb/)
         vertex(2,:,n_facets_des + 1) = (/pXe, pYs, pZb/)
         vertex(3,:,n_facets_des + 1) = (/pXe, pYn, pZb/)

         vertex(1,:,n_facets_des + 2) = (/pXe, pYn, pZb/)
         vertex(2,:,n_facets_des + 2) = (/pXw, pYn, pZb/)
         vertex(3,:,n_facets_des + 2) = (/pXw, pYs, pZb/)

         if(pZb == pZt) then
            if(pZb == 0.0d0) then
               norm_face(:,n_facets_des + 1) = (/zero, zero, one/)
               norm_face(:,n_facets_des + 2) = (/zero, zero, one/)
            elseif(pZb == zLength) then
               norm_face(:,n_facets_des + 1) = (/zero, zero,-one/)
               norm_face(:,n_facets_des + 2) = (/zero, zero,-one/)
            else
               stop 9802
            endif
         else
            norm_face(:,n_facets_des + 1) = (/zero, zero,-one/)
            norm_face(:,n_facets_des + 2) = (/zero, zero,-one/)
         endif
         n_facets_des = n_facets_des+2

! Top Face
         if(pZb /= pZt) then
            vertex(3,:,n_facets_des + 1) = (/pXw, pYs, pZb/)
            vertex(2,:,n_facets_des + 1) = (/pXe, pYs, pZb/)
            vertex(1,:,n_facets_des + 1) = (/pXe, pYn, pZb/)

            vertex(3,:,n_facets_des + 2) = (/pXe, pYn, pZb/)
            vertex(2,:,n_facets_des + 2) = (/pXw, pYn, pZb/)
            vertex(1,:,n_facets_des + 2) = (/pXw, pYs, pZb/)

            norm_face(:,n_facets_des + 1) = (/zero, zero, one/)
            norm_face(:,n_facets_des + 2) = (/zero, zero, one/)

            n_facets_des = n_facets_des + 2
         endif
      endif

      RETURN
      END SUBROUTINE GENERATE_STL_BOX


      END MODULE STL_PREPROC_DES
