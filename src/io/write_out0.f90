!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: write_out0                                              !
!  Author: P. Nicoletti, M. Syamlal                   Date: 04-DEC-91  !
!                                                                      !
!  Purpose: Echo user input.                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine write_out0(time, dt, dx, dy, dz, xlength, ylength, zlength, domlo, domhi) &
         bind(C, name="write_out0")

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      use constant, only: gravity
      use discretelement, only: des_continuum_coupled, des_coll_model_enum, hertzian, kn, kt, kn_w, kt_w, lsd
      use discretelement, only: hert_kn, hert_kt, hert_kwn, hert_kwt, des_etan, des_etat, des_etat_wall, des_etan_wall
      use fld_const, only: mw_avg, mu_g0, ro_g0

      use leqsol, only: leq_it, leq_sweep, leq_tol, leq_pc
      use machine, only: id_node, id_month, id_year, id_minute, id_hour, id_day
      use param, only: dim_ic, dim_bc
      use param, only: half, undefined, zero, is_defined
      use constant, only: mmax, ro_s0, d_p0
      use run, only: description, call_usr, dem_solids, dt_fac,  dt_min, dt_max, run_name, run_type, tstop
      use run, only: discretize, solids_model
      use scales, only: p_scale, p_ref
      use ur_facs, only: ur_fac

      use ic, only: write_out_ic
      use bc, only: write_out_bc

      implicit none

      integer(c_int), intent(in   ) :: domlo(3), domhi(3)
      real(c_real)  , intent(in   ) :: time, dt, dx, dy, dz
      real(c_real)  , intent(in   ) :: xlength, ylength, zlength
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: L, M, N

      integer :: MMAX_TOT
      real(c_real) :: TMP_DP

      real(c_real), DIMENSION(6) :: LOC

! Coefficient of restitution (old symbol)
      CHARACTER(LEN=3), DIMENSION(3) :: LEGEND
      CHARACTER(LEN=12), DIMENSION(0:2) :: DISCR_NAME
      CHARACTER(LEN=12), DIMENSION(0:2) :: DISCR_NAME1
! RUN_NAME.OUT file unit number
      integer, PARAMETER :: UNIT_OUT = 52
      integer :: ier
      integer :: i_w, j_s, k_b
      integer :: i_e, j_n, k_t
!-----------------------------------------------

!
      DATA DISCR_NAME/'FOUP', '    ', 'Superbee'/
      DATA DISCR_NAME1/'FOUP', '    ', 'Fourth Order'/

      MMAX_TOT = MMAX

      open(unit=unit_out, file=trim(run_name)//'.out', status='unknown', &
         access='sequential', form='formatted', position='append', iostat=ier)

!  Write Headers for .OUT file
!
      WRITE(UNIT_OUT,1000)
!
!  Echo input data
!
!  Run control section
!
      WRITE (UNIT_OUT, 1100)
      WRITE (UNIT_OUT, 1110) RUN_NAME
      WRITE (UNIT_OUT, 1120) DESCRIPTION
      IF (IS_DEFINED(DT)) THEN
         WRITE (UNIT_OUT, 1135) TIME, TSTOP, DT, DT_MAX, DT_MIN, DT_FAC
      ELSE
         WRITE (UNIT_OUT, 1136)
      ENDIF
      WRITE (UNIT_OUT, 1137) RUN_TYPE
      IF (RUN_TYPE == 'NEW') THEN
         WRITE (UNIT_OUT, 1138)
      ELSE IF (RUN_TYPE == 'RESTART_1') THEN
         WRITE (UNIT_OUT, 1139)
      ENDIF

      WRITE (UNIT_OUT, 1140) 'X', ' '
      WRITE (UNIT_OUT, 1140) 'Y', ' '
      WRITE (UNIT_OUT, 1140) 'Z', ' '

      IF (CALL_USR) THEN
         WRITE (UNIT_OUT, 1149) ' '
      ELSE
         WRITE (UNIT_OUT, 1149) ' NOT '
      ENDIF
!
!  Physical and numerical parameters
!
      WRITE (UNIT_OUT, 1150)
      WRITE (UNIT_OUT, 1157) P_REF, P_SCALE, GRAVITY(2)
      WRITE (UNIT_OUT, 1158)
      WRITE (UNIT_OUT, 1159) (UR_FAC(L),LEQ_IT(L),'BiCGSTAB',&
                          LEQ_SWEEP(L), LEQ_TOL(L), LEQ_PC(L),&
                          DISCR_NAME(DISCRETIZE(L)),L=1,4)

! Geometry and Discretization.
         WRITE (UNIT_OUT, 1200)

         WRITE (UNIT_OUT, 1210)
         LEGEND(1) = '  I'
         LEGEND(2) = ' DX'
         LEGEND(3) = 'X_E'
         CALL WRITE_TABLE (LEGEND, DX, 0.0d0, 1, domhi(1)+1)
         WRITE (UNIT_OUT, 1212) (domhi(1)-domlo(1)+1)
         WRITE (UNIT_OUT, 1213) xlength
         WRITE (UNIT_OUT, 1220)
         LEGEND(1) = '  J'
         LEGEND(2) = ' DY'
         LEGEND(3) = 'Y_N'
         CALL WRITE_TABLE (LEGEND, DY, ZERO, 1, domhi(2)+1)
         WRITE (UNIT_OUT, 1221) (domhi(2)-domlo(2)+1)
         WRITE (UNIT_OUT, 1222) ylength
         WRITE (UNIT_OUT, 1230)
         LEGEND(1) = '  K'
         LEGEND(2) = ' DZ'
         LEGEND(3) = 'Z_T'
         CALL WRITE_TABLE (LEGEND, DZ, ZERO, 1, domhi(3)+1)
         WRITE (UNIT_OUT, 1231) (domhi(3)-domlo(3)+1)
         WRITE (UNIT_OUT, 1232) zlength

!
!  Gas Section
!
      WRITE (UNIT_OUT, 1300)
      IF (IS_DEFINED(RO_G0)) WRITE (UNIT_OUT, 1305) RO_G0
      IF (IS_DEFINED(MU_G0)) WRITE (UNIT_OUT, 1310) MU_G0
      IF (IS_DEFINED(MW_AVG)) WRITE (UNIT_OUT, 1320) MW_AVG
!
!  Particle Section

      WRITE (UNIT_OUT, 1400)
      WRITE (UNIT_OUT, 1401) MMAX_TOT


 1400 FORMAT(//,3X,'5. SOLIDS PHASE',/)
 1401 FORMAT(7X,'Number of particulate phases (MMAX) = ',I2)

      IF(MMAX_TOT > 0) THEN

         WRITE (UNIT_OUT, 1405)
         DO M = 1, MMAX_TOT
            TMP_DP = RO_s0(M)
            WRITE (UNIT_OUT, 1406) M, SOLIDS_MODEL(M), D_P0(M),     &
               TMP_DP, .TRUE.
         END DO


 1405 FORMAT(/7x,'M',4x,'Model',5x,'Diameter',8x,'Density',6x,         &
         'Close_Packed')
 1406 FORMAT(6x,I2,4x,A3,5X,G12.5,3x,G12.5,9x,L1)


         IF(DEM_SOLIDS) THEN
            IF(.NOT.DES_CONTINUUM_COUPLED) THEN
               WRITE(UNIT_OUT,"(/7X,'Gas/Solids NOT coupled.')")
            ELSE
               WRITE(UNIT_OUT,"(/7X,'Gas/Solids Coupling Information:')")
               WRITE(UNIT_OUT,1440) 'cell averaging'
            ENDIF

 1440 FORMAT(10X,'Use ',A,' to calculate gas/particle drag.')

         ENDIF

         IF(DEM_SOLIDS) THEN

 1450 FORMAT(/7X,'Use ',A,' collsion model.',2/10X,&
         'Spring Coefficients:',T37,'Normal',7x,'Tangential')

            IF(DES_COLL_MODEL_ENUM .EQ. LSD) THEN
               WRITE(UNIT_OUT,1450) 'Linear spring-dashpot'
               WRITE(UNIT_OUT,1455) 'Particle-particle', KN, KT
               WRITE(UNIT_OUT,1455) 'Particle-wall', KN_W, KT_W

            ELSEIF(DES_COLL_MODEL_ENUM .EQ. HERTZIAN) THEN
               WRITE(UNIT_OUT,1450) 'Hertzian spring-dashpot'

               DO M = 1, MMAX
                  DO N = M, MMAX
                     IF(M==N) THEN
                       WRITE(UNIT_OUT,1456)M,N,HERT_KN(M,N),HERT_KT(M,N)
                     ELSE
                       WRITE(UNIT_OUT,1457)N,HERT_KN(M,N),HERT_KT(M,N)
                     ENDIF
                  ENDDO
                  WRITE(UNIT_OUT,1458) HERT_KWN(M),HERT_KWT(M)
               ENDDO
            ENDIF

            WRITE(UNIT_OUT,1451)
 1451 FORMAT(/10X,'Damping Coefficients:',T37,'Normal',7x,'Tangential')

            DO M = 1, MMAX
               DO N = M, MMAX
                  IF(M==N) THEN
                     WRITE(UNIT_OUT,1456)M,N,DES_ETAN(M,N),DES_ETAT(M,N)
                  ELSE
                     WRITE(UNIT_OUT,1457)N,DES_ETAN(M,N),DES_ETAT(M,N)
                  ENDIF
               ENDDO
               WRITE(UNIT_OUT,1458) DES_ETAN_WALL(M),DES_ETAT_WALL(M)
            ENDDO

 1455 FORMAT(12X,A,T35,g12.5,3x,g12.5)
 1456 FORMAT(12X,'Phase',I2,'-Phase',I2,' = ',T35,g12.5,3x,g12.5)
 1457 FORMAT(19X,'-Phase',I2,' = ',T35,g12.5,3x,g12.5)
 1458 FORMAT(19X,'-Wall',3x,' = ',T35,g12.5,3x,g12.5)

         ENDIF

      ENDIF


      call write_out_ic(unit_out, dx, dy, dz)

      call write_out_bc(unit_out, dx, dy, dz, &
         xlength, ylength, zlength, domlo, domhi)



      RETURN
 1000 FORMAT(17X,'MM      MM  FFFFFFFFFF    IIIIII    XX      XX',/17X,&
         'MM      MM  FFFFFFFFFF    IIIIII    XX      XX',/17X,&
         'MMMM  MMMM  FF              II      XX      XX',/17X,&
         'MMMM  MMMM  FF              II      XX      XX',/17X,&
         'MM  MM  MM  FF              II        XX  XX  ',/17X,&
         'MM  MM  MM  FF              II        XX  XX  ',/17X,&
         'MM      MM  FFFFFFFF        II          XX    ',/17X,&
         'MM      MM  FFFFFFFF        II          XX    ',/17X,&
         'MM      MM  FF              II        XX  XX  ',/17X,&
         'MM      MM  FF              II        XX  XX  ',/17X,&
         'MM      MM  FF              II      XX      XX',/17X,&
         'MM      MM  FF              II      XX      XX',/17X,&
         'MM      MM  FF            IIIIII    XX      XX',/17X,&
         'MM      MM  FF            IIIIII    XX      XX',2/20X,&
         'Multiphase Flow with Interphase eXchanges')

 1100 FORMAT(//,3X,'1. RUN CONTROL',/)
 1110 FORMAT(7X,'Run name(RUN_NAME): ',A60)
 1120 FORMAT(7X,'Brief description of the run (DESCRIPTION) :',/9X,A60)
 1135 FORMAT(7X,'Start-time (TIME) = ',G12.5,/7X,'Stop_time (TSTOP) = ',G12.5,/7X&
         ,'Time step (DT) = ',G12.5,/7X,'Max time step (DT_MAX) = ',G12.5,/7X&
         ,'Min time step (DT_MIN) = ',G12.5,/7X,&
         'Time step adjustment factor (DT_FAC) = ',G12.5)
 1136 FORMAT(7X,'* Steady state simulation.')
 1137 FORMAT(7X,'Type of run (RUN_TYPE) : ',A16)
 1138 FORMAT(30X,'(Initial conditions from the input (.DAT) file)')
 1139 FORMAT(30X,'(Initial conditions from the restart (.RES) file)')
 1140 FORMAT(/7X,'* Gas momentum equation-',A,' is',A,'solved.')
 1149 FORMAT(/7X,'* User-defined subroutines are',A,'called.')
!
 1150 FORMAT(//,3X,'2. PHYSICAL AND NUMERICAL PARAMETERS',/)
 1157 FORMAT(7X,'Reference pressure (P_ref) = ',G12.5,/7X,&
         'Pressure scale-factor (P_scale) = ',G12.5,/7X,&
         'Gravitational acceleration (GRAVITY) = ',G12.5)
 1158 FORMAT(7X,'Under relaxation (UR_FAC) and',&
         ' Iterations in Leq solver (LEQ_IT):'/,9X,&
         '                        UR_FAC',2X,'LEQ_IT','  LEQ_METHOD',&
         '  LEQ_SWEEP', '  LEQ_TOL', '    LEQ_PC', '  DISCRETIZE')
 1159 FORMAT(9X,&
         'Fluid cont.  and P_g  = ',F6.3,2X,I4,6X,A8,5x,A4,4X,G11.4,3X,A4,3X,A12/9X,&
         'Solids cont. and P_s  = ',F6.3,2X,I4,6X,A8,5x,A4,4X,G11.4,3X,A4,3X,A12/9X,&
         'U velocity            = ',F6.3,2X,I4,6X,A8,5x,A4,4X,G11.4,3X,A4,3X,A12/9X,&
         'V velocity            = ',F6.3,2X,I4,6X,A8,5x,A4,4X,G11.4,3X,A4,3X,A12/9X,&
         'W velocity            = ',F6.3,2X,I4,6X,A8,5x,A4,4X,G11.4,3X,A4,3X,A12/9X,&
         'User scalar           = ',F6.3,2X,I4,6X,A8,5x,A4,4X,G11.4,3X,A4,3X,A12/)
!
 1200 FORMAT(//,3X,'3. GEOMETRY AND DISCRETIZATION',/)

 1210 FORMAT(7X,'X-direction cell sizes (DX) and East face locations:')
 ! 1211 FORMAT(7X,'Minimum value of X, or R (XMIN) =',G12.5)
 1212 FORMAT(7X,'Number of cells in X, or R, direction (IMAX) = ',I4)
 1213 FORMAT(7X,'Reactor length in X, or R, direction (XLENGTH) =',G12.5//)
 1220 FORMAT(7X,'Y-direction cell sizes (DY) and North face locations:')
 1221 FORMAT(7X,'Number of cells in Y direction (JMAX) = ',I4)
 1222 FORMAT(7X,'Reactor length in Y direction (YLENGTH) =',G12.5//)
 1230 FORMAT(7X,'Z-direction cell sizes (DZ) and Top face locations:')
 1231 FORMAT(7X,'Number of cells in Z, or theta, direction (KMAX) = ',I4)
 1232 FORMAT(7X,'Reactor length in Z, or theta, direction (ZLENGTH) =',G12.5)
!
 1300 FORMAT(//,3X,'4. GAS PHASE',/)
 1305 FORMAT(7X,'Gas density (RO_g0) = ',G12.5,&
         '  (A constant value is used everywhere)')
 1310 FORMAT(7X,'Viscosity (MU_g0) = ',G12.5,&
         '  (A constant value is used everywhere)')
 1320 FORMAT(7X,'Average molecular weight (MW_avg) = ',G12.5,&
         '  (A constant value is used everywhere)')

    CONTAINS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_TABLE (LEGEND, ARRAY, DIST_MIN, LSTART, LEND)    C
!  Purpose: To write a table of DX, DY, DZ, and cell wall locations    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_TABLE(LEGEND, SCALAR, DIST_MIN, LSTART, LEND)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

!                      Legend
      CHARACTER(LEN=*)    LEGEND(3)

!                      DX, DY, or DZ Array to be written


!                      Starting array index
      integer          LSTART

!                      Ending array index
      integer          LEND

      real(c_real) SCALAR

!                      Starting value of distance
      real(c_real) DIST_MIN

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------

!                      Number of columns in the table.  When this is changed
!                      remember to change the FORMAT statement also.
      integer, PARAMETER :: NCOL = 5
!
!                      Some dimension large enough for I, J, and K.
      integer, PARAMETER :: DIMENSION_1 = 5000

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      Indices
      integer          ARRAY1(DIMENSION_1)
!
!                      Array3 to be written
      real(c_real) ARRAY3(DIMENSION_1)
!
!                      Number of rows
      integer          NROW
!
!                      Temporary storage for distance calculation
      real(c_real) DIST
!
!                      Local array indices
      integer          L, L1, L2, L3
!-----------------------------------------------
!
!
!  Fill arrays 1 and 3
!
      DIST = DIST_MIN
      DO L = LSTART, LEND
         ARRAY1(L) = L
         ARRAY3(L) = DIST
         IF (L < LEND) DIST = DIST + SCALAR
      END DO
      NROW = (LEND - LSTART + 1)/NCOL
!
      L2 = LSTART - 1
      DO L = 1, NROW
         L1 = L2 + 1
         L2 = L1 + NCOL - 1
         WRITE (UNIT_OUT, 1010) LEGEND(1), (ARRAY1(L3),L3=L1,L2)
         WRITE (UNIT_OUT, 1020) LEGEND(2), (SCALAR,L3=L1,L2)
         WRITE (UNIT_OUT, 1030) LEGEND(3), (ARRAY3(L3),L3=L1,L2)
      END DO
      IF (NROW*NCOL < LEND - LSTART + 1) THEN
         L1 = L2 + 1
         L2 = LEND
         WRITE (UNIT_OUT, 1010) LEGEND(1), (ARRAY1(L3),L3=L1,L2)
         WRITE (UNIT_OUT, 1020) LEGEND(2), (SCALAR,L3=L1,L2)
         WRITE (UNIT_OUT, 1030) LEGEND(3), (ARRAY3(L3),L3=L1,L2)
      ENDIF
      RETURN
!
 1010 FORMAT(7X,A3,2X,5(4X,I3,5X,1X))
 1020 FORMAT(7X,A3,2X,5(G12.5,1X))
 1030 FORMAT(7X,A3,2X,5(G12.5,1X),/)
      END SUBROUTINE WRITE_TABLE
      end subroutine write_out0
