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

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

      use constant, only: gravity
      use discretelement, only: des_continuum_coupled, des_coll_model_enum, hertzian, kn, kt, kn_w, kt_w, lsd
      use discretelement, only: hert_kn, hert_kt, hert_kwn, hert_kwt, des_etan, des_etat, des_etat_wall, des_etan_wall
      use fld_const, only: mw_avg, mu_g0, ro_g0

      use leqsol, only: leq_it, leq_sweep, leq_tol, leq_pc
      use param, only: dim_ic, dim_bc
      use param, only: half, undefined, zero, is_defined
      use constant, only: mmax
      use run, only: description, call_usr, dem_solids, dt_fac,  dt_min, dt_max, run_name, tstop
      use run, only: discretize
      use scales, only: p_scale, p_ref
      use ur_facs, only: ur_fac

      use ic, only: write_out_ic
      use bc, only: write_out_bc

      implicit none

      integer(c_int), intent(in   ) :: domlo(3), domhi(3)
      real(rt)  , intent(in   ) :: time, dt, dx, dy, dz
      real(rt)  , intent(in   ) :: xlength, ylength, zlength

      integer :: L, M, N
      integer :: mmax_tot

      character(LEN= 3), DIMENSION(  3) :: legend
      character(LEN=12), DIMENSION(0:2) :: DISCR_NAME

      integer, PARAMETER :: unit_out = 52
      integer :: ier
!-----------------------------------------------
!
      mmax_tot = mmax

      open(unit=unit_out, file=trim(run_name)//'.out', status='unknown', &
         access='sequential', form='formatted', position='append', iostat=ier)

!  Write Headers for .OUT file
!
      write(unit_out,1000)
!
!  Echo input data
!
!  Run control section
!
      write (unit_out, 1100)
      write (unit_out, 1110) RUN_NAME
      write (unit_out, 1120) DESCRIPTION
      IF (IS_DEFINED(DT)) THEN
         write (unit_out, 1135) time, tstop, dt, dt_max, dt_min, dt_fac
      ELSE
         write (unit_out, 1136)
      ENDIF

      write (unit_out, 1140) 'X', ' '
      write (unit_out, 1140) 'Y', ' '
      write (unit_out, 1140) 'Z', ' '

      IF (CALL_USR) THEN
         write (unit_out, 1149) ' '
      ELSE
         write (unit_out, 1149) ' NOT '
      ENDIF
!
!  Physical and numerical parameters
!
      write (unit_out, 1150)
      write (unit_out, 1157) P_REF, P_SCALE, GRAVITY(2)
      write (unit_out, 1158)
      write (unit_out, 1159) (UR_FAC(L),leq_it(L),'BiCGSTAB',&
                          leq_sweep(L), leq_tol(L), leq_pc(L),&
                          DISCR_NAME(DISCRETIZE(L)),L=1,4)

! Geometry and Discretization.
         write (unit_out, 1200)

         write (unit_out, 1210)
         legend(1) = '  I'
         legend(2) = ' DX'
         legend(3) = 'X_E'
         CALL write_table (legend, DX, 0.0d0, 1, domhi(1)+1)
         write (unit_out, 1212) (domhi(1)-domlo(1)+1)
         write (unit_out, 1213) xlength
         write (unit_out, 1220)
         legend(1) = '  J'
         legend(2) = ' DY'
         legend(3) = 'Y_N'
         CALL write_table (legend, DY, ZERO, 1, domhi(2)+1)
         write (unit_out, 1221) (domhi(2)-domlo(2)+1)
         write (unit_out, 1222) ylength
         write (unit_out, 1230)
         legend(1) = '  K'
         legend(2) = ' DZ'
         legend(3) = 'Z_T'
         CALL write_table (legend, DZ, ZERO, 1, domhi(3)+1)
         write (unit_out, 1231) (domhi(3)-domlo(3)+1)
         write (unit_out, 1232) zlength

!
!  Gas Section
!
      write (unit_out, 1300)
      IF (IS_DEFINED(RO_G0)) write (unit_out, 1305) RO_G0
      IF (IS_DEFINED(MU_G0)) write (unit_out, 1310) MU_G0
      IF (IS_DEFINED(MW_AVG)) write (unit_out, 1320) MW_AVG
!
!  Particle Section

      write (unit_out, 1400)
      write (unit_out, 1401) MMAX_TOT


 1400 FORMAT(//,3X,'5. SOLIDS PHASE',/)
 1401 FORMAT(7X,'Number of particulate phases (MMAX) = ',I2)

      IF(MMAX_TOT > 0) THEN

         IF(DEM_SOLIDS) THEN
            IF(.NOT.DES_CONTINUUM_COUPLED) THEN
               write(unit_out,"(/7X,'Gas/Solids NOT coupled.')")
            ELSE
               write(unit_out,"(/7X,'Gas/Solids Coupling Information:')")
               write(unit_out,1440) 'cell averaging'
            ENDIF

 1440 FORMAT(10X,'Use ',A,' to calculate gas/particle drag.')

         ENDIF

         IF(DEM_SOLIDS) THEN

 1450 FORMAT(/7X,'Use ',A,' collsion model.',2/10X,&
         'Spring Coefficients:',T37,'Normal',7x,'Tangential')

            IF(DES_COLL_MODEL_ENUM .EQ. LSD) THEN
               write(unit_out,1450) 'Linear spring-dashpot'
               write(unit_out,1455) 'Particle-particle', KN, KT
               write(unit_out,1455) 'Particle-wall', KN_W, KT_W

            ELSEIF(DES_COLL_MODEL_ENUM .EQ. HERTZIAN) THEN
               write(unit_out,1450) 'Hertzian spring-dashpot'

               do M = 1, MMAX
                  do N = M, MMAX
                     IF(M==N) THEN
                       write(unit_out,1456)M,N,HERT_KN(M,N),HERT_KT(M,N)
                     ELSE
                       write(unit_out,1457)N,HERT_KN(M,N),HERT_KT(M,N)
                     ENDIF
                  ENDdo
                  write(unit_out,1458) HERT_KWN(M),HERT_KWT(M)
               ENDdo
            ENDIF

            write(unit_out,1451)
 1451 FORMAT(/10X,'Damping Coefficients:',T37,'Normal',7x,'Tangential')

            do M = 1, MMAX
               do N = M, MMAX
                  IF(M==N) THEN
                     write(unit_out,1456)M,N,DES_ETAN(M,N),DES_ETAT(M,N)
                  ELSE
                     write(unit_out,1457)N,DES_ETAN(M,N),DES_ETAT(M,N)
                  ENDIF
               ENDdo
               write(unit_out,1458) DES_ETAN_WALL(M),DES_ETAT_WALL(M)
            ENDdo

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
 1140 FORMAT(/7X,'* Gas momentum equation-',A,' is',A,'solved.')
 1149 FORMAT(/7X,'* User-defined subroutines are',A,'called.')
!
 1150 FORMAT(//,3X,'2. PHYSICAL AND NUMERICAL PARAMETERS',/)
 1157 FORMAT(7X,'Reference pressure (P_ref) = ',G12.5,/7X,&
         'Pressure scale-factor (P_scale) = ',G12.5,/7X,&
         'Gravitational acceleration (GRAVITY) = ',G12.5)
 1158 FORMAT(7X,'Under relaxation (UR_FAC) and',&
         ' Iterations in Leq solver (leq_it):'/,9X,&
         '                        UR_FAC',2X,'leq_it','  leq_method',&
         '  leq_sweep', '  leq_tol', '    leq_pc', '  DISCRETIZE')
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

    contains


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: write_table (legend, ARRAY, DIST_MIN, LSTART, LEND)    C
!  Purpose: To write a table of DX, DY, DZ, and cell wall locations    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine write_table(legend, SCALAR, DIST_MIN, LSTART, LEND)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

!                      Legend
      CHARACTER(LEN=*)    legend(3)

!                      DX, DY, or DZ Array to be written


!                      Starting array index
      integer          LSTART

!                      Ending array index
      integer          LEND

      real(rt) SCALAR

!                      Starting value of distance
      real(rt) DIST_MIN

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
      real(rt) ARRAY3(DIMENSION_1)
!
!                      Number of rows
      integer          NROW
!
!                      Temporary storage for distance calculation
      real(rt) DIST
!
!                      Local array indices
      integer          L, L1, L2, L3
!-----------------------------------------------
!
!
!  Fill arrays 1 and 3
!
      DIST = DIST_MIN
      do L = LSTART, LEND
         ARRAY1(L) = L
         ARRAY3(L) = DIST
         IF (L < LEND) DIST = DIST + SCALAR
      end do
      NROW = (LEND - LSTART + 1)/NCOL
!
      L2 = LSTART - 1
      do L = 1, NROW
         L1 = L2 + 1
         L2 = L1 + NCOL - 1
         write (unit_out, 1010) legend(1), (ARRAY1(L3),L3=L1,L2)
         write (unit_out, 1020) legend(2), (SCALAR,L3=L1,L2)
         write (unit_out, 1030) legend(3), (ARRAY3(L3),L3=L1,L2)
      end do
      if (NROW*NCOL < LEND - LSTART + 1) THEN
         L1 = L2 + 1
         L2 = LEND
         write (unit_out, 1010) legend(1), (ARRAY1(L3),L3=L1,L2)
         write (unit_out, 1020) legend(2), (SCALAR,L3=L1,L2)
         write (unit_out, 1030) legend(3), (ARRAY3(L3),L3=L1,L2)
      end if

      return
!
 1010 FORMAT(7X,A3,2X,5(4X,I3,5X,1X))
 1020 FORMAT(7X,A3,2X,5(G12.5,1X))
 1030 FORMAT(7X,A3,2X,5(G12.5,1X),/)
      end subroutine write_table

      end subroutine write_out0
