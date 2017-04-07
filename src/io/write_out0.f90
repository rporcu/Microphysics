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

      use bc, only: bc_hw_g, bc_uw_g, bc_ww_g, bc_hw_g, bc_vw_s, bc_uw_s, bc_vw_g, bc_ww_s, bc_hw_s
      use bc, only: bc_t_g, bc_ep_g, bc_massflow_g, bc_volflow_g, bc_massflow_s, bc_volflow_s
      use bc, only: bc_type, delp_x, delp_y, delp_z, bc_defined, bc_p_g, bc_area, bc_rop_s
      use bc, only: bc_u_g, bc_v_g, bc_w_g
      use bc, only: bc_u_s, bc_v_s, bc_w_s
      use bc, only: bc_x_w, bc_y_n, bc_z_b, bc_x_e, bc_y_s, bc_z_t
      use constant, only: gravity, c_name, c
      use discretelement, only: des_continuum_coupled, des_coll_model_enum, hertzian, kn, kt, kn_w, kt_w, lsd
      use discretelement, only: hert_kn, hert_kt, hert_kwn, hert_kwt, des_etan, des_etat, des_etat_wall, des_etan_wall
      use fld_const, only: mw_avg, mu_g0, ro_g0
      use geometry, only: coordinates
      use geometry, only: cyclic_x, cyclic_y, cyclic_z
      use geometry, only: cyclic_x_pd, cyclic_y_pd, cyclic_z_pd
      use ic, only: ic_ep_g, ic_u_g, ic_v_g, ic_w_g, ic_p_g
      use ic, only: ic_ep_s, ic_u_s, ic_v_s, ic_w_s

      use ic, only: ic_defined
      use ic, only: ic_x_w, ic_x_e, ic_y_n, ic_y_s, ic_z_b, ic_z_t
      use leqsol, only: leq_it, leq_method, leq_sweep, leq_tol, leq_pc
      use machine, only: id_node, id_month, id_year, id_minute, id_hour, id_day
      use output, only: out_dt, res_dt
      use param, only: dimension_c, dimension_ic, dimension_bc
      use param1, only: half, undefined, zero, is_defined
      use constant, only: mmax, ro_s0, d_p0
      use run, only: description, id_version, call_usr, dem_solids, dt_fac,  dt_min, dt_max, run_name, run_type, tstop
      use run, only: units, discretize, solids_model
      use scales, only: p_scale, p_ref
      use toleranc, only: tol_com, zero_ep_s
      use ur_facs, only: ur_fac


      use ic, only: nsw_, fsw_, psw_
      use ic, only: pinf_, pout_
      use ic, only: minf_, mout_

      use calc_cell_module, only: calc_cell
      use calc_cell_module, only: calc_cell_bc_flow
      use calc_cell_module, only: calc_cell_bc_wall

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
      CHARACTER(LEN=8), DIMENSION(1:4) :: LEQ_METHOD_NAME
! RUN_NAME.OUT file unit number
      integer, PARAMETER :: UNIT_OUT = 52
      integer :: ier
      integer :: i_w, j_s, k_b
      integer :: i_e, j_n, k_t
!-----------------------------------------------

!
      DATA DISCR_NAME/'FOUP', '    ', 'Superbee'/
      DATA DISCR_NAME1/'FOUP', '    ', 'Fourth Order'/
      DATA LEQ_METHOD_NAME/'   SOR  ', 'BiCGSTAB', '  GMRES ', '   CG   '/

      MMAX_TOT = MMAX

      open(unit=unit_out, file=trim(run_name)//'.out', status='unknown', &
         access='sequential', form='formatted', position='append', iostat=ier)

!  Write Headers for .OUT file
!
      WRITE(UNIT_OUT,1000)ID_VERSION,ID_HOUR,ID_MINUTE,ID_MONTH,ID_DAY,ID_YEAR
      WRITE (UNIT_OUT, 1010) ID_NODE(1:50)
!
!  Echo input data
!
!  Run control section
!
      WRITE (UNIT_OUT, 1100)
      WRITE (UNIT_OUT, 1110) RUN_NAME
      WRITE (UNIT_OUT, 1120) DESCRIPTION
      WRITE (UNIT_OUT, 1130) UNITS
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
      WRITE (UNIT_OUT, 1159) (UR_FAC(L),LEQ_IT(L),&
                          LEQ_METHOD_NAME(LEQ_METHOD(L)),&
                          LEQ_SWEEP(L), LEQ_TOL(L), LEQ_PC(L),&
                          DISCR_NAME(DISCRETIZE(L)),L=1,9)

      DO L = 1, DIMENSION_C
         IF (IS_DEFINED(C(L))) WRITE (UNIT_OUT, 1190) C_NAME(L), L, C(L)
      END DO

! Geometry and Discretization.
         WRITE (UNIT_OUT, 1200)
         WRITE (UNIT_OUT, 1201) COORDINATES
         IF (CYCLIC_X_PD) THEN
            WRITE (UNIT_OUT, 1202) 'X', ' with pressure drop'
            WRITE (UNIT_OUT, 1203) 'X', DELP_X
         ELSE IF (CYCLIC_X) THEN
            WRITE (UNIT_OUT, 1202) 'X'
         ENDIF
         IF (CYCLIC_Y_PD) THEN
            WRITE (UNIT_OUT, 1202) 'Y', ' with pressure drop'
            WRITE (UNIT_OUT, 1203) 'Y', DELP_Y
         ELSE IF (CYCLIC_Y) THEN
            WRITE (UNIT_OUT, 1202) 'Y'
         ENDIF
         IF (CYCLIC_Z_PD) THEN
            WRITE (UNIT_OUT, 1202) 'Z', ' with pressure drop'
            WRITE (UNIT_OUT, 1203) 'Z', DELP_Z
         ELSE IF (CYCLIC_Z) THEN
            WRITE (UNIT_OUT, 1202) 'Z'
         ENDIF
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

!
!  Initial Conditions Section
!
      WRITE (UNIT_OUT, 1500)
      DO L = 1, DIMENSION_IC
         IF (IC_DEFINED(L)) THEN

            i_w = calc_cell (ic_x_w(l), dx) + 1
            i_e = calc_cell (ic_x_e(l), dx)
            j_s = calc_cell (ic_y_s(l), dy) + 1
            j_n = calc_cell (ic_y_n(l), dy)
            k_b = calc_cell (ic_z_b(l), dz) + 1
            k_t = calc_cell (ic_z_t(l), dz)

            WRITE (UNIT_OUT, 1510) L
            LOC(1) = LOCATION(I_W,DX) - HALF*DX
            LOC(2) = LOCATION(I_E,DX) + HALF*DX
            LOC(3) = LOCATION(J_S,DY) - HALF*DY
            LOC(4) = LOCATION(J_N,DY) + HALF*DY
            LOC(5) = LOCATION(K_B,DZ) - HALF*DZ
            LOC(6) = LOCATION(K_T,DZ) + HALF*DZ
            WRITE (UNIT_OUT, 1520) &
               IC_X_W(L), LOC(1), IC_X_E(L), LOC(2), &
               IC_Y_S(L), LOC(3), IC_Y_N(L), LOC(4), &
               IC_Z_B(L), LOC(5), IC_Z_T(L), LOC(6)
            WRITE (UNIT_OUT, 1530) I_W, I_E, J_S, J_N, K_B, K_T
            WRITE (UNIT_OUT, 1540) IC_EP_G(L)
            IF (IS_DEFINED(IC_P_G(L))) WRITE (UNIT_OUT, 1541) IC_P_G(L)
!
            WRITE (UNIT_OUT, 1550) IC_U_G(L), IC_V_G(L), IC_W_G(L)
            DO M = 1, MMAX_TOT
               WRITE (UNIT_OUT, 1560) M, IC_EP_S(L,M)
            END DO

            DO M = 1, MMAX_TOT
               WRITE(UNIT_OUT,1570)M,IC_U_S(L,M),M,IC_V_S(L,M),M,IC_W_S(L,M)
            END DO
         ENDIF
      END DO

! Boundary Condition Data
      write (unit_out, 1600)
      do l = 1, dimension_bc
         if (bc_defined(l)) then
            write (unit_out, 1610) l
            write (unit_out, 1611) bc_type(l)
            select case (trim(bc_type(l)))
            case ('MASS_INFLOW','MI')
               write (unit_out, 1612)
               call calc_cell_bc_flow(l, xlength, ylength, zlength, &
                  dx, dy, dz, i_w, i_e, j_s, j_n, k_b, k_t)
            case ('MASS_OUTFLOW','MO')
               write (unit_out, 1613)
               call calc_cell_bc_flow(l, xlength, ylength, zlength, &
                  dx, dy, dz, i_w, i_e, j_s, j_n, k_b, k_t)
            case ('P_INFLOW','PI')
               write (unit_out, 1614)
               call calc_cell_bc_flow(l, xlength, ylength, zlength, &
                  dx, dy, dz, i_w, i_e, j_s, j_n, k_b, k_t)
            case ('P_OUTFLOW','PO')
               write (unit_out, 1615)
               call calc_cell_bc_flow(l, xlength, ylength, zlength, &
                  dx, dy, dz, i_w, i_e, j_s, j_n, k_b, k_t)
               call calc_cell_bc_flow(l, xlength, ylength, zlength, &
                  dx, dy, dz, i_w, i_e, j_s, j_n, k_b, k_t)
            case ('OUTFLOW','OF')
               write (unit_out, 1619)
               call calc_cell_bc_flow(l, xlength, ylength, zlength, &
                  dx, dy, dz, i_w, i_e, j_s, j_n, k_b, k_t)
            case ('FREE_SLIP_WALL','FSW')
               write (unit_out, 1616)
               call calc_cell_bc_wall(l, domlo, domhi, xlength, ylength, &
                  zlength, dx, dy, dz, i_w, i_e, j_s, j_n, k_b, k_t)
            case ('NO_SLIP_WALL','NSW')
               write (unit_out, 1617)
               call calc_cell_bc_wall(l, domlo, domhi, xlength, ylength, &
                  zlength, dx, dy, dz, i_w, i_e, j_s, j_n, k_b, k_t)
            case ('PAR_SLIP_WALL','PSW')
               write (unit_out, 1618)
               call calc_cell_bc_wall(l, domlo, domhi, xlength, ylength, &
                  zlength, dx, dy, dz, i_w, i_e, j_s, j_n, k_b, k_t)
            end select

            write (unit_out, 1620) &
               bc_x_w(l), location(i_w,dx) - half*dx, &
               bc_x_e(l), location(i_e,dx) + half*dx, &
               bc_y_s(l), location(j_s,dy) - half*dy, &
               bc_y_n(l), location(j_n,dy) + half*dy, &
               bc_z_b(l), location(k_b,dz) - half*dz, &
               bc_z_t(l), location(k_t,dz) + half*dz

            write (unit_out, 1630) i_w, i_e, j_s, j_n, k_b, k_t
            write (unit_out,1635)  bc_area(l)

            IF (IS_DEFINED(BC_EP_G(L))) WRITE (UNIT_OUT, 1640) BC_EP_G(L)
            IF (IS_DEFINED(BC_P_G(L)))  WRITE (UNIT_OUT, 1641) BC_P_G(L)
            IF (IS_DEFINED(BC_T_G(L)))  WRITE (UNIT_OUT, 1642) BC_T_G(L)
            IF (IS_DEFINED(BC_MASSFLOW_G(L))) &
               WRITE (UNIT_OUT, 1648) BC_MASSFLOW_G(L)
            IF (IS_DEFINED(BC_VOLFLOW_G(L))) &
               WRITE (UNIT_OUT, 1649) BC_VOLFLOW_G(L)
            IF (IS_DEFINED(BC_U_G(L)))  WRITE (UNIT_OUT, 1650) BC_U_G(L)
            IF (IS_DEFINED(BC_V_G(L)))  WRITE (UNIT_OUT, 1651) BC_V_G(L)
            IF (IS_DEFINED(BC_W_G(L)))  WRITE (UNIT_OUT, 1652) BC_W_G(L)
            DO M = 1, MMAX_TOT
               IF (IS_DEFINED(BC_ROP_S(L,M))) THEN
                  WRITE (UNIT_OUT, "(' ')")
                  WRITE (UNIT_OUT, 1660) M, BC_ROP_S(L,M)
               ENDIF
            END DO
            DO M = 1, MMAX_TOT
               WRITE (UNIT_OUT, "(' ')")
               IF (IS_DEFINED(BC_MASSFLOW_S(L,M))) WRITE (UNIT_OUT, 1668) M, &
                  BC_MASSFLOW_S(L,M)
               IF (IS_DEFINED(BC_VOLFLOW_S(L,M))) WRITE (UNIT_OUT, 1669) M, &
                  BC_VOLFLOW_S(L,M)
               IF(IS_DEFINED(BC_U_S(L,M))) WRITE(UNIT_OUT,1670)M,BC_U_S(L,M)
               IF(IS_DEFINED(BC_V_S(L,M))) WRITE(UNIT_OUT,1671)M,BC_V_S(L,M)
               IF(IS_DEFINED(BC_W_S(L,M))) WRITE(UNIT_OUT,1672)M,BC_W_S(L,M)
            END DO
            IF (BC_TYPE(L) == 'PAR_SLIP_WALL' .or. BC_TYPE(L) == 'PSW') THEN
               WRITE (UNIT_OUT, 1675) BC_HW_G(L), BC_UW_G(L), BC_VW_G(L), &
                  BC_WW_G(L)
               DO M = 1, MMAX_TOT
                  WRITE (UNIT_OUT, 1676) M, BC_HW_S(L,M), BC_UW_S(L,M), BC_VW_S&
                     (L,M), BC_WW_S(L,M)
               END DO
            ENDIF
         ENDIF
      END DO

!
!  Print out file descriptions and write intervals.
!
      WRITE (UNIT_OUT, 1800)
      WRITE (UNIT_OUT, 1801) &
         '.OUT','This file (ASCII)',OUT_DT
      WRITE (UNIT_OUT, 1801) &
         '.LOG','Log file containing messages (ASCII)',UNDEFINED
      WRITE (UNIT_OUT, 1801) &
         '.RES','Restart file (Binary)', RES_DT
!
!  Print out tolerance values from TOLERANCE.INC
!
      WRITE (UNIT_OUT, 1900)
      WRITE (UNIT_OUT, 1901) ZERO_EP_S
      WRITE (UNIT_OUT, 1905) TOL_COM


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
         'Multiphase Flow with Interphase eXchanges'/34X,'Version: ',A,/20X,&
         'Time: ',I2,':',I2,20X,'Date: ',I2,'-',I2,'-',I4)
 1010 FORMAT(/7X,'Computer : ',A50,/,1X,79('_'))
 1100 FORMAT(//,3X,'1. RUN CONTROL',/)
 1110 FORMAT(7X,'Run name(RUN_NAME): ',A60)
 1120 FORMAT(7X,'Brief description of the run (DESCRIPTION) :',/9X,A60)
 1130 FORMAT(7X,'Units (UNITS) : ',A16)
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
 1190 FORMAT(7X,1A20,'- C(',I2,') = ',G12.5)
!
 1200 FORMAT(//,3X,'3. GEOMETRY AND DISCRETIZATION',/)
 1201 FORMAT(7X,'Coordinates: ',1A16/)
 1202 FORMAT(7X,'Cyclic boundary conditions in ',A,' direction',A)
 1203 FORMAT(7X,'Pressure drop (DELP_',A,') = ',G12.5)
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
!
!
 1500 FORMAT(//,3X,'6. INITIAL CONDITIONS')
 1510 FORMAT(/7X,'Initial condition no : ',I4)
 1520 FORMAT(9X,39X,' Specified  ',5X,' Simulated  ',/9X,&
         'X coordinate of west face   (IC_X_w) = ',G12.5,5X,G12.5/,9X,&
         'X coordinate of east face   (IC_X_e) = ',G12.5,5X,G12.5/,9X,&
         'Y coordinate of south face  (IC_Y_s) = ',G12.5,5X,G12.5/,9X,&
         'Y coordinate of north face  (IC_Y_n) = ',G12.5,5X,G12.5/,9X,&
         'Z coordinate of bottom face (IC_Z_b) = ',G12.5,5X,G12.5/,9X,&
         'Z coordinate of top face    (IC_Z_t) = ',G12.5,5X,G12.5)
 1530 FORMAT(9X,'I index of cell at west   (IC_I_w) = ',I4,/,9X,&
         'I index of cell at east   (IC_I_e) = ',I4,/,9X,&
         'J index of cell at south  (IC_J_s) = ',I4,/,9X,&
         'J index of cell at north  (IC_J_n) = ',I4,/,9X,&
         'K index of cell at bottom (IC_K_b) = ',I4,/,9X,&
         'K index of cell at top    (IC_K_t) = ',I4)
 1540 FORMAT(9X,'Void fraction (IC_EP_g) = ',G12.5)
 1541 FORMAT(9X,'Gas pressure (IC_P_g) = ',G12.5)
 1550 FORMAT(9X,'X-component of gas velocity (IC_U_g) = ',G12.5,/9X,&
         'Y-component of gas velocity (IC_V_g) = ',G12.5,/9X,&
         'Z-component of gas velocity (IC_W_g) = ',G12.5)
 1560 FORMAT(9X,'Solids phase-',I2,' Volume fr. (IC_EP_s) = ',G12.5)
 1570 FORMAT(9X,'X-component of solids phase-',I2,' velocity (IC_U_s) =',G12.5,&
         /9X,'Y-component of solids phase-',I2,' velocity (IC_V_s) =',G12.5,/9X&
         ,'Z-component of solids phase-',I2,' velocity (IC_W_s) =',G12.5)
!
 1600 FORMAT(//,3X,'7. BOUNDARY CONDITIONS')
 1610 FORMAT(/7X,'Boundary condition no : ',I4)
 1611 FORMAT(9X,'Type of boundary condition : ',A16)
 1612 FORMAT(11X,'(Inlet with specified gas and solids mass flux)')
 1613 FORMAT(11X,'(Outlet with specified gas and solids mass flux)')
 1614 FORMAT(11X,'(Inlet with specified gas pressure)')
 1615 FORMAT(11X,'(Outlet with specified gas pressure)')
 1616 FORMAT(11X,'(Gradients of parallel velocity components are zero)')
 1617 FORMAT(11X,'(Velocity is zero at wall)')
 1618 FORMAT(11X,'(Partial slip condition at wall)')
 1619 FORMAT(11X,'(Outflow condition)')
 1620 FORMAT(9X,39X,' Specified  ',5X,' Simulated  ',/9X,&
         'X coordinate of west face   (BC_X_w) = ',G12.5,5X,G12.5/,9X,&
         'X coordinate of east face   (BC_X_e) = ',G12.5,5X,G12.5/,9X,&
         'Y coordinate of south face  (BC_Y_s) = ',G12.5,5X,G12.5/,9X,&
         'Y coordinate of north face  (BC_Y_n) = ',G12.5,5X,G12.5/,9X,&
         'Z coordinate of bottom face (BC_Z_b) = ',G12.5,5X,G12.5/,9X,&
         'Z coordinate of top face    (BC_Z_t) = ',G12.5,5X,G12.5)
 1630 FORMAT(9X,'I index of cell at west   (BC_I_w) = ',I4,/,9X,&
         'I index of cell at east   (BC_I_e) = ',I4,/,9X,&
         'J index of cell at south  (BC_J_s) = ',I4,/,9X,&
         'J index of cell at north  (BC_J_n) = ',I4,/,9X,&
         'K index of cell at bottom (BC_K_b) = ',I4,/,9X,&
         'K index of cell at top    (BC_K_t) = ',I4)
 1635 FORMAT(9X,'Boundary area = ',G12.5)
 1640 FORMAT(9X,'Void fraction (BC_EP_g) = ',G12.5)
 1641 FORMAT(9X,'Gas pressure (BC_P_g) = ',G12.5)
 1642 FORMAT(9X,'Gas temperature (BC_T_g) = ',G12.5)
 1648 FORMAT(9X,'Gas mass flow rate (BC_MASSFLOW_g) = ',G12.5)
 1649 FORMAT(9X,'Gas volumetric flow rate (BC_VOLFLOW_g) = ',G12.5)
 1650 FORMAT(9X,'X-component of gas velocity (BC_U_g) = ',G12.5)
 1651 FORMAT(9X,'Y-component of gas velocity (BC_V_g) = ',G12.5)
 1652 FORMAT(9X,'Z-component of gas velocity (BC_W_g) = ',G12.5)
 1660 FORMAT(9X,'Solids phase-',I2,' Density x Volume fr. (BC_ROP_s) = ',G12.5)

 1668 FORMAT(9X,'Solids phase-',I2,' mass flow rate (BC_MASSFLOW_s) =',G12.5)
 1669 FORMAT(9X,'Solids phase-',I2,' volumetric flow rate (BC_VOLFLOW_s) =',&
         G12.5)
 1670 FORMAT(9X,'X-component of solids phase-',I2,' velocity (BC_U_s) =',G12.5)
 1671 FORMAT(9X,'Y-component of solids phase-',I2,' velocity (BC_V_s) =',G12.5)
 1672 FORMAT(9X,'Z-component of solids phase-',I2,' velocity (BC_W_s) =',G12.5)
 1675 FORMAT(9X,'Partial slip coefficient   (BC_hw_g) = ',G12.5,/,9X,&
         'Slip velociity U at wall   (BC_Uw_g) = ',G12.5,/,9X,&
         'Slip velociity V at wall   (BC_Vw_g) = ',G12.5,/,9X,&
         'Slip velociity W at wall   (BC_Ww_g) = ',G12.5)
 1676 FORMAT(9X,'Solids phase: ',I2,/,11X,&
         'Partial slip coefficient   (BC_hw_s) = ',G12.5,/,11X,&
         'Slip velociity U at wall   (BC_Uw_s) = ',G12.5,/,11X,&
         'Slip velociity V at wall   (BC_Vw_s) = ',G12.5,/,11X,&
         'Slip velociity W at wall   (BC_Ww_s) = ',G12.5)
!
 1800 FORMAT(//,3X,'9. OUTPUT DATA FILES:',/7X,'Extension',T18,&
         'Description',T59,'Interval for writing')
 1801 FORMAT(7X,A4,T18,A,T61,G12.5)
!
 1900 FORMAT(//,3X,'10. TOLERANCES',/7X,&
         'The following values are specified in the file TOLERANCE.INC.')
 1901 FORMAT(/7X,'Minimum value of EP_s tracked (ZERO_EP_s) = ',G12.5)
 1905 FORMAT(7X,'Tolerance for species and energy balances (TOL_COM) = ',G12.5)
!

    CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LOCATION(L2, DX)                                       C
!  Purpose: Find the cell center location in X, Y, or Z direction for  C
!           the given index L2.                                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      DOUBLE PRECISION function LOCATION (L2, DX)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use param1, only: half
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! Index for which the location is required
      integer :: L2
! Cell sizes (DX, DY, or DZ)
      real(c_real) :: DX
!
!-----------------------------------------------
!
      LOCATION = HALF*DX + DX*dble(L2-1)

      end function LOCATION

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
