!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR0                                                   C
!                                                                      C
!  Purpose: This routine is called before the time loop starts and is  C
!           user-definable.  The user may insert code in this routine  C
!           or call appropriate user defined subroutines.  This        C
!           can be used for setting constants and checking errors in   C
!           data.  This routine is not called from an IJK loop, hence  C
!           all indices are undefined.                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE USR0

      use compar, only: myPE, PE_IO
      use constant, only: PI

      use exit_mod, only: mfix_exit
      use usr, only: init_vel_t, init_angle

      IMPLICIT NONE

      INTEGER :: NP
      double precision :: lTMP(62,3)

      INIT_VEL_T = 0.0d0
      INIT_ANGLE = 0.0d0

! Collect initial translational velocity to IO processor.
      ltmp( 1,:) = (/ 0.000000,   0.000000,  -3.900000  /)
      ltmp( 2,:) = (/ 0.034033,   0.000000,  -3.899851  /)
      ltmp( 3,:) = (/ 0.066577,   0.014151,  -3.899406  /)
      ltmp( 4,:) = (/ 0.093264,   0.041524,  -3.898664  /)
      ltmp( 5,:) = (/ 0.110114,   0.080002,  -3.897624  /)
      ltmp( 6,:) = (/ 0.113830,   0.126421,  -3.896288  /)
      ltmp( 7,:) = (/ 0.169954,   0.294368,  -3.885159  /)
      ltmp( 8,:) = (/ 0.157306,   0.484137,  -3.866635  /)
      ltmp( 9,:) = (/ 0.070790,   0.673518,  -3.840750  /)
      ltmp(10,:) = (/-0.088234,   0.839490,  -3.807554  /)
      ltmp(11,:) = (/-0.311920,   0.959991,  -3.767111  /)
      ltmp(12,:) = (/-0.586376,   1.015634,  -3.719496  /)
      ltmp(13,:) = (/-0.892539,   0.991265,  -3.664801  /)
      ltmp(14,:) = (/-1.207430,   0.877249,  -3.603130  /)
      ltmp(15,:) = (/-1.505716,   0.670388,  -3.534600  /)
      ltmp(16,:) = (/-1.761467,   0.374411,  -3.459342  /)
      ltmp(17,:) = (/-1.950000,   0.000000,  -3.377499  /)
      ltmp(18,:) = (/-2.049677,  -0.435672,  -3.289227  /)
      ltmp(19,:) = (/-2.043554,  -0.909849,  -3.194693  /)
      ltmp(20,:) = (/-1.920744,  -1.395502,  -3.094078  /)
      ltmp(21,:) = (/-1.677425,  -1.862969,  -2.987573  /)
      ltmp(22,:) = (/-1.317401,  -2.281805,  -2.875382  /)
      ltmp(23,:) = (/-0.852181,  -2.622744,  -2.757716  /)
      ltmp(24,:) = (/-0.300559,  -2.859630,  -2.634802  /)
      ltmp(25,:) = (/ 0.312286,  -2.971207,  -2.506872  /)
      ltmp(26,:) = (/ 0.956123,  -2.942643,  -2.374170  /)
      ltmp(27,:) = (/ 1.597346,  -2.766685,  -2.236948  /)
      ltmp(28,:) = (/ 2.200922,  -2.444372,  -2.095468  /)
      ltmp(29,:) = (/ 2.732454,  -1.985244,  -1.950000  /)
      ltmp(30,:) = (/ 3.160266,  -1.407041,  -1.800820  /)
      ltmp(31,:) = (/ 3.457361,  -0.734885,  -1.648211  /)
      ltmp(32,:) = (/ 0.000000,   0.000000,  -3.900000  /)
      ltmp(33,:) = (/ 0.033847,   0.003557,  -3.899851  /)
      ltmp(34,:) = (/ 0.064733,   0.021033,  -3.899406  /)
      ltmp(35,:) = (/ 0.088413,   0.051045,  -3.898664  /)
      ltmp(36,:) = (/ 0.101148,   0.091074,  -3.897624  /)
      ltmp(37,:) = (/ 0.099991,   0.137626,  -3.896288  /)
      ltmp(38,:) = (/ 0.138253,   0.310521,  -3.885159  /)
      ltmp(39,:) = (/ 0.105838,   0.497928,  -3.866635  /)
      ltmp(40,:) = (/ 0.000000,   0.677228,  -3.840750  /)
      ltmp(41,:) = (/-0.175501,   0.825669,  -3.807554  /)
      ltmp(42,:) = (/-0.410558,   0.922128,  -3.767111  /)
      ltmp(43,:) = (/-0.689327,   0.948777,  -3.719496  /)
      ltmp(44,:) = (/-0.991265,   0.892539,  -3.664801  /)
      ltmp(45,:) = (/-1.292513,   0.746233,  -3.603130  /)
      ltmp(46,:) = (/-1.567542,   0.509325,  -3.534600  /)
      ltmp(47,:) = (/-1.790955,   0.188237,  -3.459342  /)
      ltmp(48,:) = (/-1.939318,  -0.203831,  -3.377499  /)
      ltmp(49,:) = (/-1.992909,  -0.647535,  -3.289227  /)
      ltmp(50,:) = (/-1.937254,  -1.118474,  -3.194693  /)
      ltmp(51,:) = (/-1.764352,  -1.588630,  -3.094078  /)
      ltmp(52,:) = (/-1.473502,  -2.028102,  -2.987573  /)
      ltmp(53,:) = (/-1.071670,  -2.407011,  -2.875382  /)
      ltmp(54,:) = (/-0.573361,  -2.697454,  -2.757716  /)
      ltmp(55,:) = (/ 0.000000,  -2.875382,  -2.634802  /)
      ltmp(56,:) = (/ 0.621151,  -2.922288,  -2.506872  /)
      ltmp(57,:) = (/ 1.258475,  -2.826581,  -2.374170  /)
      ltmp(58,:) = (/ 1.877793,  -2.584561,  -2.236948  /)
      ltmp(59,:) = (/ 2.444372,  -2.200922,  -2.095468  /)
      ltmp(60,:) = (/ 2.925000,  -1.688750,  -1.950000  /)
      ltmp(61,:) = (/ 3.290030,  -1.068996,  -1.800820  /)
      ltmp(62,:) = (/ 3.515237,  -0.369466,  -1.648211  /)

      ! Store the collision angle and initial tangential velocity
      IF(myPE == PE_IO) THEN
         DO NP=1, 62
            INIT_VEL_T(NP) = sqrt(lTMP(NP,1)**2 + lTMP(NP,2)**2)
            INIT_ANGLE(NP) = abs(atan(INIT_VEL_T(NP)/lTMP(NP,3)))*180.0/PI
         ENDDO
      ENDIF

      END SUBROUTINE USR0
