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


      use amrex_constants_module, only: M_PI

      use usr, only: init_vel_t, init_angle

      IMPLICIT NONE

      INTEGER :: NP
      double precision :: lTMP(62,3)

      INIT_VEL_T = 0.0d0
      INIT_ANGLE = 0.0d0

! Collect initial translational velocity to IO processor.
      ltmp( 1,:) = (/ -3.900000,   0.000000, 0.000000 /)
      ltmp( 2,:) = (/ -3.899851,   0.034033, 0.000000 /)
      ltmp( 3,:) = (/ -3.899406,   0.066577, 0.014151 /)
      ltmp( 4,:) = (/ -3.898664,   0.093264, 0.041524 /)
      ltmp( 5,:) = (/ -3.897624,   0.110114, 0.080002 /)
      ltmp( 6,:) = (/ -3.896288,   0.113830, 0.126421 /)
      ltmp( 7,:) = (/ -3.885159,   0.169954, 0.294368 /)
      ltmp( 8,:) = (/ -3.866635,   0.157306, 0.484137 /)
      ltmp( 9,:) = (/ -3.840750,   0.070790, 0.673518 /)
      ltmp(10,:) = (/ -3.807554,  -0.088234, 0.839490 /)
      ltmp(11,:) = (/ -3.767111,  -0.311920, 0.959991 /)
      ltmp(12,:) = (/ -3.719496,  -0.586376, 1.015634 /)
      ltmp(13,:) = (/ -3.664801,  -0.892539, 0.991265 /)
      ltmp(14,:) = (/ -3.603130,  -1.207430, 0.877249 /)
      ltmp(15,:) = (/ -3.534600,  -1.505716, 0.670388 /)
      ltmp(16,:) = (/ -3.459342,  -1.761467, 0.374411 /)
      ltmp(17,:) = (/ -3.377499,  -1.950000, 0.000000 /)
      ltmp(18,:) = (/ -3.289227,  -2.049677,-0.435672 /)
      ltmp(19,:) = (/ -3.194693,  -2.043554,-0.909849 /)
      ltmp(20,:) = (/ -3.094078,  -1.920744,-1.395502 /)
      ltmp(21,:) = (/ -2.987573,  -1.677425,-1.862969 /)
      ltmp(22,:) = (/ -2.875382,  -1.317401,-2.281805 /)
      ltmp(23,:) = (/ -2.757716,  -0.852181,-2.622744 /)
      ltmp(24,:) = (/ -2.634802,  -0.300559,-2.859630 /)
      ltmp(25,:) = (/ -2.506872,   0.312286,-2.971207 /)
      ltmp(26,:) = (/ -2.374170,   0.956123,-2.942643 /)
      ltmp(27,:) = (/ -2.236948,   1.597346,-2.766685 /)
      ltmp(28,:) = (/ -2.095468,   2.200922,-2.444372 /)
      ltmp(29,:) = (/ -1.950000,   2.732454,-1.985244 /)
      ltmp(30,:) = (/ -1.800820,   3.160266,-1.407041 /)
      ltmp(31,:) = (/ -1.648211,   3.457361,-0.734885 /)
      ltmp(32,:) = (/ -3.900000,   0.000000, 0.000000 /)
      ltmp(33,:) = (/ -3.899851,   0.033847, 0.003557 /)
      ltmp(34,:) = (/ -3.899406,   0.064733, 0.021033 /)
      ltmp(35,:) = (/ -3.898664,   0.088413, 0.051045 /)
      ltmp(36,:) = (/ -3.897624,   0.101148, 0.091074 /)
      ltmp(37,:) = (/ -3.896288,   0.099991, 0.137626 /)
      ltmp(38,:) = (/ -3.885159,   0.138253, 0.310521 /)
      ltmp(39,:) = (/ -3.866635,   0.105838, 0.497928 /)
      ltmp(40,:) = (/ -3.840750,   0.000000, 0.677228 /)
      ltmp(41,:) = (/ -3.807554,  -0.175501, 0.825669 /)
      ltmp(42,:) = (/ -3.767111,  -0.410558, 0.922128 /)
      ltmp(43,:) = (/ -3.719496,  -0.689327, 0.948777 /)
      ltmp(44,:) = (/ -3.664801,  -0.991265, 0.892539 /)
      ltmp(45,:) = (/ -3.603130,  -1.292513, 0.746233 /)
      ltmp(46,:) = (/ -3.534600,  -1.567542, 0.509325 /)
      ltmp(47,:) = (/ -3.459342,  -1.790955, 0.188237 /)
      ltmp(48,:) = (/ -3.377499,  -1.939318,-0.203831 /)
      ltmp(49,:) = (/ -3.289227,  -1.992909,-0.647535 /)
      ltmp(50,:) = (/ -3.194693,  -1.937254,-1.118474 /)
      ltmp(51,:) = (/ -3.094078,  -1.764352,-1.588630 /)
      ltmp(52,:) = (/ -2.987573,  -1.473502,-2.028102 /)
      ltmp(53,:) = (/ -2.875382,  -1.071670,-2.407011 /)
      ltmp(54,:) = (/ -2.757716,  -0.573361,-2.697454 /)
      ltmp(55,:) = (/ -2.634802,   0.000000,-2.875382 /)
      ltmp(56,:) = (/ -2.506872,   0.621151,-2.922288 /)
      ltmp(57,:) = (/ -2.374170,   1.258475,-2.826581 /)
      ltmp(58,:) = (/ -2.236948,   1.877793,-2.584561 /)
      ltmp(59,:) = (/ -2.095468,   2.444372,-2.200922 /)
      ltmp(60,:) = (/ -1.950000,   2.925000,-1.688750 /)
      ltmp(61,:) = (/ -1.800820,   3.290030,-1.068996 /)
      ltmp(62,:) = (/ -1.648211,   3.515237,-0.369466 /)

! Store the collision angle and initial tangential velocity
      DO NP=1, 62
         INIT_VEL_T(NP) = sqrt(lTMP(NP,2)**2 + lTMP(NP,3)**2)
         INIT_ANGLE(NP) = abs(atan(INIT_VEL_T(NP)/lTMP(NP,1)))*180.0/M_PI
      ENDDO

      END SUBROUTINE USR0
