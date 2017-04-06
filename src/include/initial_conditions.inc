!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: INITIAL_CONDITIONS_LOCKED                                 !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that CANNOT be edited during a run.     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / INITIAL_CONDITIONS_LOCKED / IC_DES_FIT_TO_REGION



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: INITIAL_CONDITIONS_UNLOCKED                               !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that may be edited during a run by      !
!  invoking a REFRESH.                                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / INITIAL_CONDITIONS_UNLOCKED /&

      IC_TYPE, &
      IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b, IC_Z_t,&
      IC_EP_g, IC_U_g, IC_V_g, IC_W_g, IC_P_g, &
      IC_EP_s, IC_U_s, IC_V_s, IC_W_s

! Currently unused:
! IC_T_g, IC_T_s, IC_X_g, IC_X_s,
