!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: RUN_CONTROL_LOCKED                                        !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that CANNOT be edited during a run.     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / RUN_CONTROL_LOCKED /&

      RUN_NAME, UNITS, RUN_TYPE, TIME


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: RUN_CONTROL_UNLOCKED                                      !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that may be edited during a run by      !
!  invoking a REFRESH.                                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / RUN_CONTROL_UNLOCKED /&

      DESCRIPTION, TSTOP, DT, DT_MAX, DT_MIN, DT_FAC,  &
      MOMENTUM_X_EQ, MOMENTUM_Y_EQ, MOMENTUM_Z_EQ,     &
      DRAG_TYPE, DRAG_C1, DRAG_D1
