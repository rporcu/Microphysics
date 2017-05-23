!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: NUMERICAL_PARAM_UNLOCKED                                  !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that may be edited during a run by      !
!  invoking a REFRESH.                                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / NUMERICAL_PARAM_UNLOCKED /&

      MAX_NIT, NORM_G, TOL_RESID, &
      TOL_DIVERGE, LEQ_IT, LEQ_SWEEP, LEQ_TOL, LEQ_PC, UR_FAC
