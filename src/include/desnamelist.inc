       NAMELIST / DES_INPUT_DATA /&

! USE discretelement
        DES_ONEWAY_COUPLED, DES_EXPLICITLY_COUPLED, &

! Particle properties
        KN, KN_W, KT_FAC, KT_W_FAC, MEW, MEW_W, &
        DES_EN_INPUT, DES_ET_INPUT, DES_ETAT_FAC, &
        DES_EN_WALL_INPUT, DES_ET_WALL_INPUT, DES_ETAT_W_FAC, &
!
! Additional properties for Hertzian model:
        DES_COLL_MODEL, EW_YOUNG, VW_POISSON, E_YOUNG, V_POISSON

! DES output control information
