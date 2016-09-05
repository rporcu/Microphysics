       NAMELIST / DES_INPUT_DATA /&

! USE discretelement
        DES_ONEWAY_COUPLED, DES_REPORT_MASS_INTERP, DES_INTG_METHOD, &
        PARTICLES, NFACTOR, &
        GENER_PART_CONFIG, DES_NEIGHBOR_SEARCH, NEIGHBOR_SEARCH_N,&
        FACTOR_RLM, NEIGHBOR_SEARCH_RAD_RATIO, DES_USR_VAR_SIZE,&

        DES_INTERP_SCHEME, DES_INTERP_WIDTH, DES_INTERP_ON, &
        DES_INTERP_MEAN_FIELDS, DES_EXPLICITLY_COUPLED, &

! Used for grid/cell based search methods
        DESGRIDSEARCH_IMAX, DESGRIDSEARCH_JMAX, DESGRIDSEARCH_KMAX, &

! Particle properties
        KN, KN_W, KT_FAC, KT_W_FAC, MEW, MEW_W, &
        DES_EN_INPUT, DES_ET_INPUT, DES_ETAT_FAC, &
        DES_EN_WALL_INPUT, DES_ET_WALL_INPUT, DES_ETAT_W_FAC, &
!
! Additional properties for Hertzian model:
        DES_COLL_MODEL, EW_YOUNG, VW_POISSON, E_YOUNG, V_POISSON,&

! DES flow boundary condition information: J.Musser
        FORCE_ORD_BC, &

! DES output control information
        PRINT_DES_DATA, DES_OUTPUT_TYPE,  &


! Flag to turn on/off optimizing the list of facets at each des grid cell
         MINIMIZE_DES_FACET_LIST
