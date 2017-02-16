module check_des_solids_module


  use bl_fort_module, only : c_real
  use iso_c_binding , only: c_int
  use param1,         only: ZERO, HALF, ONE, UNDEFINED, IS_UNDEFINED, IS_DEFINED
  use run,            only: IFILE_NAME
  use error_manager,  only: finl_err_msg, flush_err_msg, init_err_msg,  &
       & ival, ivar, err_msg

  implicit none
  private

  public check_solids_dem

contains

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  SUBROUTINE: CHECK_DES_SOLIDS                                        !
  !  Author: J.Musser                                   Date: 02-FEB-14  !
  !                                                                      !
  !  Purpose:                                                            !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_solids_dem

    ! Initialize the error manager.
    call init_err_msg("CHECK_SOLIDS_DEM")
    
    ! Particle-particle collision parameters.
    call check_solids_dem_collision

    call finl_err_msg

  end subroutine check_solids_dem


  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: CHECK_SOLIDS_DEM_COLLISION                              !
  !  Author: J.Musser                                   Date: 11-Dec-13  !
  !                                                                      !
  !  Purpose: Check user input data for DES collision calculations.      !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_solids_dem_collision

    use discretelement, only: DES_COLL_MODEL, DES_COLL_MODEL_ENUM, &
         &  LSD, HERTZIAN, MEW, MEW_W

    ! Initialize the error manager.
    call init_err_msg("CHECK_SOLIDS_DEM_COLLISION")

    ! Check coefficient friction
    if(IS_UNDEFINED(MEW)) then
       write(ERR_MSG,1000) 'MEW', trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    elseif (MEW < ZERO .or. MEW_W > ONE) then
       write(ERR_MSG,1001) 'MEW', trim(iVal(MEW)), trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    endif

    if(IS_UNDEFINED(MEW_W)) then
       write(ERR_MSG,1000) 'MEW_W', trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    elseif(MEW_w < ZERO .or. MEW_W > ONE) then
       write(ERR_MSG,1001) 'MEW_W', trim(iVal(MEW_W)), trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    endif
    
    ! Check collision model specific parameters.
    select case (trim(DES_COLL_MODEL))
       ! Linear spring-dashpot model.
    case('LSD')
       DES_COLL_MODEL_ENUM = LSD
       call check_solids_dem_coll_lsd
       ! Hertzian collision model.
    case('HERTZIAN')
       DES_COLL_MODEL_ENUM = HERTZIAN
       call check_solids_dem_coll_hertz
       ! Unknown collision model.
    case DEFAULT
       write(ERR_MSG,2000) trim(DES_COLL_MODEL), trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    end select
    
2000 format('Error 2000: Invalid particle-particle collision model:',&
         A,/'Please correct the ',A,' file.')


    call finl_err_msg


1000 format('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the ',A,' file.')

1001 format('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
         'Please correct the ',A,' file.')

  end subroutine check_solids_dem_collision

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: CHECK_SOLIDS_DEM_COLL_LSD                               !
  !  Author: J.Musser                                   Date: 11-Dec-13  !
  !                                                                      !
  !  Purpose: Check user input data for DES collision calculations.      !
  !                                                                      !
  !  References:                                                         !
  !   - Schafer et al., J. Phys. I France, 1996, 6, 5-20 (see page 7&13) !
  !   -  Van der Hoef et al., Advances in Chemical Engineering, 2006, 31,!
  !      65-149 (pages 94-95)                                            !
  !   - Silbert et al., Physical Review E, 2001, 64, 051302 1-14 (page 5)!
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_solids_dem_coll_lsd

    use constant,       only: D_P0, RO_s0, PI, MMAX
    use discretelement, only: KN, KN_W, KT, KT_W, KT_FAC, KT_W_FAC, &
         & DES_ETAN, DES_ETAN_WALL, DES_ETAT, DES_ETAT_WALL,        &
         & DES_EN_INPUT, DES_EN_WALL_INPUT, DES_ET_INPUT, DTSOLID,  &
         & DES_ET_WALL_INPUT, DES_ETAT_FAC, DES_ETAT_W_FAC

    integer      :: M, L, LC
    logical      :: FLAG_WARN
    real(c_real) :: TCOLL, TCOLL_TMP     ! Collision length scale.
    real(c_real) :: MASS_M, MASS_L, MASS_EFF
    real(c_real) :: EN     ! Alias for coefficient restitution

    ! Initialize the error manager.
    call init_err_msg("CHECK_SOLIDS_DEM_COLL_LSD")

    ! Initialize.
    TCOLL = UNDEFINED

    ! Check for particle-particle normal spring constants.
    if(IS_UNDEFINED(KN)) then
       write(ERR_MSG, 1000) 'KN', trim(IFILE_NAME)
       call flush_err_msg(ABORT=.TRUE.)
    endif

    ! Check for particle-particle tangential spring constant factors.
    if(IS_UNDEFINED(KT_FAC)) then
       write (ERR_MSG, 2100) 'KT_FAC', trim(IFILE_NAME)
       call flush_err_msg()
       KT_FAC = 2.0d0/7.0d0
    elseif(KT_FAC > ONE .or. KT_FAC < ZERO) then
       write(ERR_MSG,1001) 'KT_FAC', trim(iVal(KT_FAC)), trim(IFILE_NAME)
       call flush_err_msg(ABORT=.TRUE.)
    endif

    ! Calculate the particle-particle tangential spring factor.
    KT = KT_FAC*KN

    ! Check for particle-wall normal spring constants.
    if(IS_UNDEFINED(KN_W)) then
       write(ERR_MSG, 1000) 'KN_W', trim(IFILE_NAME)
       call flush_err_msg(ABORT=.TRUE.)
    endif

    ! Check for particle-wall tangential spring constant factors.
    if(IS_UNDEFINED(KT_W_FAC)) then
       write (ERR_MSG, 2100) 'KT_W_FAC', trim(IFILE_NAME)
       call flush_err_msg()
       KT_W_FAC = 2.0d0/7.0d0
    elseif(KT_W_FAC > ONE .or. KT_W_FAC < ZERO) then
       write(ERR_MSG,1001) 'KT_W_FAC', trim(iVal(KT_W_FAC)), trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    endif

    ! Calculate the particle-wall tangential spring factor.
    KT_W = KT_W_FAC*KN_W

2100 format('Warning 2100: Tangential spring factor ',A,' not ',      &
         'specified in ',A,'.',/'Setting to default: (2/7).')

    ! Check for particle-particle tangential damping coefficients
    if(IS_UNDEFINED(DES_ETAT_FAC)) then
       write (ERR_MSG, 2101) 'DES_ETAT_FAC', trim(IFILE_NAME)
       call flush_err_msg
       DES_ETAT_FAC = HALF
    elseif(DES_ETAT_FAC > ONE .or. DES_ETAT_FAC < ZERO) then
       write(ERR_MSG,1001) 'DES_ETAT_FAC', iVal(DES_ETAT_FAC), trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    endif

    ! Check for particle-wall tangential damping coefficients
    if(IS_UNDEFINED(DES_ETAT_W_FAC)) then
       write (ERR_MSG, 2101) 'DES_ETAT_W_FAC', trim(IFILE_NAME)
       call flush_err_msg
       DES_ETAT_W_FAC = HALF
    elseif(DES_ETAT_W_FAC > ONE .or. DES_ETAT_W_FAC < ZERO) then
       write(ERR_MSG,1001) 'DES_ETAT_W_FAC', iVal(DES_ETAT_W_FAC), trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    endif

2101 format('Warning 2101: Tangential damping factor ',A,' not ',     &
         'specified',/'in ',A,'. Setting to default: (1/2).')


    LC = 0
    do M = 1, MMAX

       ! Calculate the mass of a phase M particle.
       MASS_M = (PI/6.d0)*(D_P0(M)**3)*RO_S0(M)

       ! Particle-Particle Collision Parameters ------------------------------>
       do L = M, MMAX
          LC = LC+1

          ! Check particle-particle normal restitution coefficient
          if(IS_UNDEFINED(DES_EN_INPUT(LC))) then
             write(ERR_MSG,1000) trim(iVar('DES_EN_INPUT',LC)), trim(IFILE_NAME)
             call flush_err_msg(ABORT=.true.)
          elseif(DES_EN_INPUT(LC) > ONE .or.                         &
               DES_EN_INPUT(LC) < ZERO) then
             write(ERR_MSG,1001) trim(iVar('DES_EN_INPUT',LC)),      &
                  trim(iVal(DES_EN_INPUT(LC))), trim(IFILE_NAME)
             call flush_err_msg(ABORT=.true.)
          endif
          EN = DES_EN_INPUT(LC)

          ! Calculate masses used for collision calculations.
          MASS_L = (PI/6.d0)*(D_P0(L)**3)*RO_S0(L)
          MASS_EFF = MASS_M*MASS_L/(MASS_M+MASS_L)

          ! Calculate the M-L normal and tangential damping coefficients.
          if(abs(EN) > ZERO) then
             DES_ETAN(M,L) = 2.0D0*sqrt(KN*MASS_EFF) * abs(log(EN))
             DES_ETAN(M,L) = DES_ETAN(M,L)/sqrt(PI*PI + (log(EN)**2))
          else
             DES_ETAN(M,L) = 2.0D0*sqrt(KN*MASS_EFF)
          endif
          DES_ETAT(M,L) = DES_ETAT_FAC*DES_ETAN(M,L)

          ! Store the entries in the symmetric matrix.
          DES_ETAN(L,M) = DES_ETAN(M,L)
          DES_ETAT(L,M) = DES_ETAT(M,L)

          ! Calculate the collision time scale.
          TCOLL_TMP = PI/sqrt(KN/MASS_EFF -                          &
               ((DES_ETAN(M,L)/MASS_EFF)**2)/4.d0)
          TCOLL = min(TCOLL_TMP, TCOLL)
       end do


       ! Particle-Wall Collision Parameters ---------------------------------->
       ! Check particle-wall normal restitution coefficient.
       if(IS_UNDEFINED(DES_EN_WALL_INPUT(M))) then
          write(ERR_MSG,1000) trim(iVar('DES_EN_WALL_INPUT',M)), trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       elseif(DES_EN_WALL_INPUT(M) > ONE .or.                        &
            DES_EN_WALL_INPUT(M) < ZERO) then
          write(ERR_MSG,1001) trim(iVar('DES_EN_WALL_INPUT',M)),     &
               trim(iVal(DES_EN_WALL_INPUT(M))), trim(IFILE_NAME)
          CALL flush_err_msg(ABORT=.TRUE.)
       endif
       EN = DES_EN_WALL_INPUT(M)

       ! Calculate masses used for collision calculations.
       MASS_EFF = MASS_M

       ! Calculate the M-Wall normal and tangential damping coefficients.
       if(abs(EN) > ZERO) then
          DES_ETAN_WALL(M) = 2.d0*sqrt(KN_W*MASS_EFF)*abs(log(EN))
          DES_ETAN_WALL(M) = DES_ETAN_WALL(M)/sqrt(PI*PI+(log(EN))**2)
       else
          DES_ETAN_WALL(M) = 2.D0*sqrt(KN_W*MASS_EFF)
       endif
       DES_ETAT_WALL(M) = DES_ETAT_W_FAC*DES_ETAN_WALL(M)

       ! Calculate the collision time scale.
       TCOLL_TMP = PI/sqrt(KN_W/MASS_EFF -                           &
            ((DES_ETAN_WALL(M)/MASS_EFF)**2.d0)/4.d0)
       !         TCOLL = MIN(TCOLL_TMP, TCOLL)
    ENDDO

    ! if following are assigned warn user they are discarded
    FLAG_WARN = .false.
    do M = 1, MMAX+MMAX*(MMAX-1)/2
       if(IS_DEFINED(DES_ET_INPUT(M))) FLAG_WARN = .true.
    enddo
    if (FLAG_WARN) then
       write(ERR_MSG,2102) 'DES_ET_INPUT'
       call flush_err_msg
    endif
    
    FLAG_WARN = .false.
    do M = 1, MMAX
       if(IS_DEFINED(DES_ET_WALL_INPUT(M))) FLAG_WARN = .true.
    enddo
    if (FLAG_WARN)then
       write(ERR_MSG,2102) 'DES_ET_WALL_INPUT'
       call flush_err_msg
    endif

2102 FORMAT('Warning 2102: ',A,' values are not used ',/' with the',  &
         ' linear spring-dashpot collision model.')

    ! Store the smalled calculated collision time scale. This value is used
    ! in time-marching the DEM solids.
    DTSOLID = TCOLL/50.d0

    call finl_err_msg

1000 format('Error 1000: Required input not specified: ',A,/'Please ',&
          'correct the ',A,' file.')

1001 format('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
         'Please correct the ',A,' file.')

  end subroutine check_solids_dem_coll_lsd

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: CHECK_SOLIDS_DEM_COLL_HERTZ                             !
  !  Author: J.Musser                                   Date: 11-Dec-13  !
  !                                                                      !
  !  Purpose: Check user input data for Hertzian collisions.             !
  !                                                                      !
  !  References:                                                         !
  !   - Schafer et al., J. Phys. I France, 1996, 6, 5-20 (see page 7&13) !
  !   -  Van der Hoef et al., Advances in Chemical Engineering, 2006, 31,!
  !      65-149 (pages 94-95)                                            !
  !   - Silbert et al., Physical Review E, 2001, 64, 051302 1-14 (page 5)!
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_solids_dem_coll_hertz
    
    use constant,       only: MMAX, D_P0, RO_s0, PI
    USE param,          only: DIM_M
    use discretelement, only: DES_EN_INPUT, DES_EN_WALL_INPUT, &
                           &  DES_ET_INPUT, DES_ET_WALL_INPUT, &
                           &  DES_ETAN, DES_ETAN_WALL,         &
                           &  DES_ETAT, DES_ETAT_WALL,         &
                           &  HERT_KN, HERT_Kwn,               &
                           &  HERT_KT, HERT_Kwt,               &
                           &  DES_ETAT_FAC, DES_ETAT_W_FAC,    &
                           &  E_YOUNG, Ew_YOUNG,               &
                           &  V_POISSON, Vw_POISSON,           &
                           &  DTSOLID,  KN, KN_W, KT_FAC,      &
                           &  KT_W_FAC

    integer           :: M, L, LC
    character(len=64) :: MSG
    real(c_real)      :: TCOLL, TCOLL_TMP         ! Collision length scale.
    real(c_real)      :: MASS_M, MASS_L, MASS_EFF ! Particle and effective mass.
    ! Effective physical quantities. Radius, Youngs, Shear
    real(c_real)      :: R_EFF, E_EFF, G_MOD_EFF, RED_MASS_EFF   
    ! Alias for coefficient restitution
    real(c_real)      :: EN, ET                       
    real(c_real)      :: G_MOD(DIM_M), G_MOD_WALL     ! Shear modules for particles and wall

    ! Initialize the error manager.
    call init_err_msg("CHECK_SOLIDS_DEM_COLL_HERTZ")

    ! Initialize.
    TCOLL = UNDEFINED

    ! check young's modulus and poisson ratio
    if(IS_UNDEFINED(Ew_YOUNG)) then
       MSG='Wall value for Youngs modulus'
       write(ERR_MSG,1002) 'Ew_YOUNG', MSG, trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    endif

    if(IS_UNDEFINED(Vw_POISSON)) then
       MSG='Wall value for Poissons ratio'
       write(ERR_MSG,1002) 'Vw_POISSON', MSG, trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    elseif (Vw_POISSON > 0.5d0 .or. Vw_POISSON <= -ONE) then
       write(ERR_MSG,1001) 'Vw_POISSON',iVal(Vw_POISSON), trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    endif

    G_MOD_WALL = 0.5d0*Ew_YOUNG/(1.d0+Vw_POISSON)

    do M=1,MMAX

       if(IS_UNDEFINED(E_YOUNG(M))) then
          MSG=''; write(MSG,"('Phase ',I2,' Youngs modulus')") M
          WRITE(ERR_MSG,1002) 'E_YOUNG', MSG, trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       endif
       if(IS_UNDEFINED(V_POISSON(M))) then
          MSG=''; write(MSG,"('Phase ',I2,' Poissons ratio')") M
          write(ERR_MSG,1002) 'V_POISSON', MSG, trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       elseif(V_POISSON(M) > 0.5d0 .or.                              &
            V_POISSON(M) <= -ONE) then
          write(ERR_MSG,1001) trim(iVar('V_POISSON',M)),             &
               iVal(V_POISSON(M)), trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       ENDIF
       ! Calculate the shear modulus for phase M.
       G_MOD(M) = 0.5d0*E_YOUNG(M)/(1.d0+V_POISSON(M))
    enddo


    LC = 0
    do M=1,MMAX
       ! Calculate the mass of a phase M particle.
       MASS_M = (PI/6.d0)*(D_P0(M)**3)*RO_S0(M)

       ! Particle-Particle Collision Parameters ------------------------------>
       do L=M,MMAX
          LC = LC+1

          ! Check particle-particle normal restitution coefficient
          if(IS_UNDEFINED(DES_EN_INPUT(LC))) then
             WRITE(ERR_MSG,1000) trim(iVar('DES_EN_INPUT',LC)), trim(IFILE_NAME)
             call flush_err_msg(ABORT=.true.)

          elseif(DES_EN_INPUT(LC) > ONE .or.                         &
               DES_EN_INPUT(LC) < ZERO) then
             write(ERR_MSG,1001) trim(iVar('DES_EN_INPUT',LC)),      &
                  trim(iVal(DES_EN_INPUT(LC))), trim(IFILE_NAME)
             call flush_err_msg(ABORT=.true.)
          endif
          EN = DES_EN_INPUT(LC)

          ! Check particle-particle tangential restitution coefficient
          if(IS_UNDEFINED(DES_ET_INPUT(M))) then
             write(ERR_MSG,1000) trim(iVar('DES_ET_INPUT',M)), trim(IFILE_NAME)
             call flush_err_msg(ABORT=.true.)
          elseif(DES_ET_INPUT(M) > ONE .or.                          &
               DES_ET_INPUT(M) < ZERO) then
             write(ERR_MSG,1001) trim(iVar('DES_ET_INPUT',M)),       &
                  iVal(DES_ET_INPUT(M)), trim(IFILE_NAME)
             call flush_err_msg(ABORT=.true.)
          endif
          ET = DES_ET_INPUT(LC)


          ! Calculate masses used for collision calculations.
          MASS_L = (PI/6.d0)*(D_P0(L)**3)*RO_S0(L)
          MASS_EFF = (MASS_M*MASS_L)/(MASS_M+MASS_L)
          RED_MASS_EFF = (2.d0/7.d0)*MASS_EFF
          ! Calculate the effective radius, Youngs modulus, and shear modulus.
          R_EFF = 0.5d0*(D_P0(M)*D_P0(L)/                    &
               (D_P0(M) + D_P0(L)))
          E_EFF = E_YOUNG(M)*E_YOUNG(L) /                            &
               (E_YOUNG(M)*(1.d0 - V_POISSON(L)**2) +                  &
               E_YOUNG(L)*(1.d0 - V_POISSON(M)**2))
          G_MOD_EFF = G_MOD(M)*G_MOD(L)/                             &
               (G_MOD(M)*(2.d0 - V_POISSON(L)) +                       &
               G_MOD(L)*(2.d0 - V_POISSON(M)))

          ! Calculate the spring properties and store in symmetric matrix format.
          HERT_KN(M,L)=(4.d0/3.d0)*sqrt(R_EFF)*E_EFF
          HERT_KT(M,L)= 8.d0*sqrt(R_EFF)*G_MOD_EFF

          HERT_KN(L,M) = HERT_KN(M,L)
          HERT_KT(L,M) = HERT_KT(M,L)

          ! Calculate the normal coefficient.
          if(abs(EN) > ZERO) then
             DES_ETAN(M,L) = 2.d0*sqrt(HERT_KN(M,L)*MASS_EFF)*       &
                  abs(log(EN))
             DES_ETAN(M,L) = DES_ETAN(M,L)/                          &
                  sqrt(PI*PI + (log(EN))**2)
          else
             DES_ETAN(M,L) = 2.d0*sqrt(HERT_KN(M,L)*MASS_EFF)
          endif
          DES_ETAN(L,M) = DES_ETAN(M,L)

          ! Calculate the tangential coefficients.
          if(abs(ET) > ZERO) then
             DES_ETAT(M,L) = 2.d0*sqrt(HERT_KT(M,L)*RED_MASS_EFF)*   &
                  abs(log(ET))
             DES_ETAT(M,L) = DES_ETAT(M,L)/ sqrt(PI*PI + (log(ET))**2)
          else
             DES_ETAT(M,L) = 2.d0*sqrt(HERT_KT(M,L)*RED_MASS_EFF)
          endif
          DES_ETAT(L,M) = DES_ETAT(M,L)

          TCOLL_TMP = PI/sqrt(HERT_KN(M,L)/MASS_EFF -                &
               ((DES_ETAN(M,L)/MASS_EFF)**2)/4.d0)
          TCOLL = min(TCOLL_TMP, TCOLL)
       enddo
       
       ! Particle-Wall Collision Parameters ---------------------------------->
       ! Check particle-wall normal restitution coefficient.
       if(IS_UNDEFINED(DES_EN_WALL_INPUT(M))) then
          write(ERR_MSG,1000) trim(iVar('DES_EN_WALL_INPUT',M)), trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       elseif(DES_EN_WALL_INPUT(M) > ONE .or.                        &
            DES_EN_WALL_INPUT(M) < ZERO) then
          write(ERR_MSG,1001) trim(iVar('DES_EN_WALL_INPUT',M)),     &
               trim(iVal(DES_EN_WALL_INPUT(M))), trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       endif
       EN = DES_EN_WALL_INPUT(M)

       ! Check particle-wall tangential restitution coefficient
       if(IS_UNDEFINED(DES_ET_WALL_INPUT(M))) then
          write(ERR_MSG,1000) trim(iVar('DES_ET_WALL_INPUT',M)), trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       elseif(DES_ET_WALL_INPUT(M) > ONE .or.                        &
            DES_ET_WALL_INPUT(M) < ZERO) then
          write(ERR_MSG,1001) trim(iVar('DES_ET_WALL_INPUT',M)),     &
               trim(iVal(DES_ET_WALL_INPUT(M))), trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       endif
       ET = DES_ET_WALL_INPUT(M)

       ! Calculate masses used for collision calculations.
       MASS_EFF = MASS_M
       RED_MASS_EFF = (2.d0/7.d0)*MASS_EFF
       ! Calculate the effective radius, Youngs modulus, and shear modulus.
       R_EFF = 0.5d0*D_P0(M)
       E_EFF = E_YOUNG(M)*Ew_YOUNG /                                 &
            (E_YOUNG(M)*(1.d0-Vw_POISSON**2) +                         &
            Ew_YOUNG  *(1.d0-V_POISSON(M)**2))
       G_MOD_EFF = G_MOD(M)*G_MOD_WALL /                             &
            (G_MOD(M)*(2.d0 - Vw_POISSON) +                            &
            G_MOD_WALL*(2.d0 - V_POISSON(M)))

       ! Calculate the spring properties.
       HERT_Kwn(M) = (4.d0/3.d0)*sqrt(R_EFF)*E_EFF
       HERT_Kwt(M) = 8.0*sqrt(R_EFF)*G_MOD_EFF
       
       ! Calculate the tangential coefficients.
       if(abs(EN) > ZERO) then
          DES_ETAN_WALL(M) = 2.d0*sqrt(HERT_Kwn(M)*MASS_EFF)*&
               abs(log(EN))
          DES_ETAN_WALL(M) = DES_ETAN_WALL(M)/&
               sqrt(PI*PI + (log(EN))**2)
       else
          DES_ETAN_WALL(M) = 2.d0*sqrt(HERT_Kwn(M)*MASS_EFF)
       endif

       if(abs(ET) > ZERO) then
          DES_ETAT_WALL(M) = 2.d0*sqrt(HERT_Kwt(M)*RED_MASS_EFF)*    &
               abs(log(ET))
          DES_ETAT_WALL(M) = DES_ETAT_WALL(M)/sqrt(PI*PI+(log(ET))**2)
       else
          DES_ETAT_WALL(M) = 2.d0*sqrt(HERT_Kwt(M)*RED_MASS_EFF)
       endif

       ! Calculate the collision time scale.
       TCOLL_TMP = PI/sqrt(HERT_Kwn(M)/MASS_EFF -                    &
            ((DES_ETAN_WALL(M)/MASS_EFF)**2)/4.d0)
    enddo
    

    ! If following are assigned warn user they are discarded.
    if(IS_DEFINED(KN)) then
       write(ERR_MSG, 2200) 'KN'
       call flush_err_msg
    endif
    if(IS_DEFINED(KN_W)) then
       write(ERR_MSG, 2200) 'KN_W'
       call flush_err_msg
    endif
    if(IS_DEFINED(KT_FAC)) then
       write(ERR_MSG, 2200) 'KT_FAC'
       call flush_err_msg
    endif
    if(IS_DEFINED(KT_W_FAC)) then
       write(ERR_MSG, 2200) 'KT_W_FAC'
       call flush_err_msg
    endif
    if(IS_DEFINED(DES_ETAT_FAC)) then
       write(ERR_MSG, 2200) 'DES_ETAT_FAC'
       call flush_err_msg
    endif
    if(IS_DEFINED(DES_ETAT_W_FAC)) then
       WRITE(ERR_MSG, 2200) 'DES_ETAT_W_FAC'
       call flush_err_msg
    endif
    
2200 format('Warning 2200: ',A,' values are not used ',/' with the',  &
          ' linear spring-dashpot collision model.')
    
    
    ! Store the smalled calculated collision time scale. This value is used
    ! in time-marching the DEM solids.
    DTSOLID = TCOLL/50.d0
    
    
    call finl_err_msg
    
1000 format('Error 1000: Required input not specified: ',A,/'Please ',&
          'correct the ',A,' file.')
    
1001 format('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
          'Please correct the ',A,' file.')
    
1002 format('Error 1002: Required input not specified: ',A,/          &
         'Description:',A,/'Please correct the ',A,' file.')
    
  end subroutine check_solids_dem_coll_hertz

end module check_des_solids_module
