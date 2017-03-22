module mod_bc

  use compar,         only: mype
  use run,            only: IFILE_NAME
  use open_files_mod, only: open_pe_log
  use ic,             only: nsw_, fsw_, psw_, fluid_ 
  use bc,             only: bc_i_w, bc_i_e, bc_j_s, bc_j_n, bc_k_b, &
                          &  bc_k_t, bc_plane
  use error_manager,  only: finl_err_msg, err_msg, flush_err_msg,   & 
                          & init_err_msg, ivar

  implicit none
  private 

  public mod_bc_i
  public mod_bc_j
  public mod_bc_k

contains
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
  !                                                                      C
  !  Subroutine: MOD_BC_I                                                C
  !  Author: P. Nicoletti                               Date: 10-DEC-91  C
  !                                                                      C
  !  Purpose: modify the "I" values for the b.c. plane                   C
  !     This is a yz plane                                               C
  !                                                                      C
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  subroutine mod_bc_i(bcv, flag, slo, shi)

    integer, intent(in) :: slo(3), shi(3), BCV
    integer, intent(in) :: flag( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3), 4 )
    integer             :: i, j, k, owner, ier, i_fluid, i_wall
    logical             :: error

    call init_err_msg("MOD_BC_I")

    k = bc_k_b(bcv)
    j = bc_j_s(bcv)
    i = bc_i_w(bcv)

    ! establish the owner of the bc
    owner = mype
    ! call global_all_sum(owner)

    if(mype == owner) then

       ! flow on west boundary (fluid cell on east).
       if(flag(i+1,j,k,1) == fluid_ .and. (&
            flag(i,j,k,1) == nsw_ .or. &
            flag(i,j,k,1) == fsw_ .or. &
            flag(i,j,k,1) == psw_)) then

          bc_plane(bcv) = 'E'

          ! flow on east boundary (fluid cell on west).
       elseif(flag(i,j,k,1) == fluid_ .and. (&
            flag(i+1,j,k,1) == nsw_ .or. &
            flag(i+1,j,k,1) == fsw_ .or. &
            flag(i+1,j,k,1) == psw_)) then

          bc_i_w(bcv) = bc_i_w(bcv) + 1
          bc_i_e(bcv) = bc_i_e(bcv) + 1
          bc_plane(bcv) = 'W'

          ! set the plane of a value we know to be wrong so we can detect the error.
       else
          bc_plane(bcv) = '.'
       endif
    endif



    ! if there is an error, send i,j,k to all ranks. report and exit.
    if(bc_plane(bcv) == '.') then
       write(err_msg, 1100) bcv, bc_i_w(bcv), bc_i_e(bcv), &
            bc_j_s(bcv), bc_k_b(bcv)
       !           flag(i,j,k,1), flag(i+1,j,k,1)
       call flush_err_msg(abort=.true.)
    endif

1100 FORMAT('Error 1100: Cannot locate flow plane for boundary ',     &
         'condition ',I3,'.',2/3x,'I West   =  ',I6,' I East   = ',I6,/&
         3x,'J South  =  ',I6,' K Bottom = ',I6)

    ! Set up the I-indices for checking the entire BC region.
    i_wall = bc_i_w(bcv)
    i_fluid = merge(i_wall-1, i_wall+1, bc_plane(bcv)=='W')


    ! first pass through all of the bc region and verify that you have
    ! inflow/outflow specified against a wall. flag any errors.
    error = .false.
    do k = bc_k_b(bcv), bc_k_t(bcv)
       do j = bc_j_s(bcv), bc_j_n(bcv)

          ! verify that the the fluid and wall cells match the flag.
          ! only check cells that you own and contain fluid.
          if(flag(i_fluid,j,k,1) /= fluid_ .and. (&
               flag(i_wall,j,k,1) /= nsw_ .or. &
               flag(i_wall,j,k,1) /= fsw_ .or. &
               flag(i_wall,j,k,1) /= psw_)) error = .true.

       enddo
    enddo

    ! Sync up the error flag across all processes.
    ! CALL GLOBAL_ALL_OR(ERROR)

    ! If an error is detected, have each rank open a log file and write
    ! it's own message. Otherwise, we need to send all the data back to
    ! PE_IO and that's too much work!
    if(error) then

       call open_pe_log(ier)

       write(err_msg, 1200) bcv
       call flush_err_msg(footer=.false.)

       do k = bc_k_b(bcv), bc_k_t(bcv)
          do j = bc_j_s(bcv), bc_j_n(bcv)

             if(flag(i_fluid,j,k,1) /= fluid_ .and. (&
                  flag(i_wall,j,k,1) /= nsw_ .or. &
                  flag(i_wall,j,k,1) /= fsw_ .or. &
                  flag(i_wall,j,k,1) /= psw_)) then

                write(err_msg, 1201) i_wall, j, k, flag(i_wall,j,k,1),  &
                     i_fluid, j, k, flag(i_fluid,j,k,1)
                call flush_err_msg(header=.false., footer=.false.)
             endif
          enddo
       enddo

       write(err_msg,"('please correct the ',A,' file.')") trim(IFILE_NAME)
       call flush_err_msg(header=.false., abort=.true.)

    endif


1200 format('Error 1200: Illegal geometry for boundary condition:',I3)

1201 format(' ',/14X,'I',7X,'J',7X,'K',7X,'IJK',4x,'FLAG',/3x,        &
         'WALL ',3(2x,I6),3x,I3,/3x,'FLUID',3(2x,I6),3x,I3)

    call finl_err_msg

  end subroutine mod_bc_i

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Module name: MOD_BC_J(BC, I_w, J_s, J_n, K_b, PLANE)                !
  !  Author: P. Nicoletti                               Date: 10-DEC-91  !
  !                                                                      !
  !  Purpose: modify the "J" values for the b.c. plane                   !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine mod_bc_j(bcv, flag, slo, shi)


    integer, intent(in) :: slo(3), shi(3), BCV 
    integer, intent(in) :: flag( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3), 4 )
    integer             :: i, j, k, owner, ier, j_fluid, j_wall
    logical             :: error


    call init_err_msg("MOD_BC_J")

    k = bc_k_b(bcv)
    j = bc_j_s(bcv)
    i = bc_i_w(bcv)

    ! establish the owner of the bc
    owner = mype
    ! call global_all_sum(owner)

    if(mype == owner)then

       if(flag(i,j+1,k,1) == fluid_ .and. (&
            flag(i,j  ,k,1) == nsw_ .or. &
            flag(i,j  ,k,1) == fsw_ .or. &
            flag(i,j  ,k,1) == psw_)) then

          bc_plane(bcv) = 'N'

       elseif(flag(i,j,k,1) == fluid_ .and. (&
            flag(i,j+1,k,1) == nsw_ .or. &
            flag(i,j+1,k,1) == fsw_ .or. &
            flag(i,j+1,k,1) == psw_)) then

          bc_j_s(bcv) = bc_j_s(bcv) + 1
          bc_j_n(bcv) = bc_j_n(bcv) + 1
          bc_plane(bcv) = 'S'

       else
          bc_plane(bcv) = '.'
       endif
    endif

    !CALL BCAST(J_S,OWNER)
    !CALL BCAST(J_N,OWNER)
    !CALL BCAST(BC_PLANE(BCV),OWNER)

    ! If there is an error, send i,j,k to all ranks. Report and exit.
    if(bc_plane(bcv) == '.') then
       write(err_msg, 1100) bcv, bc_j_s(bcv), bc_j_n(bcv), &
            bc_i_w(bcv), bc_k_b(bcv)
       call flush_err_msg(abort=.true.)
    endif


1100 FORMAT('Error 1100: Cannot locate flow plane for boundary ',     &
         'condition ',I3,'.',2/3x,'J South   =  ',I6,' J North   = ',I6,/&
         3x,'I West  =  ',I6,' K Bottom = ',I6)


    j_wall = bc_j_s(bcv)
    j_fluid = merge(j_wall-1, j_wall+1, bc_plane(bcv)=='S')


    ! First pass through all of the BC region and verify that you have
    ! inflow/outflow specified against a wall. Flag any errors.
    error = .false.
    do k = bc_k_b(bcv), bc_k_t(bcv)
       do i = bc_i_w(bcv), bc_i_e(bcv)
          if(flag(i,j_fluid,k,1) /= fluid_ .and. (&
               flag(i,j_wall,k,1) /= nsw_ .or. &
               flag(i,j_wall,k,1) /= fsw_ .or. &
               flag(i,j_wall,k,1) /= psw_)) error = .true.
       enddo
    enddo


    ! Sync up the error flag across all processes.
    ! CALL GLOBAL_ALL_OR(ERROR)

    ! If an error is detected, have each rank open a log file and write
    ! it's own message. Otherwise, we need to send all the data back to
    ! PE_IO and that's too much work!
    if(error) then

       call open_pe_log(ier)

       write(err_msg, 1200) bcv
       call flush_err_msg(footer=.false.)

1200   format('Error 1200: Illegal geometry for boundary condition:',I3)

       do k = bc_k_b(bcv), bc_k_t(bcv)
          do i = bc_i_w(bcv), bc_i_e(bcv)

             if(flag(i,j_fluid,k,1) /= fluid_ .and. (&
                  flag(i,j_wall,k,1) /= nsw_ .or. &
                  flag(i,j_wall,k,1) /= fsw_ .or. &
                  flag(i,j_wall,k,1) /= psw_)) then

                write(err_msg, 1201) &
                     i, j_wall,  k, flag(i,j_wall, k,1), &
                     i, j_fluid, k, flag(i,j_fluid,k,1)
                call flush_err_msg(header=.false., footer=.false.)
             endif

1201         format(' ',/14X,'I',7X,'J',7X,'K',7X,'FLAG',/3x,        &
                  'WALL ',3(2x,I6),3x,I3,/3x,'FLUID',3(2x,I6),3x,I3)

          enddo
       enddo

       write(err_msg,"('Please correct the ',A,' file.')") trim(IFILE_NAME)
       call flush_err_msg(header=.false., abort=.true.)

    endif

    call finl_err_msg


  end subroutine mod_bc_j


  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Module name: MOD_BC_K (BC, I_w, J_s, K_b, K_t, PLANE)               !
  !  Author: P. Nicoletti                               Date: 10-DEC-91  !
  !                                                                      !
  !  Purpose: modify the "K" values for the b.c. plane                   !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine mod_bc_k(bcv,flag,slo,shi)


    integer, intent(in) :: slo(3), shi(3), BCV
    integer, intent(in) ::flag( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3), 4 )
    integer             :: i, j, k, owner, ier, k_fluid, k_wall
    logical             :: error


    call init_err_msg("MOD_BC_K")

    k = bc_k_b(bcv)
    j = bc_j_s(bcv)
    i = bc_i_w(bcv)


    ! establish the owner of the bc
    owner = mype
    ! call global_all_sum(owner)

    if(mype == owner) then

       if(flag(i,j,k+1,1) ==fluid_ .and. (&
            flag(i,j,k,1) == nsw_ .or. &
            flag(i,j,k,1) == fsw_ .or. &
            flag(i,j,k,1) == psw_)) then
          bc_plane(bcv) = 'T'

       elseif(flag(i,j,k,1) == fluid_ .and. (&
            flag(i,j,k+1,1) == nsw_ .or. &
            flag(i,j,k+1,1) == fsw_ .or. &
            flag(i,j,k+1,1) == psw_)) then
          bc_k_b(bcv) = bc_k_b(bcv) + 1
          bc_k_t(bcv) = bc_k_t(bcv) + 1
          bc_plane(bcv) = 'B'
       else
          bc_plane(bcv) = '.'
       endif
    endif

    ! The owner distributes the new Iw/Ie coordinates to the other ranks.
    !CALL BCAST(K_B,OWNER)
    !CALL BCAST(K_T,OWNER)
    !CALL BCAST(BC_PLANE(BCV),OWNER)

    ! If there is an error, send i,j,k to all ranks. Report and exit.
    if(bc_plane(bcv) == '.') then

       write(err_msg, 1100) bcv, bc_k_b(bcv), bc_k_t(bcv), &
            bc_i_w(bcv), bc_j_s(bcv)
       call flush_err_msg(abort=.true.)
    endif

1100 format('Error 1100: Cannot locate flow plane for boundary ',     &
         'condition ',I3,'.',2/3x,'K Bottom =  ',I6,' K Top    = ',I6,/&
         3x,'I West   =  ',I6,' J South  = ',I6)

    ! Set up the I-indices for checking the entire BC region.
    k_wall = bc_k_b(bcv)
    k_fluid = merge(k_wall-1, k_wall+1, bc_plane(bcv)=='B')


    error = .false.
    do j = bc_j_s(bcv), bc_j_n(bcv)
       do i = bc_i_w(bcv), bc_i_e(bcv)

          ! only check cells that you own and contain fluid.
          if(flag(i,j,k_fluid,1) /= fluid_ .and. (&
               flag(i,j,k_wall,1) /= nsw_ .or. &
               flag(i,j,k_wall,1) /= fsw_ .or. &
               flag(i,j,k_wall,1) /= psw_)) error = .true.

       enddo
    enddo

    ! Sync up the error flag across all processes.
    ! CALL GLOBAL_ALL_OR(ERROR)
    ! If an error is detected, have each rank open a log file and write
    ! it's own message. Otherwise, we need to send all the data back to
    ! PE_IO and that's too much work!
    if(error) then

       call open_pe_log(ier)

       write(err_msg, 1200) bcv
       call flush_err_msg(footer=.false.)

1200   format('Error 1200: Illegal geometry for boundary condition:',I3)

       do j = bc_j_s(bcv), bc_j_n(bcv)
          do i = bc_i_w(bcv), bc_i_e(bcv)

             ! only check cells that you own and contain fluid.
             if(flag(i,j,k_fluid,1) /= fluid_ .and. (&
                  flag(i,j,k_wall,1) /= nsw_ .or. &
                  flag(i,j,k_wall,1) /= fsw_ .or. &
                  flag(i,j,k_wall,1) /= psw_)) then

                write(err_msg, 1201) i, j, k_wall,  flag(i,j,k_wall,1),  &
                     i, j, k_fluid, flag(i,j,k_fluid,1)
                call flush_err_msg(header=.false., footer=.false.)
             endif

1201         format(' ',/14X,'I',7X,'J',7X,'K',7X,'FLAG',/3x,        &
                  'WALL ',3(2x,I6),3x,I3,/3x,'FLUID',3(2x,I6),3x,I3)

          enddo
       enddo

       write(err_msg,"('Please correct the ',A,' file.')") trim(IFILE_NAME)
       call flush_err_msg(header=.false., abort=.true.)

    endif

    call finl_err_msg

  end subroutine mod_bc_k

end module mod_bc
