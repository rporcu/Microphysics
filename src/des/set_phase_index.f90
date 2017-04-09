module set_phase_index_module

  use amrex_fort_module,  only: c_real => amrex_real
  use iso_c_binding ,     only: c_int

  implicit none
  private 

  public set_phase_index

contains
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: SET_PHASE_INDEX                                         !
  !                                                                      !
  !  Purpose: Set the index of a particle based on its diameter and      !
  !  density.                                                            !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine set_phase_index( phase, radius, density, state, pidx ) &
       bind(C, name="mfix_set_phase_index")

    use error_manager,  only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar
    use discretelement, only: nonexistent, normal_ghost, entering_ghost, exiting_ghost
    use param1,         only: small_number
    use constant,       only: MMAX, D_P0, RO_s0
    use run,            only: IFILE_NAME

    real(c_real),   intent(in)  :: radius, density, state
    integer(c_int), intent(in)  :: pidx
    real(c_real),   intent(out) :: phase ! Got to set it to real because of Particle Container
    integer                     :: m     ! Phase index
    integer                     :: istate! Integer phase
    real(c_real)                :: dDp   ! Diameter minus Diameter of the phase specified
    real(c_real)                :: dRho  ! Density  minus density of the phase specified

    phase = 0

    ! Proceed only for certain particle states ( do we really need this???)
    istate = int( state )
    if ( (istate == NONEXISTENT) .or. (istate == NORMAL_GHOST)  .or. &
         & (istate == ENTERING_GHOST) .or. (istate == EXITING_GHOST) )  return

    ! Determining the solids phase of each particle by matching the diameter
    ! and density to those specified in the data file.
    do m = 1, MMAX
       dDp  = abs( 2.0d0 * radius - D_P0(M) )
       dRho = abs( density - RO_S0(M) )
       if( dDp < SMALL_NUMBER .and. dRho < SMALL_NUMBER) then
          phase  = M
          exit  ! Match is found
       endif
    end do

    ! If no error (match found) return
    if ( phase > 0 ) return

    ! Point of no return: Report errors and abort
    call INIT_ERR_MSG("SET_PHASE_INDEX")
    write(ERR_MSG, '(A,/8X,A,4X,A,6X,A,/)') &
         & 'Error 1100: Unable to determine the phase of the following particle: ', 'ID', &
         & 'Diameter','Density'

    call flush_err_msg (FOOTER=.false.)

    write(ERR_MSG,'(I10,2(2x,g12.5))') pidx,  2.0 * radius, density
    call flush_err_msg(HEADER=.false., FOOTER=.false.)

    write(ERR_MSG, '(A,/,A,/8X,A,4X,A,6X,A)')  &
         & ' ','Defined phase parameters from '// trim(IFILE_NAME) // ':', &
         & 'ID', 'Diameter', 'Density'
    call flush_err_msg(HEADER=.false., FOOTER=.false.)

    do M = 1, MMAX
       write(ERR_MSG, '(I10,2(2x,g12.5))') M, D_P0(M), RO_S0(M)
       call flush_err_msg(HEADER=.false., FOOTER=.false.)
    enddo

    write(ERR_MSG, '(/,A)') 'Please correct the '//trim(IFILE_NAME) &
         & // ' or particle_input.dat files.'
    call flush_err_msg(HEADER=.false., ABORT=.true.)

  end subroutine set_phase_index
end module set_phase_index_module
