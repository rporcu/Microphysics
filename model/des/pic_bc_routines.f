!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: PIC_APPLY_WALLBC_STL                                    !
!  Author: R. Garg                                                     !
!                                                                      !
!  Purpose: Detect collisions between PIC particles and STLs.          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE PIC_APPLY_WALLBC_STL

      RETURN
      END SUBROUTINE PIC_APPLY_WALLBC_STL



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: Mass_OUTFLOW_PIC                                        !
!  Author: R. Garg                                   Date: 23-Jun-14   !
!                                                                      !
!  Purpose:  Routine to delete out of domain parcels for PIC           !
!  implementation                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE PIC_MO_BC

      RETURN
      END SUBROUTINE PIC_MO_BC




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DELETE_PARCEL                                           !
!  Author: R. Garg                                    Date: 23-Jun-14  !
!                                                                      !
!  Purpose:  Routine to delete parcel                                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DELETE_PARCEL(NP)

      END SUBROUTINE DELETE_PARCEL

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: MPPIC_MI_BC                                             C
!  Purpose:                                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE  PIC_MI_BC

      END SUBROUTINE PIC_MI_BC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PIC_FIND_EMPTY_SPOT                                   C
!  Purpose:                                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PIC_FIND_EMPTY_SPOT(LAST_INDEX, EMPTY_SPOT)
      END SUBROUTINE PIC_FIND_EMPTY_SPOT

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PIC_REFLECT_PARTICLE                                  C
!  Purpose:                                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PIC_REFLECT_PART(LL, WALL_NORM)

      END SUBROUTINE PIC_REFLECT_PART




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PIC_FIND_NEW_CELL                                     C
!  Purpose:                                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PIC_FIND_NEW_CELL(LL)
      END SUBROUTINE PIC_FIND_NEW_CELL

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: CHECK_IF_PARCEL_OVERLAPS_STL                            C
!  Authors: Rahul Garg                               Date: 21-Mar-2014 C
!                                                                      C
!  Purpose: This subroutine is special written to check if a particle  C
!          overlaps any of the STL faces. The routine exits on         C
!          detecting an overlap. It is called after initial            C
!          generation of lattice configuration to remove out of domain C
!          particles                                                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CHECK_IF_PARCEL_OVERLAPS_STL(POSITION, &
      OVERLAP_EXISTS)


      END SUBROUTINE CHECK_IF_PARCEL_OVERLAPS_STL

      SUBROUTINE write_this_facet_and_parcel(FID, position, velocity)
    end SUBROUTINE write_this_facet_and_parcel

