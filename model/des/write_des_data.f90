MODULE WRITE_DES_DATA_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: WRITE_DES_DATA                                          !
!  Purpose: Writing DES output in Paraview format                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_DES_DATA(max_pip, particle_state, des_radius, &
         des_pos_new, des_vel_new, des_usr_var)

      use vtp, only: vtp_open_file, vtp_close_file
      use vtp, only: add_vtp_to_pvd, vtp_write_element, vtp_write_data
      use discretelement, only: s_time
      use discretelement, only: des_usr_var_size

      implicit none

! Dummy arguments ....................................................//
      integer     , intent(in   ) :: max_pip
      integer     , intent(in   ) :: particle_state(max_pip)
      real(c_real), intent(in   ) :: des_radius(max_pip)
      real(c_real), intent(in   ) :: des_vel_new(max_pip,3)
      real(c_real), intent(in   ) :: des_pos_new(max_pip,3)
      real(c_real), intent(in   ) :: des_usr_var(max_pip,1)

! Local variables ....................................................//
      real(c_real), allocatable :: des_diameter(:)
      character(len=10) :: lnop
      character(len=24) :: stimec

      sTIMEc=''; WRITE(sTIMEc,"(ES24.16)") S_TIME

! This routine opens the VTP file and calculates send/recv information.
! It returns back the number of points as a string.
      CALL VTP_OPEN_FILE(max_pip, particle_state, lNoP)

! Standard VTP header information:
!----------------------------------------------------------------------/
      CALL VTP_WRITE_ELEMENT('<?xml version="1.0"?>')
      CALL VTP_WRITE_ELEMENT('<!-- Time ='//sTIMEc//'s -->')
      CALL VTP_WRITE_ELEMENT('<VTKFile type="PolyData" version="0.1" &
         &byte_order="LittleEndian">')
      CALL VTP_WRITE_ELEMENT('<PolyData>')

      CALL VTP_WRITE_ELEMENT('<Piece NumberOfPoints="'//lNoP//'" &
         &NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" &
         &NumberOfPolys="0">')

! Points are the particle identified by position:
!----------------------------------------------------------------------/
      CALL VTP_WRITE_ELEMENT('<Points>')
      CALL VTP_WRITE_DATA('Position', DES_POS_NEW)
      CALL VTP_WRITE_ELEMENT('</Points>')

! PointData are individual particle properties:
!----------------------------------------------------------------------/
      CALL VTP_WRITE_ELEMENT('<PointData Scalars="Diameter" &
         &Vectors="Velocity">')

      ALLOCATE(DES_DIAMETER(SIZE(DES_RADIUS)))
      DES_DIAMETER(:) = 2.0d0*DES_RADIUS(:)
      CALL VTP_WRITE_DATA('Diameter', DES_DIAMETER)
      DEALLOCATE(DES_DIAMETER)

      CALL VTP_WRITE_DATA('Velocity', DES_VEL_NEW)

      IF(DES_USR_VAR_SIZE > 0) &
         CALL VTP_WRITE_DATA('User Defined Var', DES_USR_VAR)

      CALL VTP_WRITE_ELEMENT('</PointData>')

! Open/Close the unused VTP tags.
!----------------------------------------------------------------------/
      CALL VTP_WRITE_ELEMENT('<CellData></CellData>')
      CALL VTP_WRITE_ELEMENT('<Verts></Verts>')
      CALL VTP_WRITE_ELEMENT('<Lines></Lines>')
      CALL VTP_WRITE_ELEMENT('<Strips></Strips>')
      CALL VTP_WRITE_ELEMENT('<Polys></Polys>')

! Close all the opened tags:
!----------------------------------------------------------------------/
      CALL VTP_WRITE_ELEMENT('</Piece>')
      CALL VTP_WRITE_ELEMENT('</PolyData>')
      CALL VTP_WRITE_ELEMENT('</VTKFile>')

      CALL VTP_CLOSE_FILE

! Add the new VTP file to the PVD file for time association.
      CALL ADD_VTP_TO_PVD

      RETURN
   END SUBROUTINE WRITE_DES_DATA
END MODULE WRITE_DES_DATA_MODULE
