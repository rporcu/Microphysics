MODULE WRITE_DES_DATA_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: WRITE_DES_DATA                                          !
!  Purpose: Writing DES output in Paraview format                      !
!                                                                      !
!                                                                      !
!  Author: Jay Boyalakuntla                           Date: 26-Jul-06  !
!  Reviewer: Sreekanth Pannala                        Date: 31-Oct-06  !
!                                                                      !
!  Reviewer: J. Musser                                Date: 20-Apr-10  !
!  Comments: Split original subroutine into one for ParaView *.vtp     !
!  files, and a second for TECPLOT files *.dat.                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE WRITE_DES_DATA(des_radius, des_pos_new, des_vel_new, des_usr_var)

         IMPLICIT NONE

         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: des_radius
         DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: des_vel_new, des_pos_new, des_usr_var

         CALL WRITE_DES_VTP(des_radius, des_pos_new, des_vel_new, des_usr_var)

      RETURN
      END SUBROUTINE WRITE_DES_DATA
!-----------------------------------------------

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: WRITE_DES_VTP                                           !
!  Purpose: Writing DES output in Paraview format                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE WRITE_DES_VTP(des_radius, des_pos_new, des_vel_new, des_usr_var)

      use vtp
      use discretelement, only: S_TIME
      use discretelement, only: DES_USR_VAR_SIZE

      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: des_radius
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: des_vel_new, des_pos_new, des_usr_var

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_DIAMETER !(PARTICLES)
      CHARACTER(len=10) :: lNoP
      CHARACTER(len=24) :: sTIMEc

      sTIMEc=''; WRITE(sTIMEc,"(ES24.16)") S_TIME

! This routine opens the VTP file and calculates send/recv information.
! It returns back the number of points as a string.
      CALL VTP_OPEN_FILE(lNoP)

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
      END SUBROUTINE WRITE_DES_VTP
END MODULE WRITE_DES_DATA_MODULE
