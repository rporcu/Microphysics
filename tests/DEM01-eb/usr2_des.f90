!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS2_DES                                               !
!                                                                      !
!  Purpose: This routine is called within the discrete phase time loop !
!  after the source terms are applied and the time step updated. The   !
!  The user may insert code in this routine or call user defined       !
!  subroutines.                                                        !
!                                                                      !
!  This routien is called from the time loop, but no indicies (fluid   !
!  cell or particle) are defined.                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 06-Nov-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine USR2_DES( np, particle )

  use amrex_fort_module, only: c_real => amrex_real
  use particle_mod,      only: particle_t

  implicit none

  integer,          intent(in   ) :: np
  type(particle_t), intent(inout) :: particle(np)

  integer, save :: passes = 0

  character(len=4) :: cpass
  integer :: lc1, lc2


  return

  passes = passes + 1

  write(cpass,"(I4.4)") passes

  open(unit=100, file='out.'//cpass//'.vtp', status='unknown')

  ! Write the necessary header information for a PolyData file type
  write(100,"(A)")'<?xml version="1.0"?>'
  write(100,"(2A)") '<VTKFile type="PolyData"',&
       ' version="0.1" byte_order="LittleEndian">'
  write(100,"(3x,A)") '<PolyData>'

  ! Write Piece tag and identify the number of particles in the system.
  write(100,"(6x,a,i10.10,a,a)") &
       '<Piece NumberOfPoints="',np, '" NumberOfVerts="0" ', &
       'NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">'

  write(100,"(9x,a)")'<PointData>'

  write(100,"(12x,a)") '<DataArray type="Float32" Name="radius" &
       &NumberOfComponents="1" format="ascii">'
  do lc1 = 1, np
     write (100,"(15x,es13.6)") real(particle(lc1)%radius)
  end do
  write(100,"(12x,a)") '</DataArray>'

  write(100,"(12x,a)") '<DataArray type="Float32" Name="velocity" &
       &NumberOfComponents="3" format="ascii">'
  do lc1 = 1, np
     write (100,"(15x,3(es13.6,1x))") real(particle(lc1) % vel)
  end do
  write(100,"(12x,a)") '</DataArray>'
  write(100,"( 9x,a)") '</PointData>'

  write(100,"(9x,a)") '<Points>'
  write(100,"(12x,a,a)") '<DataArray type="Float32" ',&
       'Name="Position" NumberOfComponents="3" format="ascii">'
  do lc1 = 1,np
     write (100,"(15x,3(es13.6,3x))") real(particle(lc1) % pos)
  enddo
  write(100,"(12x,a,/9x,a)")'</DataArray>','</Points>'

  ! Write tags for data not included (vtp format style)
  write(100,"(9x,a,/9x,a,/9x,a,/9x,a)")'<Verts></Verts>',&
       '<Lines></Lines>','<Strips></Strips>','<Polys></Polys>'
  write(100,"(6x,a,/3x,a,/a)")&
       '</Piece>','</PolyData>','</VTKFile>'

  close(100)

  return
end subroutine usr2_des
