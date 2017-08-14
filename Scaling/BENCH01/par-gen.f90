! Program to generate the ASCII particle_input.dat file for BENCH01
!
! Usage:
!
!   pargen --size INTEGER --enable-vtp
!
!  --enable-vtp    write vtp output file for visualization
!
module data

  logical :: enable_vtp = .false.

  ! base domain length and particle count
  double precision, parameter :: length = 0.004
  integer,          parameter :: pcount = 1222

  ! Particle diameter, radius density
  double precision, parameter :: dp = 100.0d-6 ! (m)
  double precision, parameter :: rp = dp*0.5d0 ! (m)
  double precision, parameter :: rho = 1.0d3   ! (kg/m^3)

  ! Initial granular energy
  double precision, parameter :: T0 = 1.0d-1   ! (m^2/sec^2)

  double precision, parameter :: pi = 4.0d0*atan(1.0d0)

contains

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
subroutine new_pos(lo, hi, pc, problo, probhi, dx, np, pos, pinc, pbin, fails)

  implicit none

  integer,          intent(in   ) :: lo, hi, pc
  double precision, intent(in   ) :: problo, probhi, dx

  double precision, intent(inout) :: pos(pc,3)

  integer,          intent(inout) :: np, fails
  integer,          intent(inout) :: pinc(hi,hi,hi)

  integer, allocatable, intent(inout) :: pbin(:,:,:,:)


  double precision :: Oodx

  integer :: i, j, k, pi, pj, pk
  integer :: l, ll, overlaps
  double precision :: min_dist, dist, lpos(3), rand3(3)

  integer :: ob, nb, tenp
  integer, allocatable :: int_tmp(:,:,:,:)

  Oodx = 1.0d0/dx
  min_dist = (1.05d0*dp)**2

  do

     call random_number(rand3)
     lpos(:) = problo + probhi*rand3(:)

     ! Grid containing the new particle
     pi = floor(lpos(1)*Oodx)+1
     pj = floor(lpos(2)*Oodx)+1
     pk = floor(lpos(3)*Oodx)+1

! Local grid search for collisions.
     overlaps=0
     do k=max(lo,pk-1), min(pk+1,hi)
        do i=max(lo,pi-1), min(pi+1,hi)
           do j=max(lo,pj-1), min(pj+1,hi)

              do l=1, pinc(i,j,k)

                 ll = pbin(i,j,k,l)

                 dist=(pos(ll,1) - lpos(1))**2 + &
                      (pos(ll,2) - lpos(2))**2 + &
                      (pos(ll,3) - lpos(3))**2

                 if(dist < min_dist) overlaps = overlaps+1

              enddo
           enddo
        enddo
     enddo
     if(overlaps == 0) then
        exit
     else
        fails = fails + 1
        if(fails > 10*pc) then
           write(*,*)'Max seed failures exceeded!'
           stop 3223
        endif
     endif
  enddo

  np = np + 1

  pos(np,:) = lpos(:)
  pinc(pi,pj,pk) = pinc(pi,pj,pk) + 1
  pbin(pi,pj,pk,pinc(pi,pj,pk)) = np

  ob = ubound(pbin,4)

  if(pinc(pi,pj,pk) + 1 >= ob) then
     nb = ob + 2
     ! write(*,"('Growing pbin to ',I2)") nb
     allocate(int_tmp(hi,hi,hi,nb))
     int_tmp(:,:,:,1:ob) = pbin(:,:,:,1:ob)
     call move_alloc(int_tmp, pbin)
  endif

  if((mod(np, pc/10) == 0 .and. np < pc*0.95) .or. np==pc) &
       write(*,"(2x,'Seeded: ',I9,3x,'(',f5.0,'%)')") np,100*dble(np)/pc

end subroutine new_pos

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
subroutine scale_vel(pc, meanVel, vel)

  implicit none

  integer, intent(in   ) :: pc
  double precision, intent(inout) :: meanVel(3), vel(pc,3)

  integer :: lc1
  double precision :: gTemp

  ! Calc the average mean velocity in each direction
  meanVel = meanVel/dble(pc)

  ! Subtract mean velocity from the random velocities to get a zero
  ! mean velocity. Also, calculate the mean granular temperature.
  gTemp = 0.0d0
  do lc1 = 1, pc
     vel(lc1,:) = vel(lc1,:) - meanVel
     gTemp = gTemp + dot_product (vel(lc1,:),vel(lc1,:))
  enddo
  gtemp = gtemp/(3.0d0*dble(pc))

  ! Scale velocities so the mean granular temperature is equal to
  ! the targeted valued.
  gTemp = dsqrt(T0/gTemp)

  do lc1 = 1, pc
     vel(lc1,:) = vel(lc1,:)*gTemp
  enddo

end subroutine scale_vel


!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
subroutine write_dat(size, pc, pos, vel)

  implicit none

  integer, intent(in   ) :: size, pc

  double precision, intent(in   ) :: pos(pc,3)
  double precision, intent(in   ) :: vel(pc,3)

  integer :: lc1
  integer, parameter :: fUnit = 2227
  character(len=4) :: char4

  write(char4,"(I4.4)") size**3

  open(file='Size'//char4//'/particle_input.dat',unit=funit,status='unknown')
  write(fUnit,"(i9)") pc
  do lc1 = 1, pc
     write(fUnit,"(i1,8(1x,es13.6))") &
          1, pos(lc1,:), rp, rho, vel(lc1,:)
  enddo
  close(funit)
end subroutine write_dat


!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
subroutine write_vtp(pc, pos, vel)

  implicit none

  integer, intent(in) :: pc
  double precision, intent(in   ) :: pos(pc,3)
  double precision, intent(in   ) :: vel(pc,3)

  integer :: lc1, lc2

  open(unit=100, file='setup.vtp', status='unknown')

! Write the necessary header information for a PolyData file type
  write(100,"(A)")'<?xml version="1.0"?>'
  write(100,"(2A)") '<VTKFile type="PolyData"',&
       ' version="0.1" byte_order="LittleEndian">'
  write(100,"(3x,A)") '<PolyData>'

! Write Piece tag and identify the number of particles in the system.
  write(100,"(6x,a,i10.10,a,a)") &
       '<Piece NumberOfPoints="',pc, '" NumberOfVerts="0" ', &
       'NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">'

  write(100,"(9x,a)")'<PointData>'
  write(100,"(12x,a)") '<DataArray type="Float32" Name="Diameter" &
       &NumberOfComponents="1" format="ascii">'
  do lc1 = 1, pc
     write (100,"(15x,es13.6)") real(dp)
  end do
  write(100,"(12x,a)") '</DataArray>'

  write(100,"(12x,a)") '<DataArray type="Float32" Name="Velocity" &
       &NumberOfComponents="3" format="ascii">'
  do lc1 = 1, pc
     write (100,"(15x,3(es13.6,3x))") (real(vel(lc1,lc2)),lc2=1,3)
  end do
  write(100,"(12x,a)") '</DataArray>'
  write(100,"( 9x,a)") '</PointData>'

  write(100,"(9x,a)") '<Points>'
  write(100,"(12x,a,a)") '<DataArray type="Float32" ',&
       'Name="Position" NumberOfComponents="3" format="ascii">'
  do lc1 = 1,pc
     write (100,"(15x,3(es13.6,3x))") (real(pos(lc1,lc2)),lc2=1,3)
  enddo
  write(100,"(12x,a,/9x,a)")'</DataArray>','</Points>'

! Write tags for data not included (vtp format style)
  write(100,"(9x,a,/9x,a,/9x,a,/9x,a)")'<Verts></Verts>',&
       '<Lines></Lines>','<Strips></Strips>','<Polys></Polys>'
  write(100,"(6x,a,/3x,a,/a)")&
       '</Piece>','</PolyData>','</VTKFile>'

  close(100)

end subroutine write_vtp


!-----------------------------------------------------------------------!
!                                                                       !
!                                                                       !
!                                                                       !
!-----------------------------------------------------------------------!
subroutine read_inputs ()

  implicit none

  integer        :: length
  character(512) :: val1, val2

  integer :: lc, nargs, ios, cbrt

  nargs = command_argument_count()

  lc = 1
  do while(lc <= nargs)

     call get_command_argument(lc, val1, length); lc = lc+1

     select case (trim (val1))
     case ( "--enable-vtp" )
        enable_vtp = .true.

     case default
        write(*,*) "Option "//trim(val1)//" not recognized"
        stop 2299
     end select

  end do

end subroutine read_inputs

end module data


!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
program main

  use data

  implicit none

  integer :: lc
  integer, parameter :: sizes(6) = (/1,2,3,4,6,10/)

  call read_inputs()

  do lc=1,6
     call set_case(sizes(lc))
  enddo

contains

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
   subroutine set_case(size)

  use data

  implicit none

  integer, intent(in) :: size

  integer :: lo, hi, np, fails, pc
  double precision :: problo, probhi, dx, meanVel(3)
  double precision, allocatable :: pos(:,:), vel(:,:)

  integer, allocatable :: pinc(:,:,:), pbin(:,:,:,:)

  double precision :: ru(4)
  character(len=4) :: char4

  write(char4,"(I4.4)") size**3

  lo =  1
  hi = 10*size

  pc = pcount*(size**3)

  problo = 0.00d0
  probhi = length * dble(size)
  dx = probhi/dble(hi)

  write(*,"(/2x,'Size',A4,/2x,36('-'))")char4
  write(*,"( 2x,'Domain length:  ',f9.3)") probhi
  write(*,"( 2x,'Particle count: ',i9  )") pc

  if(allocated(pinc)) deallocate(pinc)
  if(allocated(pbin)) deallocate(pbin)
  if(allocated(pos )) deallocate(pos )
  if(allocated(vel )) deallocate(vel )

  allocate(pinc(hi,hi,hi))
  allocate(pbin(hi,hi,hi,2))
  allocate(pos(pc,3))
  allocate(vel(pc,3))

  write(*,"(/2x,'Generating particles.')")
  pinc = 0
  np = 0
  fails = 0
  do while(np < pc)

     call new_pos(lo, hi, pc, problo, probhi, dx, np, pos, pinc, pbin, fails)

     call random_number(ru)

     vel(np,1) = dsqrt(-2.0d0*dlog(ru(1)))*cos(2.0d0*pi*ru(2))
     vel(np,2) = dsqrt(-2.0d0*dlog(ru(1)))*sin(2.0d0*pi*ru(2))
     vel(np,3) = dsqrt(-2.0d0*dlog(ru(3)))*cos(2.0d0*pi*ru(4))

     meanVel = meanVel + vel(np,:)

  enddo

  write(*,"(/2x,'Scaling velocities.')")
  call scale_vel(pc, meanVel, vel)

  write(*,"( 2x,'Writing dat file.')")
  call write_dat(size, pc, pos, vel)

  if(enable_vtp) then
     write(*,"( 2x,'Writing vtp file.')")
     call write_vtp(pc, pos, vel)
  endif
end subroutine set_case

end program main
