program main

integer :: scale

integer :: pcount, m, lc, lc1, lc2

double precision, allocatable :: particles(:,:)

integer :: narg

character(len=256) :: fname

narg = command_argument_count()

if(narg /=1 ) stop 3322

call get_command_argument(narg, value = fname)

read(fname,*) scale

open(file='../Size001/particle_input.dat',unit=1000)

read(1000,*) pcount

write(*,*) 'particle count',pcount

allocate(particles(pcount,8))

do lc=1,pcount

   read(1000,*)m, particles(lc,1:8)

enddo

close(1000)

open(file='particle_input.dat',unit=1000,status='unknown')

write(1000,*) 2500*scale*scale

do lc2=0,(scale-1)

  do lc1=0,(scale-1)

     do lc=1,pcount

       write(1000,*) 1, &

          particles(lc,1)+dble(lc1)*0.0008d0, &

          particles(lc,2), &

          particles(lc,3)+dble(lc2)*0.0008d0, &

          particles(lc,4:8)

     enddo

  enddo

enddo

close(1000)

end program main
