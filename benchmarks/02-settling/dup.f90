program main

implicit none

integer :: scale

integer :: pcount, m, lc, lc1, lc2
integer :: np0, np

double precision, allocatable :: particles(:,:)
character(len=4) :: char4

write(*,"(/2x,'Reading initial particles data')")

open(file='Size0001/particle_input.dat',unit=1000)
read(1000,*) pcount
allocate(particles(pcount,8))
do lc=1,pcount
   read(1000,*)m, particles(lc,1:8)
enddo
close(1000)

write(*,"(/13x,'Size0001.  Particle Count: ',I8)") pcount

do scale=1,5
   np = 0
   np0 = pcount*(4**scale)
   write(char4,"(I4.4)") 4**scale

   write(*,"(2x,'Generating Size',A4,'.  Particle Count: ',I8)") char4, np0

   open(file='Size'//char4//'/particle_input.dat',unit=1000,status='unknown')
   write(1000,*) np0

   do lc2=0,2**scale-1
      do lc1=0,2**scale-1
         do lc=1,pcount
            write(1000,*) 1, &
               particles(lc,1)+dble(lc1)*0.0015d0, &
               particles(lc,2), &
               particles(lc,3)+dble(lc2)*0.0015d0, &
               particles(lc,4:8)
            np = np + 1
         enddo
      enddo
   enddo
   close(1000)

   if(np /= np0) then
      write(*,*) 'particle mismatch :: ',np0, np
      stop 3322
   endif
enddo

end program main
