      MODULE time_cpu

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

! cpu time/second
      real(c_real) :: CPUos
! old cpu time and time for calculating CPUos
      real(c_real) :: CPU_NLOG, TIME_NLOG
! Initial value of CPU time.
      real(c_real) :: CPU0
! Time for IO
      real(c_real) :: CPU_IO = 0.0d0

! Initial value of CPU time at the begin of MFIX, prior any I/O
      real(c_real) :: CPU00
      real(c_real) :: WALL0

! Time at start of simulation
      real(c_real) :: TIME_START
! Wall time at the beginning
      real(c_real) :: WALL_START

      END MODULE time_cpu
