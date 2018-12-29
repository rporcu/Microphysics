!
! Module:   interpolation_m
!
! Purpose:  collection of routines for general-purpose interpolation
!
! Author: Michele Rosso, LBNL
!
! Date:   November 15, 2018
!
module interpolation_m

   use amrex_error_module,      only: amrex_abort
   use amrex_fort_module,       only: rt => amrex_real
   use iso_c_binding,           only: c_int
   use param,                   only: half, zero, one

   implicit none
   private

   public trilinear_interp
   public trilinear_interp_eb
   public interp_stencil_is_valid


   interface trilinear_interp
      procedure trilinear_interp_single_var
      procedure trilinear_interp_double_var
   end interface trilinear_interp


   interface trilinear_interp_eb
      procedure trilinear_interp_single_var_eb
      procedure trilinear_interp_double_var_eb
   end interface trilinear_interp_eb


contains

   !
   ! Check whether interpolation stencil around x_i is valid
   ! 
   function interp_stencil_is_valid ( x_i, x_0, dx, flags, flo, fhi ) result(res)

      use amrex_ebcellflag_module, only: is_covered_cell

      real(rt),       intent(in   ) :: x_i(3), x_0(3), dx(3)
      integer(c_int), intent(in   ) :: flo(3), fhi(3)
      integer(c_int), intent(in   ) :: flags(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))     
      logical                       :: res
      integer                       :: i, j, k  

      ! Pick upper cell in the stencil
      i = floor((x_i(1) - x_0(1))/dx(1) + half)
      j = floor((x_i(2) - x_0(2))/dx(2) + half)
      k = floor((x_i(3) - x_0(3))/dx(3) + half)

      res = .not. any(is_covered_cell(flags(i-1:i,j-1:j,k-1:k)))     

   end function interp_stencil_is_valid

   !
   ! Single variable interpolation
   !
   function trilinear_interp_single_var (var, vlo, vhi, nc, x_i, x_0, dx) &
    &       result(var_i)

      ! Array bounds
      integer(c_int), intent(in   ) :: vlo(3), vhi(3)

      ! Number of components
      integer(c_int), intent(in   ) :: nc

      ! Field variable
      real(rt),       intent(in   ) :: var(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),1:nc)

      ! Coordinates of domain lower corner
      real(rt),       intent(in   ) :: x_0(3)

      ! Coordinates of interpolation point
      real(rt),       intent(in   ) :: x_i(3)

      ! Grid spacing
      real(rt),       intent(in   ) :: dx(3)

      ! Interpolated value
      real(rt)                      :: var_i(nc)

      ! Local variables
      integer                       :: i, j, k
      real(rt)                      :: odx, ody, odz
      real(rt)                      :: lx, ly, lz
      real(rt)                      :: sx_lo, sy_lo, sz_lo
      real(rt)                      :: sx_hi, sy_hi, sz_hi


      odx = one/dx(1)
      ody = one/dx(2)
      odz = one/dx(3)

      ! Pick upper cell in the stencil
      lx = (x_i(1) - x_0(1))*odx + half
      ly = (x_i(2) - x_0(2))*ody + half
      lz = (x_i(3) - x_0(3))*odz + half

      i = floor(lx)
      j = floor(ly)
      k = floor(lz)

      ! Weigths
      sx_hi = lx - i;  sx_lo = one - sx_hi
      sy_hi = ly - j;  sy_lo = one - sy_hi
      sz_hi = lz - k;  sz_lo = one - sz_hi

      var_i(1:nc) = sx_lo*sy_lo*sz_lo*var(i-1, j-1, k-1,1:nc) + &
       &            sx_lo*sy_lo*sz_hi*var(i-1, j-1, k  ,1:nc) + &
       &            sx_lo*sy_hi*sz_lo*var(i-1, j  , k-1,1:nc) + &
       &            sx_lo*sy_hi*sz_hi*var(i-1, j  , k  ,1:nc) + &
       &            sx_hi*sy_lo*sz_lo*var(i  , j-1, k-1,1:nc) + &
       &            sx_hi*sy_lo*sz_hi*var(i  , j-1, k  ,1:nc) + &
       &            sx_hi*sy_hi*sz_lo*var(i  , j  , k-1,1:nc) + &
       &            sx_hi*sy_hi*sz_hi*var(i  , j  , k  ,1:nc)


   end function trilinear_interp_single_var



   !
   ! Interpolate the sum of two variables
   !
   function trilinear_interp_double_var (var1, var2, vlo, vhi, nc, x_i, x_0, dx) &
    &       result(var_i)

      ! Array bounds
      integer(c_int), intent(in   ) :: vlo(3), vhi(3)

      ! Number of components
      integer(c_int), intent(in   ) :: nc

      ! Field variable
      real(rt),       intent(in   ) :: var1(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),1:nc)
      real(rt),       intent(in   ) :: var2(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),1:nc)

      ! Coordinates of domain lower corner
      real(rt),       intent(in   ) :: x_0(3)

      ! Coordinates of interpolation point
      real(rt),       intent(in   ) :: x_i(3)

      ! Grid spacing
      real(rt),       intent(in   ) :: dx(3)

      ! Interpolated value
      real(rt)                      :: var_i(nc)

      ! Local variables
      integer                       :: i, j, k
      real(rt)                      :: odx, ody, odz
      real(rt)                      :: lx, ly, lz
      real(rt)                      :: sx_lo, sy_lo, sz_lo
      real(rt)                      :: sx_hi, sy_hi, sz_hi


      odx = one/dx(1)
      ody = one/dx(2)
      odz = one/dx(3)

      ! Pick upper cell in the stencil
      lx = (x_i(1) - x_0(1))*odx + half
      ly = (x_i(2) - x_0(2))*ody + half
      lz = (x_i(3) - x_0(3))*odz + half

      i = floor(lx)
      j = floor(ly)
      k = floor(lz)

      ! Weigths
      sx_hi = lx - i;  sx_lo = one - sx_hi
      sy_hi = ly - j;  sy_lo = one - sy_hi
      sz_hi = lz - k;  sz_lo = one - sz_hi


      var_i(1:nc) = sx_lo*sy_lo*sz_lo*(var1(i-1, j-1, k-1,1:nc) + var2(i-1, j-1, k-1,1:nc)) + &
       &            sx_lo*sy_lo*sz_hi*(var1(i-1, j-1, k  ,1:nc) + var2(i-1, j-1, k  ,1:nc)) + &
       &            sx_lo*sy_hi*sz_lo*(var1(i-1, j  , k-1,1:nc) + var2(i-1, j  , k-1,1:nc)) + &
       &            sx_lo*sy_hi*sz_hi*(var1(i-1, j  , k  ,1:nc) + var2(i-1, j  , k  ,1:nc)) + &
       &            sx_hi*sy_lo*sz_lo*(var1(i  , j-1, k-1,1:nc) + var2(i  , j-1, k-1,1:nc)) + &
       &            sx_hi*sy_lo*sz_hi*(var1(i  , j-1, k  ,1:nc) + var2(i  , j-1, k  ,1:nc)) + &
       &            sx_hi*sy_hi*sz_lo*(var1(i  , j  , k-1,1:nc) + var2(i  , j  , k-1,1:nc)) + &
       &            sx_hi*sy_hi*sz_hi*(var1(i  , j  , k  ,1:nc) + var2(i  , j  , k  ,1:nc))


   end function trilinear_interp_double_var



   !
   ! Single variable interpolation -- EB version
   !
   function trilinear_interp_single_var_eb (var, vlo, vhi, nc, flags, flo, fhi, x_i, x_0, dx) &
    &       result(var_i)

      use amrex_ebcellflag_module, only: is_covered_cell

      ! Array bounds
      integer(c_int), intent(in   ) :: vlo(3), vhi(3)
      integer(c_int), intent(in   ) :: flo(3), fhi(3)

      ! Number of components
      integer(c_int), intent(in   ) :: nc

      ! Field variable
      real(rt),       intent(in   ) :: var(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),1:nc)

      ! EB Flags
      integer(c_int), intent(in   ) :: flags(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

      ! Coordinates of domain lower corner
      real(rt),       intent(in   ) :: x_0(3)

      ! Coordinates of interpolation point
      real(rt),       intent(in   ) :: x_i(3)

      ! Grid spacing
      real(rt),       intent(in   ) :: dx(3)

      ! Interpolated value
      real(rt)                      :: var_i(nc)

      ! Local variables
      integer                       :: i, j, k
      real(rt)                      :: odx, ody, odz
      real(rt)                      :: lx, ly, lz
      real(rt)                      :: sx_lo, sy_lo, sz_lo
      real(rt)                      :: sx_hi, sy_hi, sz_hi

      odx = one/dx(1)
      ody = one/dx(2)
      odz = one/dx(3)

      ! Pick upper cell in the stencil
      lx = (x_i(1) - x_0(1))*odx + half
      ly = (x_i(2) - x_0(2))*ody + half
      lz = (x_i(3) - x_0(3))*odz + half

      i = floor(lx)
      j = floor(ly)
      k = floor(lz)

      if ( any(is_covered_cell(flags(i-1:i,j-1:j,k-1:k))) ) then
         call amrex_abort("trilinear_interp_single_var_eb(): "// &
          & "interpolation stencil include covered cells.")
      end if

      ! Weigths
      sx_hi = lx - i;  sx_lo = one - sx_hi
      sy_hi = ly - j;  sy_lo = one - sy_hi
      sz_hi = lz - k;  sz_lo = one - sz_hi

      var_i(1:nc) = sx_lo*sy_lo*sz_lo*var(i-1, j-1, k-1,1:nc) + &
       &            sx_lo*sy_lo*sz_hi*var(i-1, j-1, k  ,1:nc) + &
       &            sx_lo*sy_hi*sz_lo*var(i-1, j  , k-1,1:nc) + &
       &            sx_lo*sy_hi*sz_hi*var(i-1, j  , k  ,1:nc) + &
       &            sx_hi*sy_lo*sz_lo*var(i  , j-1, k-1,1:nc) + &
       &            sx_hi*sy_lo*sz_hi*var(i  , j-1, k  ,1:nc) + &
       &            sx_hi*sy_hi*sz_lo*var(i  , j  , k-1,1:nc) + &
       &            sx_hi*sy_hi*sz_hi*var(i  , j  , k  ,1:nc)


   end function trilinear_interp_single_var_eb


   !
   ! Interpolate the sum of two variables -- EB version
   !
   function trilinear_interp_double_var_eb (var1, var2, vlo, vhi, nc, &
    &  flags, flo, fhi, x_i, x_0, dx) result(var_i)

      use amrex_ebcellflag_module, only: is_covered_cell

      ! Array bounds
      integer(c_int), intent(in   ) :: vlo(3), vhi(3)
      integer(c_int), intent(in   ) :: flo(3), fhi(3)

      ! Number of components
      integer(c_int), intent(in   ) :: nc

      ! Field variable
      real(rt),       intent(in   ) :: var1(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),1:nc)
      real(rt),       intent(in   ) :: var2(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),1:nc)

      ! EB Flags
      integer(c_int), intent(in   ) :: flags(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

      ! Coordinates of domain lower corner
      real(rt),       intent(in   ) :: x_0(3)

      ! Coordinates of interpolation point
      real(rt),       intent(in   ) :: x_i(3)

      ! Grid spacing
      real(rt),       intent(in   ) :: dx(3)

      ! Interpolated value
      real(rt)                      :: var_i(nc)

      ! Local variables
      integer                       :: i, j, k
      real(rt)                      :: odx, ody, odz
      real(rt)                      :: lx, ly, lz
      real(rt)                      :: sx_lo, sy_lo, sz_lo
      real(rt)                      :: sx_hi, sy_hi, sz_hi


      odx = one/dx(1)
      ody = one/dx(2)
      odz = one/dx(3)

      ! Pick upper cell in the stencil
      lx = (x_i(1) - x_0(1))*odx + half
      ly = (x_i(2) - x_0(2))*ody + half
      lz = (x_i(3) - x_0(3))*odz + half

      i = floor(lx)
      j = floor(ly)
      k = floor(lz)

      if ( any(is_covered_cell(flags(i-1:i,j-1:j,k-1:k))) ) then
         call amrex_abort("trilinear_interp_double_var_eb(): "// &
          & "interpolation stencil include covered cells.")
      end if

      ! Weigths
      sx_hi = lx - i;  sx_lo = one - sx_hi
      sy_hi = ly - j;  sy_lo = one - sy_hi
      sz_hi = lz - k;  sz_lo = one - sz_hi

      var_i(1:nc) = sx_lo*sy_lo*sz_lo*(var1(i-1, j-1, k-1,1:nc) + var2(i-1, j-1, k-1,1:nc)) + &
       &            sx_lo*sy_lo*sz_hi*(var1(i-1, j-1, k  ,1:nc) + var2(i-1, j-1, k  ,1:nc)) + &
       &            sx_lo*sy_hi*sz_lo*(var1(i-1, j  , k-1,1:nc) + var2(i-1, j  , k-1,1:nc)) + &
       &            sx_lo*sy_hi*sz_hi*(var1(i-1, j  , k  ,1:nc) + var2(i-1, j  , k  ,1:nc)) + &
       &            sx_hi*sy_lo*sz_lo*(var1(i  , j-1, k-1,1:nc) + var2(i  , j-1, k-1,1:nc)) + &
       &            sx_hi*sy_lo*sz_hi*(var1(i  , j-1, k  ,1:nc) + var2(i  , j-1, k  ,1:nc)) + &
       &            sx_hi*sy_hi*sz_lo*(var1(i  , j  , k-1,1:nc) + var2(i  , j  , k-1,1:nc)) + &
       &            sx_hi*sy_hi*sz_hi*(var1(i  , j  , k  ,1:nc) + var2(i  , j  , k  ,1:nc))

   end function trilinear_interp_double_var_eb


end module interpolation_m
