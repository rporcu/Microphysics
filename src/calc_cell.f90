MODULE CALC_CELL_MODULE

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine name: CALC_CELL                                          !
!  Purpose: calculate the i, j or k cell index for the corresponding   !
!     x y or z reactor location. the index returned depends on which   !
!     half of the i, j or k cell that the x, y, or z position          !
!     intersects                                                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   pure integer function calc_cell(location, dx)

      use param1, only : half

      implicit none

      ! the x, y or z location
      real(c_real), intent(in) :: location

      ! the cell lengths along the corresponding axis (dx, dy or dz)
      real(c_real), intent(in) :: dx

      calc_cell = floor(location/dx + half) - 1

   end function calc_cell

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_cell_bc_wall                                       !
!  Purpose: calculate the i, j or k cell index for wall BCs.           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine calc_cell_bc_wall(bcv, domlo, domhi, xlength, ylength, &
      zlength, dx, dy, dz, i_w, i_e, j_s, j_n, k_b, k_t)

      use bc, only: bc_x_e, bc_x_w
      use bc, only: bc_y_n, bc_y_s
      use bc, only: bc_z_t, bc_z_b

      use param1, only: zero, equal

      implicit none

      integer,      intent(in   ) :: bcv, domlo(3), domhi(3)
      real(c_real), intent(in   ) :: dx, dy, dz
      real(c_real), intent(in   ) :: xlength, ylength, zlength

      integer,      intent(  out) :: i_w, j_s, k_b, i_e, j_n, k_t

      i_w = calc_cell (bc_x_w(bcv), dx) + 1
      i_e = calc_cell (bc_x_e(bcv), dx)
      if(equal(bc_x_w(bcv), bc_x_e(bcv))) then
         if(equal(bc_x_w(bcv),0.0d0)) then
            i_w = domlo(1)-1
            i_e = domlo(1)-1
         elseif(equal(bc_x_w(bcv),xlength)) then
            i_w = domhi(1)+1
            i_e = domhi(1)+1
         endif
      endif

      j_s = calc_cell (bc_y_s(bcv), dy) + 1
      j_n = calc_cell (bc_y_n(bcv), dy)
      if(equal(bc_y_s(bcv), bc_y_n(bcv))) then
         if(equal(bc_y_s(bcv),zero)) then
            j_s = domlo(2)-1
            j_n = domlo(2)-1
         else if (equal(bc_y_s(bcv),ylength)) then
            j_s = domhi(2)+1
            j_n = domhi(2)+1
         endif
      endif

      k_b = calc_cell (bc_z_b(bcv), dz) + 1
      k_t = calc_cell (bc_z_t(bcv), dz)
      if(equal(bc_z_b(bcv), bc_z_t(bcv))) then
         if(equal(bc_z_b(bcv),zero)) then
            k_b = domlo(3)-1
            k_t = domlo(3)-1
         elseif(equal(bc_z_b(bcv),zlength)) then
            k_b = domhi(3)+1
            k_t = domhi(3)+1
         endif
      endif

   end subroutine calc_cell_bc_wall

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_cell_bc_wall                                       !
!  Purpose: calculate the i, j or k cell index for wall BCs.           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine calc_cell_bc_flow(bcv, xlength, ylength, zlength, &
      dx, dy, dz, i_w, i_e, j_s, j_n, k_b, k_t)

      use bc, only: bc_x_e, bc_x_w
      use bc, only: bc_y_n, bc_y_s
      use bc, only: bc_z_t, bc_z_b

      use param1, only: equal

      implicit none

      integer,      intent(in   ) :: bcv
      real(c_real), intent(in   ) :: dx, dy, dz
      real(c_real), intent(in   ) :: xlength, ylength, zlength

      integer,      intent(  out) :: i_w, j_s, k_b, i_e, j_n, k_t

      i_w = calc_cell (bc_x_w(bcv), dx)
      i_e = calc_cell (bc_x_e(bcv), dx)
      if (.not.equal(bc_x_w(bcv), bc_x_e(bcv))) then
         i_w = i_w + 1
      else if(equal(bc_x_w(bcv),xlength)) then
         i_w = i_w + 1
         i_e = i_w
      endif

      j_s = calc_cell (bc_y_s(bcv), dy)
      j_n = calc_cell (bc_y_n(bcv), dy)
      if(.not.equal(bc_y_s(bcv), bc_y_n(bcv))) then
         j_s = j_s + 1
      else if(equal(bc_y_s(bcv),ylength)) then
         j_s = j_s + 1
         j_n = j_s
      endif

      k_b = calc_cell (bc_z_b(bcv), dz)
      k_t = calc_cell (bc_z_t(bcv), dz)
      if(.not.equal(bc_z_b(bcv), bc_z_t(bcv))) then
         k_b = k_b + 1
      else if(equal(bc_z_b(bcv),zlength)) then
         k_b = k_b + 1
         k_t = k_b
      endif

   end subroutine calc_cell_bc_flow



END MODULE CALC_CELL_MODULE
