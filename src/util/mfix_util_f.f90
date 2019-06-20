  subroutine get_gravity ( grav_out) bind(C)

      use amrex_fort_module, only: rt => amrex_real

      use constant, only: gravity
      real(rt), intent(out)  :: grav_out(3)
      grav_out(:) = gravity(:)

  end subroutine get_gravity

  subroutine get_bc_list ( bc_list_out) bind(C)

      use iso_c_binding
      use bc, only: pinf_, pout_, minf_

      type, bind(C) :: bc_list
        integer(c_int) :: minf
        integer(c_int) :: pinf
        integer(c_int) :: pout
      end type bc_list

      type(bc_list), intent(out)  :: bc_list_out
      
      bc_list_out%minf = minf_
      bc_list_out%pinf = pinf_
      bc_list_out%pout = pout_

  end subroutine get_bc_list

  subroutine get_domain_bc ( domain_bc_out ) bind(C)

      use iso_c_binding, only: c_int
      use bc, only: bc_type, cyclic_x, cyclic_y, cyclic_z

      integer(c_int), intent(out)  :: domain_bc_out(6)

      integer :: bcv

      ! Default is that we reflect particles off domain boundaries if not periodic
      domain_bc_out(1:6) = 1
      if (cyclic_x) domain_bc_out(1:2) = 0
      if (cyclic_y) domain_bc_out(3:4) = 0
      if (cyclic_z) domain_bc_out(5:6) = 0

      do bcv = 1,6
         select case (trim(bc_type(bcv)))
           case ('P_OUTFLOW','PO')
              domain_bc_out(bcv) = 0
         end select
      end do

  end subroutine get_domain_bc
