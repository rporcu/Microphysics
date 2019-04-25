  subroutine mfix_sum_mf(lo,hi,rho,r_lo,r_hi,dx, &
                         vol,v_lo,v_hi,mass) bind(c,name='mfix_sum_mf')

      use amrex_fort_module, only: rt => amrex_real, amrex_add

      implicit none

      integer,  intent(in   ) :: lo(3), hi(3)
      integer,  intent(in   ) :: r_lo(3), r_hi(3)
      integer,  intent(in   ) :: v_lo(3), v_hi(3)
      real(rt), intent(in   ) :: dx(3)
      real(rt), intent(in   ) :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
      real(rt), intent(in   ) :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
      real(rt), intent(inout) :: mass

      integer  :: i, j, k
      real(rt) :: dm

      !$gpu

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
  
               dm = rho(i,j,k) * vol(i,j,k)
  
               call amrex_add(mass, dm)
  
            enddo
         enddo
      enddo

  end subroutine mfix_sum_mf

  subroutine get_gravity ( grav_out) bind(C)

      use amrex_fort_module, only: rt => amrex_real

      use constant, only: gravity
      real(rt), intent(out)  :: grav_out(3)
      grav_out(:) = gravity(:)

  end subroutine get_gravity

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
