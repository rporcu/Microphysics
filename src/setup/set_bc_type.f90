module set_bc_type_module

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc_type                                             C
!                                                                      C
!  Author: J. Musser                                  Date: 05-FEB-17  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   subroutine set_bc_type(bc_ilo_type, bc_ihi_type, &
                          bc_jlo_type, bc_jhi_type, &
                          bc_klo_type, bc_khi_type, &
                          domlo, domhi, dx, dy, dz, &
                          xlength, ylength, zlength,&
                          ng) &
               bind(c,name='set_bc_type')

      use bc, only: bc_defined, bc_type, bc_plane

      use bc, only: nsw_, pinf_, pout_, minf_, ignore_
      use bc, only: undef_cell
      use bc, only: cyclic_x, cyclic_y, cyclic_z

      use param, only: dim_bc
      use param, only: equal
      use calc_cell_module, only: calc_cell_bc_flow
      use calc_cell_module, only: calc_cell_bc_wall

      implicit none

      integer(c_int), intent(in   ) :: ng
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)

      real(rt), intent(in) :: dx, dy, dz
      real(rt)  , intent(in) :: xlength, ylength, zlength

      integer(c_int), intent(inout) :: bc_ilo_type&
         (domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2)
      integer(c_int), intent(inout) :: bc_ihi_type&
         (domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2)
      integer(c_int), intent(inout) :: bc_jlo_type&
         (domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2)
      integer(c_int), intent(inout) :: bc_jhi_type&
         (domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2)
      integer(c_int), intent(inout) :: bc_klo_type&
         (domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)
      integer(c_int), intent(inout) :: bc_khi_type&
         (domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

      ! Local index for boundary condition
      integer :: type, bcv

      ! bc_ilo_type(:,:,1) = ignore_
      ! bc_ihi_type(:,:,1) = ignore_
      ! bc_jlo_type(:,:,1) = ignore_
      ! bc_jhi_type(:,:,1) = ignore_
      ! bc_klo_type(:,:,1) = ignore_
      ! bc_khi_type(:,:,1) = ignore_
      bc_ilo_type(:,:,1) = merge(undef_cell, nsw_, cyclic_x)
      bc_ihi_type(:,:,1) = merge(undef_cell, nsw_, cyclic_x)
      bc_jlo_type(:,:,1) = merge(undef_cell, nsw_, cyclic_y)
      bc_jhi_type(:,:,1) = merge(undef_cell, nsw_, cyclic_y)
      bc_klo_type(:,:,1) = merge(undef_cell, nsw_, cyclic_z)
      bc_khi_type(:,:,1) = merge(undef_cell, nsw_, cyclic_z)

      do bcv = 1, dim_bc
         if (bc_defined(bcv)) then

            select case (trim(bc_type(bcv)))
            case('IGNORE'        ,'IG' ); type = ignore_
            case('NO_SLIP_WALL'  ,'NSW'); type = ignore_
            case('P_INFLOW'      ,'PI' ); type = pinf_
            case('P_OUTFLOW'     ,'PO' ); type = pout_
            case('MASS_INFLOW'   ,'MI' ); type = minf_
            case default
               write(6,*) 'unknown bc type'
               stop 7655
            end select

            if (bc_plane(bcv) == 'E') then
               bc_ilo_type(:,:,1) = type
               bc_ilo_type(:,:,2) = bcv

            else if(bc_plane(bcv) == 'W') then
               bc_ihi_type(:,:,1) = type
               bc_ihi_type(:,:,2) = bcv

            else if(bc_plane(bcv) == 'N') then
               bc_jlo_type(:,:,1) = type
               bc_jlo_type(:,:,2) = bcv

            else if(bc_plane(bcv) == 'S') then
               bc_jhi_type(:,:,1) = type
               bc_jhi_type(:,:,2) = bcv

            else if(bc_plane(bcv) == 'T') then
               bc_klo_type(:,:,1) = type
               bc_klo_type(:,:,2) = bcv

            else if(bc_plane(bcv) == 'B') then
               bc_khi_type(:,:,1) = type
               bc_khi_type(:,:,2) = bcv

            endif

         endif
      enddo

    end subroutine set_bc_type


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc_type                                             C
!                                                                      C
!  Author: J. Musser                                  Date: 05-FEB-17  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
subroutine mfix_set_bc_mod(pID, pType, pLo, pHi, pLoc, pPg, pVel) &
     bind(c,name='mfix_set_bc_mod')

  use bc, only: bc_defined
  use bc, only: bc_type, bc_plane

  use bc, only: nsw_, pinf_, pout_, minf_, ignore_

  use bc, only: bc_center
  use bc, only: bc_normal

  use bc, only: bc_ep_g, bc_ep_s
  use bc, only: bc_p_g
  use bc, only: bc_u_g, bc_v_g, bc_w_g

  implicit none

  integer(c_int), intent(in   ) :: pID, pType
  real(rt),       intent(in   ) :: pLo(3), pHi(3), pLoc, pPg, pVel

  real(rt), parameter :: offset = 1.0d-15

  bc_ep_g(pID)   = 1.0_rt;
  bc_ep_s(pID,:) = 0.0_rt;

  select case(pID)

  case(1); bc_plane(pID) = 'E'

     bc_center(pID,1) = pLoc + offset
     bc_center(pID,2) = pLo(2) + 0.5_rt*(pHi(2) - pLo(2))
     bc_center(pID,3) = pLo(3) + 0.5_rt*(pHi(3) - pLo(3))

     bc_normal(pID,:) = (/ 1.0d0, 0.0d0, 0.0d0/)

  case(2); bc_plane(pID) = 'W'

     bc_center(pID,1) = pLoc - offset
     bc_center(pID,2) = pLo(2) + 0.5_rt*(pHi(2) - pLo(2))
     bc_center(pID,3) = pLo(3) + 0.5_rt*(pHi(3) - pLo(3))

     bc_normal(pID,:) = (/-1.0d0, 0.0d0, 0.0d0/)

  case(3); bc_plane(pID) = 'N'

     bc_center(pID,1) = pLo(1) + 0.5_rt*(pHi(1) - pLo(1))
     bc_center(pID,2) = pLoc + offset
     bc_center(pID,3) = pLo(3) + 0.5_rt*(pHi(3) - pLo(3))

     bc_normal(pID,:) = (/ 0.0d0, 1.0d0, 0.0d0/)

  case(4); bc_plane(pID) = 'S'

     bc_center(pID,1) = pLo(1) + 0.5_rt*(pHi(1) - pLo(1))
     bc_center(pID,2) = pLoc - offset
     bc_center(pID,3) = pLo(3) + 0.5_rt*(pHi(3) - pLo(3))

     bc_normal(pID,:) = (/ 0.0d0,-1.0d0, 0.0d0/)

  case(5); bc_plane(pID) = 'T'

     bc_center(pID,1) = pLo(1) + 0.5_rt*(pHi(1) - pLo(1))
     bc_center(pID,2) = pLo(2) + 0.5_rt*(pHi(2) - pLo(2))
     bc_center(pID,3) = pLoc + offset

     bc_normal(pID,:) = (/ 0.0d0, 0.0d0, 1.0d0/)

  case(6); bc_plane(pID) = 'B'

     bc_center(pID,1) = pLo(1) + 0.5_rt*(pHi(1) - pLo(1))
     bc_center(pID,2) = pLo(2) + 0.5_rt*(pHi(2) - pLo(2))
     bc_center(pID,3) = pLoc - offset

     bc_normal(pID,:) = (/ 0.0d0, 0.0d0,-1.0d0/)

  end select


  select case(pType)

  case(minf_)

     bc_type(pID) = 'MI'

     bc_p_g(pID) =   pPg;

     bc_u_g(pID) = 0.0d0;
     bc_v_g(pID) = 0.0d0;
     bc_w_g(pID) = 0.0d0;

     select case(pId)
     case(1,2); bc_u_g(pID) = pVel;
     case(3,4); bc_v_g(pID) = pVel;
     case(5,6); bc_w_g(pID) = pVel;
     end select

     bc_defined(pID) = .true.

  case(pinf_)

     bc_type(pID) = 'PI'

     bc_p_g(pID) =   pPg;

     bc_defined(pID) = .true.

  case(pout_)

     bc_type(pID) = 'PO'

     bc_p_g(pID) =   pPg;

     bc_defined(pID) = .true.

  case(nsw_)

     bc_type(pID) = 'NSW'

     bc_defined(pID) = .true.

  case(ignore_)

     bc_type(pID) = 'IG'

     bc_defined(pID) = .true.

  case default

     bc_defined(pID) = .false.

  end select


end subroutine mfix_set_bc_mod

end module set_bc_type_module
