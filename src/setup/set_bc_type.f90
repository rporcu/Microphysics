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
      use bc, only: bclo => bc_center
      use bc, only: bchi => bc_normal

      use bc, only: nsw_, pinf_, pout_, minf_, ignore_
      use bc, only: undef_cell
      use bc, only: cyclic_x, cyclic_y, cyclic_z

      use param, only: dim_bc
      use param, only: equal, half
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
      integer :: type, bcv, ilo, jlo, klo, ihi, jhi, khi

      bc_ilo_type(:,:,1) = merge(undef_cell, nsw_, cyclic_x)
      bc_ihi_type(:,:,1) = merge(undef_cell, nsw_, cyclic_x)
      bc_jlo_type(:,:,1) = merge(undef_cell, nsw_, cyclic_y)
      bc_jhi_type(:,:,1) = merge(undef_cell, nsw_, cyclic_y)
      bc_klo_type(:,:,1) = merge(undef_cell, nsw_, cyclic_z)
      bc_khi_type(:,:,1) = merge(undef_cell, nsw_, cyclic_z)

      ! Cover the domain extents
      do bcv = 1, 6
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

      do bcv = 10, dim_bc
         if (bc_defined(bcv)) then

            type = minf_

            ilo = floor(bclo(bcv,1)/dx + half)
            jlo = floor(bclo(bcv,2)/dy + half)
            klo = floor(bclo(bcv,3)/dz + half)

            ihi = floor(bchi(bcv,1)/dx + half)-1
            jhi = floor(bchi(bcv,2)/dy + half)-1
            khi = floor(bchi(bcv,3)/dz + half)-1


            if (bc_plane(bcv) == 'E') then
               bc_ilo_type(jlo:jhi,klo:khi,1) = type
               bc_ilo_type(jlo:jhi,klo:khi,2) = bcv

            else if(bc_plane(bcv) == 'W') then
               bc_ihi_type(jlo:jhi,klo:khi,1) = type
               bc_ihi_type(jlo:jhi,klo:khi,2) = bcv

            else if(bc_plane(bcv) == 'N') then
               bc_jlo_type(ilo:ihi,klo:khi,1) = type
               bc_jlo_type(ilo:ihi,klo:khi,2) = bcv

            else if(bc_plane(bcv) == 'S') then
               bc_jhi_type(ilo:ihi,klo:khi,1) = type
               bc_jhi_type(ilo:ihi,klo:khi,2) = bcv

            else if(bc_plane(bcv) == 'T') then
               bc_klo_type(ilo:ihi,jlo:jhi,1) = type
               bc_klo_type(ilo:ihi,jlo:jhi,2) = bcv

            else if(bc_plane(bcv) == 'B') then
               bc_khi_type(ilo:ihi,jlo:jhi,1) = type
               bc_khi_type(ilo:ihi,jlo:jhi,2) = bcv

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
  use bc, only: bc_u_g, bc_v_g, bc_w_g, bc_vel_g

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

     bc_vel_g(pID,1) = bc_u_g(pID)
     bc_vel_g(pID,2) = bc_v_g(pID)
     bc_vel_g(pID,3) = bc_w_g(pID)

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


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc_type                                             C
!                                                                      C
!  Author: J. Musser                                  Date: 05-FEB-17  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
subroutine mfix_set_bc_mod_add_mi(pPlane, xLo, yLo, zLo, xHi, yHi, zHi, pPg, pVel) &
     bind(c,name='mfix_set_bc_mod_add_mi')

  use bc, only: bc_defined
  use bc, only: bc_type, bc_plane

  use bc, only: nsw_, pinf_, pout_, minf_, ignore_

  use bc, only: bclo => bc_center
  use bc, only: bchi => bc_normal

  use bc, only: bc_ep_g, bc_ep_s
  use bc, only: bc_p_g
  use bc, only: bc_u_g, bc_v_g, bc_w_g, bc_vel_g

  implicit none

  integer(c_int), intent(in   ) :: pPlane
  real(rt),       intent(in   ) :: pPg, pVel
  real(rt),       intent(in   ) :: xlo, ylo, zlo
  real(rt),       intent(in   ) :: xhi, yhi, zhi

  real(rt), parameter :: offset = 1.0d-15

  integer :: pID
  integer, save :: next = 10

  pID  = next;
  next = next + 1

  bc_defined(pID) = .true.

  bc_ep_g(pID)   = 1.0_rt;
  bc_ep_s(pID,:) = 0.0_rt;

  bcLo(pID,:) = (/xlo, ylo, zlo/)
  bcHi(pID,:) = (/xhi, yhi, zhi/)

  bc_type(pID) = 'MI'

  bc_p_g(pID) =   pPg;

  bc_u_g(pID) = 0.0d0;
  bc_v_g(pID) = 0.0d0;
  bc_w_g(pID) = 0.0d0;

  select case(pPlane)
  case(1)
     bc_plane(pID) = 'E'
     bc_u_g(pID) = pVel;
  case(2)
     bc_plane(pID) = 'W'
     bc_u_g(pID) = pVel;
  case(3)
     bc_plane(pID) = 'N'
     bc_v_g(pID) = pVel;
  case(4)
     bc_plane(pID) = 'S'
     bc_v_g(pID) = pVel;
  case(5)
     bc_plane(pID) = 'T'
     bc_w_g(pID) = pVel;
  case(6)
     bc_plane(pID) = 'B'
     bc_w_g(pID) = pVel;
  end select

  bc_vel_g(pID,1) = bc_u_g(pID)
  bc_vel_g(pID,2) = bc_v_g(pID)
  bc_vel_g(pID,3) = bc_w_g(pID)

end subroutine mfix_set_bc_mod_add_mi


end module set_bc_type_module
