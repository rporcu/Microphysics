module set_bc_type_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc_type                                             C
!                                                                      C
!  Author: J. Musser                                  Date: 05-FEB-17  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   subroutine set_bc_type(slo, shi,&
                          bc_i_type, bc_j_type, bc_k_type, &
                          bc_i_ptr , bc_j_ptr , bc_k_ptr, flag)

      use bc, only: bc_type, bc_defined
      use bc, only: bc_k_b, bc_k_t
      use bc, only: bc_j_s, bc_j_n
      use bc, only: bc_i_w, bc_i_e

      use param, only: dimension_bc

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)

      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      integer(c_int), intent(  out) :: bc_i_type&
         (2,slo(2):shi(2),slo(3):shi(3))
      integer(c_int), intent(  out) :: bc_j_type&
         (2,slo(1):shi(1),slo(3):shi(3))
      integer(c_int), intent(  out) :: bc_k_type&
         (2,slo(1):shi(1),slo(2):shi(2))

      integer(c_int), intent(  out) :: bc_i_ptr&
         (2,slo(2):shi(2),slo(3):shi(3))
      integer(c_int), intent(  out) :: bc_j_ptr&
         (2,slo(1):shi(1),slo(3):shi(3))
      integer(c_int), intent(  out) :: bc_k_ptr&
         (2,slo(1):shi(1),slo(2):shi(2))

! local index for boundary condition
      integer, parameter :: lo_side = 1
      integer, parameter :: hi_side = 2
      integer :: lc, bcv, i, j, k, lo

      bc_i_type(1,:,:) = flag(slo(1),:,:,1)
      bc_i_type(2,:,:) = flag(shi(1),:,:,1)
      bc_j_type(1,:,:) = flag(:,slo(2),:,1)
      bc_j_type(2,:,:) = flag(:,shi(2),:,1)
      bc_k_type(1,:,:) = flag(:,:,slo(3),1)
      bc_k_type(2,:,:) = flag(:,:,shi(3),1)

      bc_i_ptr = -1
      bc_j_ptr = -1
      bc_k_ptr = -1

      do bcv = 1, dimension_bc
         if (bc_defined(bcv)) then

            if (bc_i_w(bcv) == bc_i_e(bcv)) then
               if(bc_i_w(bcv) == slo(1)) then
                  bc_i_ptr(lo_side, &
                           bc_j_s(bcv):bc_j_n(bcv),&
                           bc_k_b(bcv):bc_k_t(bcv)) = bcv
               else if(bc_i_w(bcv) == shi(1)) then
                  bc_i_ptr(hi_side, &
                           bc_j_s(bcv):bc_j_n(bcv),&
                           bc_k_b(bcv):bc_k_t(bcv)) = bcv
               endif
            endif

            if (bc_j_s(bcv) == bc_j_n(bcv)) then
               if(bc_j_s(bcv) == slo(2)) then
                  bc_j_ptr(lo_side, &
                           bc_i_w(bcv):bc_i_e(bcv),&
                           bc_k_b(bcv):bc_k_t(bcv)) = bcv
               else if(bc_j_s(bcv) == shi(2)) then
                  bc_j_ptr(hi_side, &
                           bc_i_w(bcv):bc_i_e(bcv),&
                           bc_k_b(bcv):bc_k_t(bcv)) = bcv
               endif
            endif

            if (bc_k_b(bcv) == bc_k_t(bcv)) then
               if(bc_k_b(bcv) == slo(3)) then
                  bc_k_ptr(lo_side, &
                           bc_i_w(bcv):bc_i_e(bcv),&
                           bc_j_s(bcv):bc_j_n(bcv)) = bcv
               elseif(bc_k_b(bcv) == shi(3)) then
                  bc_k_ptr(hi_side, &
                           bc_i_w(bcv):bc_i_e(bcv),&
                           bc_j_s(bcv):bc_j_n(bcv)) = bcv
               endif
            endif

         endif
      enddo

   end subroutine set_bc_type

   end module set_bc_type_module
