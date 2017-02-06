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
   subroutine set_bc_type(slo, shi, lo, hi, &
                          bc_ilo_type, bc_ihi_type, &
                          bc_jlo_type, bc_jhi_type, &
                          bc_klo_type, bc_khi_type, flag) &
               bind(c,name='set_bc_type')

      use bc, only: bc_defined
      use bc, only: bc_k_b, bc_k_t
      use bc, only: bc_j_s, bc_j_n
      use bc, only: bc_i_w, bc_i_e

      use param, only: dimension_bc

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      integer(c_int), intent(  out) :: bc_ilo_type&
         (slo(2):shi(2),slo(3):shi(3),2)
      integer(c_int), intent(  out) :: bc_ihi_type&
         (slo(2):shi(2),slo(3):shi(3),2)
      integer(c_int), intent(  out) :: bc_jlo_type&
         (slo(1):shi(1),slo(3):shi(3),2)
      integer(c_int), intent(  out) :: bc_jhi_type&
         (slo(1):shi(1),slo(3):shi(3),2)
      integer(c_int), intent(  out) :: bc_klo_type&
         (slo(1):shi(1),slo(2):shi(2),2)
      integer(c_int), intent(  out) :: bc_khi_type&
         (slo(1):shi(1),slo(2):shi(2),2)

      ! Local index for boundary condition
      integer :: bcv

      bc_ilo_type(slo(2):shi(2),slo(3):shi(3),1) = flag(slo(1),slo(2):shi(2),slo(3):shi(3),1)
      bc_ihi_type(slo(2):shi(2),slo(3):shi(3),1) = flag(shi(1),slo(2):shi(2),slo(3):shi(3),1)
      bc_jlo_type(slo(1):shi(1),slo(3):shi(3),1) = flag(slo(1):shi(1),slo(2),slo(3):shi(3),1)
      bc_jhi_type(slo(1):shi(1),slo(3):shi(3),1) = flag(slo(1):shi(1),shi(2),slo(3):shi(3),1)
      bc_klo_type(slo(1):shi(1),slo(2):shi(2),1) = flag(slo(1):shi(1),slo(2):shi(2),slo(3),1)
      bc_khi_type(slo(1):shi(1),slo(2):shi(2),1) = flag(slo(1):shi(1),slo(2):shi(2),shi(3),1)

      do bcv = 1, dimension_bc
         if (bc_defined(bcv)) then

            if (bc_i_w(bcv) == bc_i_e(bcv)) then
               if(bc_i_w(bcv) == slo(1)) then
                  bc_ilo_type(&
                           bc_j_s(bcv):bc_j_n(bcv),&
                           bc_k_b(bcv):bc_k_t(bcv),2) = bcv
               else if(bc_i_w(bcv) == shi(1)) then
                  bc_ihi_type(&
                           bc_j_s(bcv):bc_j_n(bcv),&
                           bc_k_b(bcv):bc_k_t(bcv),2) = bcv
               endif
            endif

            if (bc_j_s(bcv) == bc_j_n(bcv)) then
               if(bc_j_s(bcv) == slo(2)) then
                  bc_jlo_type(&
                           bc_i_w(bcv):bc_i_e(bcv),&
                           bc_k_b(bcv):bc_k_t(bcv),2) = bcv
               else if(bc_j_s(bcv) == shi(2)) then
                  bc_jhi_type(&
                           bc_i_w(bcv):bc_i_e(bcv),&
                           bc_k_b(bcv):bc_k_t(bcv),2) = bcv
               endif
            endif

            if (bc_k_b(bcv) == bc_k_t(bcv)) then
               if(bc_k_b(bcv) == slo(3)) then
                  bc_klo_type(&
                           bc_i_w(bcv):bc_i_e(bcv),&
                           bc_j_s(bcv):bc_j_n(bcv),2) = bcv
               elseif(bc_k_b(bcv) == shi(3)) then
                  bc_khi_type(&
                           bc_i_w(bcv):bc_i_e(bcv),&
                           bc_j_s(bcv):bc_j_n(bcv),2) = bcv
               endif
            endif

         endif
      enddo

   end subroutine set_bc_type

   end module set_bc_type_module
