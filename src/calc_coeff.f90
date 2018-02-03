module calc_coeff_module

      use calc_tau_u_g_module, only: calc_tau_u_g
      use calc_tau_v_g_module, only: calc_tau_v_g
      use calc_tau_w_g_module, only: calc_tau_w_g
      use calc_trd_g_module, only: calc_trd_g
      use physical_prop_module, only: physical_prop

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

  contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_COEFF                                              !
!  Purpose: This routine directs the calculation of all physical and   !
!           transport properties, and exchange rates.                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
    subroutine calc_coeff(slo, shi, lo, hi, plevel, &
      ro_g, p_g, p0_g, ep_g, rop_g) bind(C, name="calc_coeff")

      implicit none

! Dummy arguments
!-----------------------------------------------------------------------
! Level to calculate physical properties.
! 0) Only density
! 1) Everything but density
! 2) All physical properties
      integer(c_int), intent(in   ) :: slo(3), shi(3), lo(3), hi(3)

      integer(c_int), intent(in   ) :: plevel

      real(c_real), intent(in   ) ::  &
         p_g (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
        p0_g (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
        ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(inout) :: &
        ro_g (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
        rop_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Calculate physical properties: (density, specific heat, diameter)
      call physical_prop(slo, shi, lo, hi, plevel, ro_g, p_g, p0_g, ep_g, rop_g)

      end subroutine calc_coeff

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_tau_g                                              !
!  Purpose: This routine directs the calculation of all physical and   !
!           transport properties, and exchange rates.                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine calc_tau_g(slo,shi,ulo,uhi,vlo,vhi,wlo,whi,lo,hi,&
         tau_u_g, tau_v_g, tau_w_g,&
         u_g,v_g,w_g,trd_g,lambda_g,mu_g,dx,dy,dz) &
        bind(C, name="calc_tau_g")

      implicit none

      integer(c_int), intent(in ) :: slo(3),shi(3)
      integer(c_int), intent(in ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      integer(c_int), intent(in ) ::  lo(3), hi(3)

      ! Stress tensor cross terms.
      real(c_real), intent(inout) :: &
           tau_u_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           tau_v_g(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
           tau_w_g(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: &
           u_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           v_g(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
           w_g(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)), &
           trd_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           lambda_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           mu_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: dx,dy,dz

      ! Calculate the cross terms of the stress tensor (gas phase; m=0)
      call calc_tau_u_g (slo,shi,ulo,uhi,vlo,vhi,wlo,whi,lo,hi,&
         tau_u_g,trd_g,u_g,v_g,w_g,lambda_g,mu_g,dx,dy,dz)
      call calc_tau_v_g (slo,shi,ulo,uhi,vlo,vhi,wlo,whi,lo,hi,&
         tau_v_g,trd_g,u_g,v_g,w_g,lambda_g,mu_g,dx,dy,dz)
      call calc_tau_w_g (slo,shi,ulo,uhi,vlo,vhi,wlo,whi,lo,hi,&
         tau_w_g,trd_g,u_g,v_g,w_g,lambda_g,mu_g,dx,dy,dz)

    end subroutine calc_tau_g

end module calc_coeff_module
