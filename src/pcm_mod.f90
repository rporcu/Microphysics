!
! PCM = Predictor Corrector Method
!
! This module implements the second-order predictor-corrector algorithm
! (Heun's method) in two-steps
!
! STEP 1: Using I order Euler, predict the value of the velocity field at
!         the next time step:
!
!          u_t = u^n + dt * RU(u^n)
!
!        RU  = velocity acceleration term ( convection + diffusion only )
!        u_t = predicted velocity field (WITHOUT external forcing)
!
! STEP 1a: The calling subroutine should add any external forcing, F, to
!          u_t to obtain u*
!
!         u*  = u_t + dt * F
!
!         u* = predicted velocity field (WITH external forcing)
!
! STEP 2: Evaluate RU(u*) and correct u*
!
!          u*_t = 0.5 * ( u^n + u* ) + 0.5 * dt * RU(u*) =
!               = 0.5 * ( u^n + u* + dt * RU(u*) )
!
! STEP 2a: the calling subroutine should add once again any external forcing,
!          F, to u*_t to obtain u**, the corrected velocity at the next time
!          step
!
!          u**  = u*_t + 0.5 * dt * F 
!
!  
! Author: Michele Rosso
!
! Date: October 26, 2017
!
! NOTE: this module makes use of automatic arrays for the time being.
!       If this proves to be too slow, I will replace this at a later
!       time.
!
! 
module pcm_mod
   
   use amrex_fort_module, only: ar => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one
   use convection_mod
   use diffusion_mod

   implicit none
   private

contains


   !
   ! Performs the prediction step (STEP 1) of the PCM algorithm
   ! along direction "dir" for component "comp" of the velocity
   ! field.
   !
   ! uo, vo, wo are the velocity components at the previous
   ! time step
   !
   ! ulo,uhi,vlo,vhi,wlo,whi are the array bounds for uo,vo,wo
   ! respectively
   !
   ! This routine is direction-agnostic because the array bounds
   ! for "comp", namely clo and chi, are passed explicity instead
   ! of using the array bounds for uo, vo, wo.
   ! 
   ! sl contains the slopes of comp in the 3 directions 
   ! 
   subroutine apply_pcm_prediction ( lo, hi, comp, clo, chi, sl,&
        & uo, ulo, uhi, vo, vlo, vhi, wo, wlo, whi,             &
        & mu, slo, shi, rop, dx, dt, dir ) bind(C)
      
      ! Loop bounds
      integer(c_int), intent(in   ) :: lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: clo(3), chi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3)
      integer(c_int), intent(in   ) :: vlo(3), vhi(3)
      integer(c_int), intent(in   ) :: wlo(3), whi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)

      ! Grid and time step
      real(ar),       intent(in   ) :: dx(3), dt

      ! Direction
      integer(c_int), intent(in   ) :: dir
      
      ! Arrays
      real(ar),       intent(in   ) ::                       &
           &  uo(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           &  vo(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
           &  wo(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)), &
           &  mu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           & rop(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           &  sl(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
      
      real(ar),       intent(inout) ::                       &
           & comp(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
      
      ! Local working arrays
      real(ar)                      ::                        &
           & conv(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)), &
           & diff(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))

      ! Local variables
      integer                       :: i, j, k
      
      ! Compute convection term
      select case ( dir )
      case (1)
         call compute_divuu_x ( lo, hi, uo, ulo, uhi, sl, vo, vlo, vhi, &
              & wo, wlo, whi, conv, dx )  

         ! call compute_ugradu_x ( lo, hi, uo, ulo, uhi, vo, vlo, vhi, &
         !      & wo, wlo, whi, conv, dx )  
 
         
         ! No diffusion term for the time being
         diff(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = zero
         
         ! Perform Euler step
         call apply_euler ( lo, hi, comp, uo, clo, chi, conv, diff, dt )      

      case(2)
         call compute_divuu_y ( lo, hi, uo, ulo, uhi, vo, vlo, vhi, sl, &
              & wo, wlo, whi, conv, dx )
         
         ! call compute_ugradu_y ( lo, hi, uo, ulo, uhi, vo, vlo, vhi,  &
         !      & wo, wlo, whi, conv, dx )
         
         ! No diffusion term for the time being
         diff(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = zero

         ! Perform Euler step
         call apply_euler ( lo, hi, comp, vo, clo, chi, conv, diff, dt )      

         
      case(3)         
         call compute_divuu_z ( lo, hi, uo, ulo, uhi, vo, vlo, vhi, &
              & wo, wlo, whi, sl, conv, dx )

         ! call compute_ugradu_z ( lo, hi, uo, ulo, uhi, vo, vlo, vhi,  &
         !      & wo, wlo, whi, conv, dx )
         
         
         ! No diffusion term for the time being
         diff(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = zero
         
         ! Perform Euler step
         call apply_euler ( lo, hi, comp, wo, clo, chi, conv, diff, dt )      

      end select
      

   end subroutine apply_pcm_prediction


   !
   ! Performs the correction step (STEP 2) of the PCM algorithm:
   !
   !           u** = 0.5 * ( uo + u* + dt * RU(u*) )
   !
   ! "comp" is the component of the velocity vector to solve for and must
   ! be consistent with "dir". For example, comp=v and compo=v^n for dir=2
   ! 
   ! IMPORTANT: upon entry into this routine, it is required that:
   !
   !          1) us,vs,ws be set to u*, v*,w*, namely the  the predictors
   !             computed via apply_pcm_prediction
   !
   !          2) us,vs,ws must have the halo cells up-to-date and the
   !             correct BCs must be in place
   !
   !          3) compo must be set to the old value of comp, i.e. the value
   !             at the previous time step
   ! 
   subroutine apply_pcm_correction ( lo, hi,        &
        & comp, clo, chi, compo, sl,                &
        & us, ulo, uhi, vs, vlo, vhi, ws, wlo, whi, &
        & mu, slo, shi, rop, dx, dt, dir ) bind(C)
      
      ! Loop bounds
      integer(c_int), intent(in   ) :: lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: clo(3), chi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3)
      integer(c_int), intent(in   ) :: vlo(3), vhi(3)
      integer(c_int), intent(in   ) :: wlo(3), whi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)

      ! Grid and time step
      real(ar),       intent(in   ) :: dx(3), dt

      ! Direction
      integer(c_int), intent(in   ) :: dir
      
      ! Arrays
      real(ar),       intent(in   ) ::                        &
           & compo(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)),&  
           &  us(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),  &
           &  vs(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)),  &
           &  ws(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)),  &
           &  mu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),  &
           & rop(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),  &
           &  sl(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
      
      real(ar),       intent(inout) ::                       &
           & comp(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
      
      ! Local variables
      integer                       :: i, j, k

      ! u* + dt * R(u*) is like applying the prediction on u*
      call apply_pcm_prediction ( lo, hi, comp, clo, chi, sl, &
           & us, ulo, uhi, vs, vlo, vhi, ws, wlo, whi,        &
           & mu, slo, shi, rop, dx, dt, dir ) 

      ! Compute corrected value
      ! Before the following loop, comp = u* + dt * R(u*)
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               comp(i,j,k) = half * ( comp(i,j,k) +  compo(i,j,k) )      
            end do
         end do
      end do

   end subroutine apply_pcm_correction

   
   !
   ! Compute a single Euler time integration step for the x-component
   ! of velocity:
   ! 
   !      unew = uold + dt * RU(u)  with  RU(u) = CONV(u) + DIFF(u)
   !
   subroutine apply_euler ( lo, hi, unew, uold, ulo, uhi, conv, diff, dt )

      ! Loop bounds
      integer(c_int), intent(in   ) :: lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: ulo(3), uhi(3)

      ! Time step width
      real(ar),       intent(in   ) :: dt
      
      ! Arrays
      real(ar),       intent(in   ) ::                        &
           & uold(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           & conv(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           & diff(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      real(ar),       intent(  out) ::                        &
           & unew(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      ! Local variables
      integer                       :: i, j, k
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               unew(i,j,k) = uold(i,j,k) + dt * ( -conv(i,j,k) + diff(i,j,k) )      
            end do
         end do
      end do

   end subroutine apply_euler

   


end module pcm_mod
