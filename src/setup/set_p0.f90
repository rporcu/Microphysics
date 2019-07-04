!                                                                      !
!  Subroutine: set_p0                                                  !
!                                                                      !
!  Purpose: Set the pressure field inside the bed assuming gravity     !
!           is acting in the negative y-direction.                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine set_p0(lo, hi, domlo, domhi, &
 p0_g, slo, shi, &
 dx, dy, dz, xlength, ylength, zlength, delp_dir_in, &
 bct_ilo, bct_ihi, bct_jlo, bct_jhi, &
 bct_klo, bct_khi, ng) &
 bind(C, name="set_p0")

   use bc       , only: delp_x, delp_y, delp_z
   use bc       , only: dim_bc, bc_type, bc_p_g, bc_defined
   use bc       , only: pinf_, pout_, minf_
   use constant , only: gravity
   use fld_const, only: ro_g0
   use ic       , only: ic_p_g, ic_defined
   use scales   , only: scale_pressure

   use amrex_fort_module, only : ar => amrex_real
   use iso_c_binding , only: c_int
   use param   , only: zero, undefined
   use param   , only: is_defined, is_undefined
   use param, only: dim_ic

   implicit none

   integer, intent(in) ::  lo(3),  hi(3)
   integer, intent(in) :: slo(3), shi(3)
   integer, intent(in) :: domlo(3), domhi(3), ng

   real(ar), intent(inout) :: p0_g&
    (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

   real(ar), intent(in) :: dx, dy, dz
   real(ar), intent(in) :: xlength, ylength, zlength
   integer , intent(in) :: delp_dir_in

   integer(c_int), intent(in   ) :: &
    bct_ilo(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
    bct_ihi(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
    bct_jlo(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
    bct_jhi(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
    bct_klo(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
    bct_khi(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

   real(ar) :: offset = - 0.5_ar

   !-----------------------------------------------
   ! Local variables
   !-----------------------------------------------
   ! indices
   integer :: i, j, k, ibc, jbc, kbc
   integer :: icv, bcv, bcv_lo, bcv_hi
   integer :: nlft, nbot, ndwn, nrgt, ntop, nup
   integer :: delp_dir

   ! Gas pressure at the axial location j
   real(ar) :: pj, p_lo, p_hi

   ! Average pressure drop per unit length
   real(ar) :: dpodx, dpody, dpodz

   delp_dir = delp_dir_in

   nlft = max(0,domlo(1)-slo(1)+1)
   nbot = max(0,domlo(2)-slo(2)+1)
   ndwn = max(0,domlo(3)-slo(3)+1)

   nrgt = max(0,shi(1)-domhi(1))
   ntop = max(0,shi(2)-domhi(2))
   nup  = max(0,shi(3)-domhi(3))

   ! ---------------------------------------------------------------->>>
   !     If the bc's are pressure inflow/outflow then be sure to capture that in p0 
   ! ---------------------------------------------------------------->>>

   if ( (bct_ilo(domlo(2),domlo(3),1) .eq. pinf_)   .and. &
        (bct_ihi(domlo(2),domlo(3),1) .eq. pout_) ) then

      delp_dir = 0

      bcv_lo = bct_ilo(domlo(2),domlo(3),2)
      p_lo   = scale_pressure(bc_p_g(bcv_lo))

      bcv_hi = bct_ihi(domlo(2),domlo(3),2)
      p_hi   = scale_pressure(bc_p_g(bcv_hi))

      delp_x = p_lo - p_hi

      pj  = p_hi

   else if ( bct_ihi(domlo(2),domlo(3),1) .eq. pinf_  .and. &
             bct_ilo(domlo(2),domlo(3),1) .eq. pout_) then

      delp_dir = 0

      bcv_lo = bct_ilo(domlo(2),domlo(3),2)
      p_lo   = scale_pressure(bc_p_g(bcv_lo))

      bcv_hi = bct_ihi(domlo(2),domlo(3),2)
      p_hi   = scale_pressure(bc_p_g(bcv_hi))

      delp_x = p_lo - p_hi

      pj  = p_hi

   else if ( bct_jlo(domlo(1),domlo(3),1) .eq. pinf_  .and. &
             bct_jhi(domlo(1),domlo(3),1) .eq. pout_) then

      delp_dir = 1

      bcv_lo = bct_jlo(domlo(1),domlo(3),2)
      p_lo   = scale_pressure(bc_p_g(bcv_lo))

      bcv_hi = bct_jhi(domlo(1),domlo(3),2)
      p_hi   = scale_pressure(bc_p_g(bcv_hi))

      delp_y = p_lo - p_hi

      pj  = p_hi

   else if ( bct_jhi(domlo(1),domlo(3),1) .eq. pinf_  .and. &
             bct_jlo(domlo(1),domlo(3),1) .eq. pout_) then

      delp_dir = 1

      bcv_lo = bct_jlo(domlo(1),domlo(3),2)
      p_lo   = scale_pressure(bc_p_g(bcv_lo))

      bcv_hi = bct_jhi(domlo(1),domlo(3),2)
      p_hi   = scale_pressure(bc_p_g(bcv_hi))

      delp_y = p_lo - p_hi

      pj = p_hi

   else if ( bct_klo(domlo(1),domlo(2),1) .eq. pinf_  .and. &
             bct_khi(domlo(1),domlo(2),1) .eq. pout_) then

      delp_dir = 2

      bcv_lo = bct_klo(domlo(1),domlo(2),2)
      p_lo   = scale_pressure(bc_p_g(bcv_lo))

      bcv_hi = bct_khi(domlo(1),domlo(2),2)
      p_hi   = scale_pressure(bc_p_g(bcv_hi))

      delp_z = p_lo - p_hi

      pj  = p_hi


   else if ( bct_khi(domlo(1),domlo(2),1) .eq. pinf_  .and. &
             bct_klo(domlo(1),domlo(2),1) .eq. pout_) then

      delp_dir = 2

      bcv_lo = bct_klo(domlo(1),domlo(2),2)
      p_lo   = scale_pressure(bc_p_g(bcv_lo))

      bcv_hi = bct_khi(domlo(1),domlo(2),2)
      p_hi   = scale_pressure(bc_p_g(bcv_hi))

      delp_z = p_lo - p_hi

      pj = p_hi

   end if

   ! ---------------------------------------------------------------->>>
   !  Set default value of pj to zero in case no initial conditions are set
   pj = zero

   !  Make sure that ic_p_g is set if using delp pressure conditions
   do icv = 1, dim_ic
      if (ic_defined(icv)) then
         if ( (delp_dir .ge. 0) .and. (delp_dir .eq. delp_dir_in) ) then
            if (.not. is_defined(ic_p_g(icv))) then
               print *,'MUST DEFINE ic_p_g if using the DELP pressure condition'
               stop
            end if
            pj = ic_p_g(icv)

         else if ( (delp_dir .ge. 0) .and. (delp_dir .ne. delp_dir_in) ) then
            if (is_defined(ic_p_g(icv))) then
               print *,'MUST not define ic_p_g if setting p_inflow and p_outflow'
               stop
            end if

         else
            if (is_undefined(ic_p_g(icv))) goto 60
            if (gravity(1).ne.0.d0 .or. gravity(2).ne.0.d0 .or. gravity(3).ne.0.d0) goto 60
            p0_g(:,:,:) = ic_p_g(icv)
         end if
      end if
   end do

   ! ---------------------------------------------------------------->>>

   ! Here the pressure in each cell is determined from a specified pressure
   ! drop across the domain length. This section requires that the pressure
   ! is already defined in all initial condition regions (otherwise this
   ! section would be skipped)

   !  This hack allows to set the IC pressure  at L-dx/2 for both
   !  nodal and CC pressure -> reference value for pressure, AKA IC_P_G,
   !  is set at the last cell center location.
   if (delp_dir .ne. delp_dir_in) offset = -1.0_ar

   if (abs(delp_x) > epsilon(zero)) then
      dpodx = delp_x/xlength
      pj = pj - dpodx*dx*(hi(1)-domhi(1)+ng+2 + offset )
      do i = shi(1), slo(1), -1
         pj = pj + dpodx*dx
         p0_g(i,slo(2):shi(2),slo(3):shi(3)) = scale_pressure(pj)
      enddo
   endif

   if (abs(delp_y) > epsilon(zero)) then
      dpody = delp_y/ylength
      pj = pj - dpody*dy*(hi(2)-domhi(2)+ng+2 + offset )
      do j = shi(2), slo(2), -1
         pj = pj + dpody*dy
         p0_g(slo(1):shi(1),j,slo(3):shi(3)) = scale_pressure(pj)
      enddo
   endif

   if (abs(delp_z) > epsilon(zero)) then
      dpodz = delp_z/zlength
      pj = pj - dpodz*dz*(hi(3)-domhi(3)+ng+2 + offset )
      do k = shi(3), slo(3), -1
         pj = pj + dpodz*dz
         p0_g(slo(1):shi(1),slo(2):shi(2),k) = scale_pressure(pj)
      end do
   endif

   GOTO 100   ! pressure in all initial condition region cells was defined

   ! ----------------------------------------------------------------<<<

60 CONTINUE   ! pressure in an initial condition region cell was undefined

   ! ---------------------------------------------------------------->>>

   ! Search for an outflow boundary condition where pressure is specified
   pj = undefined
   do icv = 1, dim_bc
      if (bc_defined(icv)) then
         if(bc_type(icv)=='P_OUTFLOW' .or. bc_type(icv)=='PO') &
          pj = bc_p_g(icv)
      endif
   enddo

   ! Either a PO was not specified and/or a PO was specified but not the
   ! pressure at the outlet
   if (is_undefined(pj)) then
      p0_g = zero
      goto 100
   endif

   ! ----------------------------------------------------------------<<<

   ! Set an approximate pressure field assuming that the pressure drop
   ! balances the weight of the bed, if the initial pressure-field is not
   ! specified

   if (abs(gravity(1)) > epsilon(0.0d0)) then

      ! Find the average weight per unit area over an x-z slice
      dpodx = -gravity(1)*ro_g0

      if (gravity(1) <= 0.0d0) then
         do i = domhi(1)+1, domlo(1), -1
            if (i <= shi(1) .and. i >= slo(1)) &
             p0_g(i,:,:) = scale_pressure(pj)
            pj = pj + dpodx*dx
         enddo
      else
         do i = domlo(1), domhi(1)+1
            if (i <= shi(1) .and. i >= slo(1)) &
             p0_g(i,:,:) = scale_pressure(pj)
            pj = pj - dpodx*dx
         enddo
      endif

   else if (abs(gravity(2)) > epsilon(0.0d0)) then

      dpody = -gravity(2)*ro_g0

      if (gravity(2) <= 0.0d0) then
         do j = domhi(2)+1, domlo(2), -1
            if (j <= shi(2) .and. j >= slo(2)) &
             p0_g(:,j,:) = scale_pressure(pj)
            pj = pj + dpody*dy
         enddo
      else
         do j = domlo(2),domhi(2)+1
            if (j <= shi(2) .and. j >= slo(2)) &
             p0_g(:,j,:) = scale_pressure(pj)
            pj = pj - dpody*dy
         enddo
      endif

   else if (abs(gravity(3)) > epsilon(0.0d0)) then

      dpodz = -gravity(3)*ro_g0

      if(gravity(3) <= 0.0d0) then
         do k = domhi(3)+1, domlo(3), -1
            if (k <= shi(3) .and. k >= slo(3)) &
             p0_g(:,:,k) = scale_pressure(pj)
            pj = pj + dpodz*dz
         enddo
      else
         do k = domlo(3),domhi(3)+1
            if (k <= shi(3) .and. k >= slo(3)) &
             p0_g(:,:,k) = scale_pressure(pj)
            pj = pj - dpodz*dz
         enddo
      endif

   endif

   ! ----------------------------------------------------------------<<<

100 continue

   ! ---------------------------------------------------------------->>>
   if (nlft .gt. 0) then
      do k=slo(3),shi(3)
         do j=slo(2),shi(2)

            kbc = k
            jbc = j
            if (k .gt. domhi(3)+ng) kbc = kbc - 1
            if (j .gt. domhi(2)+ng) jbc = jbc - 1

            select case ( bct_ilo(jbc,kbc,1) )

            case (pinf_, pout_)

               bcv = bct_ilo(jbc,kbc,2)
               p0_g(slo(1):domlo(1)  ,j,k) = scale_pressure(bc_p_g(bcv))

            end select
         end do
      end do
   endif

   if (nrgt .gt. 0) then
      do k=slo(3),shi(3)
         do j=slo(2),shi(2)

            kbc = k
            jbc = j
            if (k .gt. domhi(3)+ng) kbc = kbc - 1
            if (j .gt. domhi(2)+ng) jbc = jbc - 1

            select case ( bct_ihi(jbc,kbc,1) )

            case (pinf_, pout_)

               bcv = bct_ihi(jbc,kbc,2)
               p0_g(domhi(1)+1:shi(1),j,k) = scale_pressure(bc_p_g(bcv))

            end select
         end do
      end do
   endif

   if (nbot .gt. 0) then
      do k=slo(3),shi(3)
         do i=slo(1),shi(1)

            kbc = k
            ibc = i
            if (k .gt. domhi(3)+ng) kbc = kbc - 1
            if (i .gt. domhi(1)+ng) ibc = ibc - 1

            select case ( bct_jlo(ibc,kbc,1) )

            case (pinf_, pout_)

               bcv = bct_jlo(ibc,kbc,2)
               p0_g(i,slo(2):domlo(2)  ,k) = scale_pressure(bc_p_g(bcv))

            end select
         end do
      end do
   endif

   if (ntop .gt. 0) then
      do k = slo(3),shi(3)
         do i = slo(1),shi(1)

            kbc = k
            ibc = i
            if (k .gt. domhi(3)+ng) kbc = kbc - 1
            if (i .gt. domhi(1)+ng) ibc = ibc - 1

            select case ( bct_jhi(ibc,kbc,1) )

            case (pinf_, pout_)

               bcv = bct_jhi(ibc,kbc,2)
               p0_g(i,domhi(2)+1:shi(2),k) = scale_pressure(bc_p_g(bcv))

            end select
         end do
      end do
   endif

   if (ndwn .gt. 0) then
      do j=slo(2),shi(2)
         do i=slo(1),shi(1)

            jbc = j
            ibc = i
            if (j .gt. domhi(2)+ng) jbc = jbc - 1
            if (i .gt. domhi(1)+ng) ibc = ibc - 1

            select case ( bct_klo(ibc,jbc,1) )

            case (pinf_, pout_)

               bcv = bct_klo(ibc,jbc,2)
               p0_g(i,j,slo(3):domlo(3)  ) = scale_pressure(bc_p_g(bcv))

            end select
         end do
      end do
   endif

   if (nup .gt. 0) then
      do j=slo(2),shi(2)
         do i=slo(1),shi(1)

            jbc = j
            ibc = i
            if (j .gt. domhi(2)+ng) jbc = jbc - 1
            if (i .gt. domhi(1)+ng) ibc = ibc - 1

            select case ( bct_khi(ibc,jbc,1) )

            case (pinf_, pout_)

               bcv = bct_khi(ibc,jbc,2)
               p0_g(i,j,domhi(3)+1:shi(3)) = scale_pressure(bc_p_g(bcv))

            end select
         end do
      end do
   endif

end subroutine set_p0
