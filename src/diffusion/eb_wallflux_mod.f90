module eb_wallflux_mod

   use amrex_fort_module,  only: rt=>amrex_real, c_int
   use amrex_error_module, only: amrex_abort
   use param,              only: zero, half, one, two_thirds

   use amrex_mlebabeclap_3d_module

   implicit none

   private
   public      :: compute_diff_wallflux

contains

   !
   ! We use no-slip boundary for velocities.
   !
   subroutine compute_diff_wallflux (divw, dx, i, j, k, &
        vel, vlo, vhi,         &
        mu, slo, shi,          &
        bcent, blo, bhi,       &
        flag,  flo, fhi,       &
        apx,  axlo, axhi,      &
        apy,  aylo, ayhi,      &
        apz,  azlo, azhi,      &
        vfrac, vflo, vfhi,     &
        do_explicit_diffusion, &
        eb_ho_dirichlet)

      ! Wall divergence operator
      real(rt),       intent(  out) :: divw(3)

      ! Cell indices
      integer(c_int), intent(in   ) :: i, j, k

      ! Grid spacing
      real(rt),       intent(in   ) :: dx(3)

      ! Array bounds
      integer(c_int), intent(in   ) ::  vlo(3),  vhi(3)
      integer(c_int), intent(in   ) ::  slo(3),  shi(3)
      integer(c_int), intent(in   ) :: axlo(3), axhi(3)
      integer(c_int), intent(in   ) :: aylo(3), ayhi(3)
      integer(c_int), intent(in   ) :: azlo(3), azhi(3)
      integer(c_int), intent(in   ) :: vflo(3), vfhi(3)
      integer(c_int), intent(in   ) ::  blo(3),  bhi(3)
      integer(c_int), intent(in   ) ::  flo(3),  fhi(3)

      ! Arrays
      real(rt),       intent(in   ) ::                               &
           &   vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3),     &
           &    mu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),       &
           & bcent(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),3),     &
           &   apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3)), &
           &   apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3)), &
           &   apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3)), &
           & vfrac(vflo(1):vfhi(1),vflo(2):vfhi(2),vflo(3):vfhi(3))

      integer(c_int),  intent(in   ) :: flag(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

      ! If true  then we include all the diffusive terms in this explicit result
      ! If false then we include all only the off-diagonal terms here -- we do this
      !     by computing the full tensor then subtracting the diagonal terms
      integer(c_int),  intent(in   ) :: do_explicit_diffusion
      integer(c_int),  intent(in   ) :: eb_ho_dirichlet

      ! Local variable
      real(rt)   :: dxinv(3)
      real(rt)   :: dapx, dapy, dapz
      real(rt)   :: apnorm, apnorminv, anrmx, anrmy, anrmz
      real(rt)   :: bct(3)
      real(rt)   :: dudn, dvdn, dwdn
      real(rt)   :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, divu
      real(rt)   :: tauxx, tauyy, tauzz, tauxy, tauxz, tauyx, tauyz, tauzx, tauzy, tautmp
      real(rt)   :: phib
      integer    :: ixit, iyit, izit, is, index
      real(rt)   :: dphidn(3)

      divw  = zero
      dxinv = one / dx

      dapx = apx(i+1,j,k)-apx(i,j,k)
      dapy = apy(i,j+1,k)-apy(i,j,k)
      dapz = apz(i,j,k+1)-apz(i,j,k)

      apnorm = sqrt(dapx**2+dapy**2+dapz**2)

      if ( apnorm == zero ) then
         call amrex_abort("compute_diff_wallflux: we are in trouble.")
      end if

      apnorminv = one/apnorm
      anrmx = -dapx * apnorminv  ! unit vector pointing toward the wall
      anrmy = -dapy * apnorminv
      anrmz = -dapz * apnorminv

      ! The center of the wall
      bct = bcent(i,j,k,:)

      ! Value on wall -- here we enforce no-slip therefore 0 for all components
      phib = 0.d0

      dphidn = [dudn, dvdn, dwdn]
      do index = 1, 3
         if (eb_ho_dirichlet .eq. 1) then
            call compute_dphidn_3d_ho(dphidn(index), dxinv, i, j, k, &
               vel(:,:,:,index), vlo, vhi, &
               flag, flo, fhi, &
               bcent(i,j,k,:), phib,  &
               anrmx, anrmy, anrmz)
         else
            call compute_dphidn_3d(dphidn(index), dxinv, i, j, k, &
               vel(:,:,:,index), vlo, vhi, &
               flag, flo, fhi, &
               bcent(i,j,k,:), phib,  &
               anrmx, anrmy, anrmz, vfrac(i,j,k))
         end if
      end do

      !
      ! transform them to d/dx, d/dy and d/dz given transverse derivatives are zero
      dudx = dudn * anrmx
      dudy = dudn * anrmy
      dudz = dudn * anrmz
      !
      dvdx = dvdn * anrmx
      dvdy = dvdn * anrmy
      dvdz = dvdn * anrmz
      !
      dwdx = dwdn * anrmx
      dwdy = dwdn * anrmy
      dwdz = dwdn * anrmz

      divu = dudx+dvdy+dwdz
      tautmp = -two_thirds*mu(i,j,k)*divu  ! This MUST be verified

      tauxx = mu(i,j,k) * (dudx + dudx) + tautmp
      tauxy = mu(i,j,k) * (dudy + dvdx)
      tauxz = mu(i,j,k) * (dudz + dwdx)

      tauyx = mu(i,j,k) * (dvdx + dudy)
      tauyy = mu(i,j,k) * (dvdy + dvdy) + tautmp
      tauyz = mu(i,j,k) * (dvdz + dwdy)

      tauzx = mu(i,j,k) * (dwdx + dudz)
      tauzy = mu(i,j,k) * (dwdy + dvdz)
      tauzz = mu(i,j,k) * (dwdz + dwdz) + tautmp

      if (do_explicit_diffusion .eq. 0) then
         !
         ! Subtract diagonal terms of stress tensor, to be obtained through
         ! implicit solve instead.
         !
         tauxx = tauxx - mu(i,j,k) * dudx
         tauxy = tauxy - mu(i,j,k) * dudy
         tauxz = tauxz - mu(i,j,k) * dudz

         tauyx = tauyx - mu(i,j,k) * dvdx
         tauyy = tauyy - mu(i,j,k) * dvdy
         tauyz = tauyz - mu(i,j,k) * dvdz

         tauzx = tauzx - mu(i,j,k) * dwdx
         tauzy = tauzy - mu(i,j,k) * dwdy
         tauzz = tauzz - mu(i,j,k) * dwdz
      end if

      divw(1) = dxinv(1) * (dapx*tauxx + dapy*tauxy + dapz*tauxz)
      divw(2) = dxinv(2) * (dapx*tauyx + dapy*tauyy + dapz*tauyz)
      divw(3) = dxinv(3) * (dapx*tauzx + dapy*tauzy + dapz*tauzz)

   end subroutine compute_diff_wallflux

end module eb_wallflux_mod
