module eb_wallflux_mod
   
   use amrex_fort_module,  only: rt=>amrex_real, c_int
   use amrex_error_module, only: amrex_abort
   use param,              only: zero, half, one
   
   implicit none
   
   private
   public      :: compute_diff_wallflux

contains

   !
   ! We use no-slip boundary for velocities.
   !
   subroutine compute_diff_wallflux (divw, dx, i, j, k, &
        vel, vlo, vhi,     &
        lam, mu, slo, shi, &
        bcent, blo, bhi,   &
        apx, axlo, axhi,   &
        apy, aylo, ayhi,   &
        apz, azlo, azhi,   &
        do_explicit_diffusion, &
        covered_val)

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
      integer(c_int), intent(in   ) ::  blo(3),  bhi(3)

      ! Arrays
      real(rt),       intent(in   ) ::                               &
           &   vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3),     &
           &   lam(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),       &
           &    mu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),       &
           & bcent(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),3),     &
           & apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3)),   &
           & apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3)),   &
           & apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))

      ! If true  then we include all the diffusive terms in this explicit result
      ! If false then we include all only the off-diagonal terms here -- we do this
      !     by computing the full tensor then subtracting the diagonal terms
      integer(c_int),  intent(in   ) :: do_explicit_diffusion

      real(rt),        intent(in   ) :: covered_val

      ! Local variable
      real(rt)   :: dxinv(3)
      real(rt)   :: dapx, dapy, dapz
      real(rt)   :: apnorm, apnorminv, anrmx, anrmy, anrmz
      real(rt)   :: xit, yit, zit, s
      real(rt)   :: bct(3), d1, d2, ddinv
      real(rt)   :: u1, v1, w1, u2, v2, w2, dudn, dvdn, dwdn
      real(rt)   :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, divu
      real(rt)   :: tauxx, tauyy, tauzz, tauxy, tauxz, tauyx, tauyz, tauzx, tauzy, tautmp
      integer    :: ixit, iyit, izit, is
      
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

      if (abs(anrmx).ge.abs(anrmy) .and. abs(anrmx).ge.abs(anrmz)) then
         ! y-z plane: x = const
         ! the equation for the line:  x = bct(1) - d*anrmx
         !                             y = bct(2) - d*anrmy
         !                             z = bct(3) - d*anrmz
         s = sign(one,-anrmx)
         is = nint(s)

         !
         ! the line intersects the y-z plane (x = s) at ...
         !
         d1 = (bct(1) - s) * (one/anrmx)  ! this is also the distance from wall to intersection
         yit = bct(2) - d1*anrmy
         zit = bct(3) - d1*anrmz
         iyit = j + nint(yit)
         izit = k + nint(zit)
         yit = yit - nint(yit)  ! shift so that the center of the nine cells are (0.,0.)
         zit = zit - nint(zit)
         !
         !
         ! interpolation
         call interp2d(u1,v1,w1,yit,zit,vel(i+is,iyit-1:iyit+1,izit-1:izit+1,1:3),covered_val)

         !
         ! the line intersects the y-z plane (x = 2*s) at ...
         !
         d2 = (bct(1) - 2.d0*s) * (one/anrmx)
         yit = bct(2) - d2*anrmy
         zit = bct(3) - d2*anrmz
         iyit = j + nint(yit)
         izit = k + nint(zit)
         yit = yit - nint(yit)  ! shift so that the center of the nine cells are (0.,0.)
         zit = zit - nint(zit)
         !
         ! interpolation
         call interp2d(u2,v2,w2,yit,zit,vel(i+2*is,iyit-1:iyit+1,izit-1:izit+1,1:3),covered_val)

      else if (abs(anrmy).ge.abs(anrmx) .and. abs(anrmy).ge.abs(anrmz)) then
         ! z-x plane
         s = sign(one,-anrmy)
         is = nint(s)

         d1 = (bct(2) - s) * (one/anrmy)
         xit = bct(1) - d1*anrmx
         zit = bct(3) - d1*anrmz
         ixit = i + nint(xit)
         izit = k + nint(zit)
         xit = xit - nint(xit)
         zit = zit - nint(zit)

         call interp2d(u1,v1,w1,xit,zit,vel(ixit-1:ixit+1,j+is,izit-1:izit+1,1:3),covered_val)

         d2 = (bct(2) - 2.d0*s) * (one/anrmy)
         xit = bct(1) - d2*anrmx
         zit = bct(3) - d2*anrmz
         ixit = i + nint(xit)
         izit = k + nint(zit)
         xit = xit - nint(xit)
         zit = zit - nint(zit)

         call interp2d(u2,v2,w2,xit,zit,vel(ixit-1:ixit+1,j+2*is,izit-1:izit+1,1:3),covered_val)

      else
         ! x-y plane
         s = sign(one,-anrmz)
         is = nint(s)

         d1 = (bct(3) - s) * (one/anrmz)
         xit = bct(1) - d1*anrmx
         yit = bct(2) - d1*anrmy
         ixit = i + nint(xit)
         iyit = j + nint(yit)
         xit = xit - nint(xit)
         yit = yit - nint(yit)

         call interp2d(u1,v1,w1,xit,yit,vel(ixit-1:ixit+1,iyit-1:iyit+1,k+is,1:3),covered_val)

         d2 = (bct(3) - 2.d0*s) * (one/anrmz)
         xit = bct(1) - d2*anrmx
         yit = bct(2) - d2*anrmy
         ixit = i + nint(xit)
         iyit = j + nint(yit)
         xit = xit - nint(xit)
         yit = yit - nint(yit)

         call interp2d(u2,v2,w2,xit,yit,vel(ixit-1:ixit+1,iyit-1:iyit+1,k+2*is,1:3),covered_val)

      end if

      !
      ! compute derivatives on the wall given that velocity is zero on the wall.
      !
      ddinv = one/(d1*d2*(d2-d1))
      dudn = -ddinv*(d2*d2*u1-d1*d1*u2)  ! note that the normal vector points toward the wall
      dvdn = -ddinv*(d2*d2*v1-d1*d1*v2)
      dwdn = -ddinv*(d2*d2*w1-d1*d1*w2)
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
      tautmp = lam(i,j,k)*divu  ! This MUST be verified

      if (do_explicit_diffusion .eq. 1) then
         !
         ! Stress tensor is fully explicit 
         !
         tauxx = mu(i,j,k) * (dudx + dudx) + tautmp
         tauxy = mu(i,j,k) * (dudy + dvdx)
         tauxz = mu(i,j,k) * (dudz + dwdx)

         tauyx = mu(i,j,k) * (dvdx + dudy)
         tauyy = mu(i,j,k) * (dvdy + dvdy) + tautmp
         tauyz = mu(i,j,k) * (dvdz + dwdy)

         tauzx = mu(i,j,k) * (dwdx + dudz)
         tauzy = mu(i,j,k) * (dwdy + dvdz)
         tauzz = mu(i,j,k) * (dwdz + dwdz) + tautmp

      else 
         !
         ! Stress tensor is partially explicit and partically implicit
         ! Here we neglect diagonal terms of stress tensor -- they will be 
         ! obtained through implicit solve.                   
         !
         tauxx = mu(i,j,k) * dudx + tautmp
         tauxy = mu(i,j,k) * dvdx
         tauxz = mu(i,j,k) * dwdx

         tauyx = mu(i,j,k) * dudy
         tauyy = mu(i,j,k) * dvdy + tautmp
         tauyz = mu(i,j,k) * dwdy

         tauzx = mu(i,j,k) * dudz
         tauzy = mu(i,j,k) * dvdz
         tauzz = mu(i,j,k) * dwdz + tautmp
      end if

      divw(1) = dxinv(1) * (dapx*tauxx + dapy*tauyx + dapz*tauzx)
      divw(2) = dxinv(2) * (dapx*tauxy + dapy*tauyy + dapz*tauzy)
      divw(3) = dxinv(3) * (dapx*tauxz + dapy*tauyz + dapz*tauzz)

   end subroutine compute_diff_wallflux

   subroutine interp2d(u1,v1,w1,yit,zit,v,covered_val)

      real(rt), intent(in   ) :: yit,zit,v(3,3,3),covered_val
      real(rt), intent(  out) :: u1,v1,w1

      real(rt)             :: cym,cy0,cyp,czm,cz0,czp
      real(rt)             :: val(3,3),val_c(3,3),val_l(3,3),val_r(3,3)
      real(rt)             :: interp2d_l(3), interp2d_c(3), interp2d_r(3)
      real(rt)             :: covered_eps

      covered_eps = 1.d-10 * covered_val

      ! Coefficents for quadratic interpolation
      cym = half*yit*(yit-one)
      cy0 = one-yit*yit
      cyp = half*yit*(yit+one)

      val_c(1:3,1:3) =  cym*v(1,1:3,1:3) +       cy0*v(2,1:3,1:3) + cyp*v(3,1:3,1:3)
      val_l(1:3,1:3) = -yit*v(1,1:3,1:3) + (one+yit)*v(2,1:3,1:3)
      val_r(1:3,1:3) =                     (one-yit)*v(2,1:3,1:3) + yit*v(3,1:3,1:3)

      if (any(abs(val_l(1:3,1)) .gt. covered_eps)) then
         val(1:3,1:3) = val_r(1:3,1:3)
      else if (any(abs(val_r(1:3,1)) .gt. covered_eps)) then
         val(1:3,1:3) = val_l(1:3,1:3)
      else
         val(1:3,1:3) = val_c(1:3,1:3)
      end if

      ! Coefficents for quadratic interpolation
      czm = half*zit*(zit-one)
      cz0 = one-zit*zit
      czp = half*zit*(zit+one)

      interp2d_c(1:3) =  czm*val(1,1:3) +       cz0*val(2,1:3) + czp*val(3,1:3)
      interp2d_l(1:3) = -yit*val(1,1:3) + (one+yit)*val(2,1:3)
      interp2d_r(1:3) =                   (one-yit)*val(2,1:3) + yit*val(3,1:3)

      if (abs(interp2d_l(1)) .gt. covered_eps) then
         u1 = interp2d_r(1)
         v1 = interp2d_r(2)
         w1 = interp2d_r(3)
      else if (abs(interp2d_r(1)) .gt. covered_eps) then
         u1 = interp2d_l(1)
         v1 = interp2d_l(2)
         w1 = interp2d_l(3)
      else
         u1 = interp2d_c(1)
         v1 = interp2d_c(2)
         w1 = interp2d_c(3)
      end if

      ! interp2d = czm*(cym*v(1,1) + cy0*v(2,1) + cyp*v(3,1)) &
      !      +     cz0*(cym*v(1,2) + cy0*v(2,2) + cyp*v(3,2)) &
      !      +     czp*(cym*v(1,3) + cy0*v(2,3) + cyp*v(3,3))

   end subroutine interp2d

end module eb_wallflux_mod
