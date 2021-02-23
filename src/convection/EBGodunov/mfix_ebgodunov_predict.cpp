#include <mfix_ebgodunov_corner_couple.H>
#include <mfix_godunov_trans_bc.H>
#include <mfix_ebgodunov_transverse_K.H>

#include <Godunov.H>
#include <EBGodunov.H>

using namespace amrex;

void ebgodunov::predict_godunov (Real /*time*/,
                                 MultiFab& u_mac,
                                 MultiFab& v_mac,
                                 MultiFab& w_mac,
                                 MultiFab const& vel,
                                 MultiFab const& vel_forces,
                                 Vector<BCRec> const& h_bcrec,
                                        BCRec  const* d_bcrec,
                                 EBFArrayBoxFactory const* ebfact,
                                 Geometry& geom,
                                 Real l_dt,
                                 MultiFab const& gmacphi_x,
                                 MultiFab const& gmacphi_y,
                                 MultiFab const& gmacphi_z,
                                 bool use_mac_phi_in_godunov)
{
    Box const& domain = geom.Domain();
    const Real* dx    = geom.CellSize();

    auto const& flags = ebfact->getMultiEBCellFlagFab();
    auto const& fcent = ebfact->getFaceCent();
    auto const& ccent = ebfact->getCentroid();
    auto const& vfrac = ebfact->getVolFrac();

    auto const& areafrac = ebfact->getAreaFrac();

    // Since we don't fill all the ghost cells in the mac vel arrays
    // we need to initialize to something which won't make the code crash
    u_mac.setVal(1.e40);
    v_mac.setVal(1.e40);
    w_mac.setVal(1.e40);

    const int ncomp = AMREX_SPACEDIM;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FArrayBox scratch;
        for (MFIter mfi(vel,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Box const& bxg2 = amrex::grow(bx,2);

            EBCellFlagFab const& flagfab = flags[mfi];
            Array4<EBCellFlag const> const& flagarr = flagfab.const_array();

            Array4<Real      > const& a_umac      = u_mac.array(mfi);
            Array4<Real      > const& a_vmac      = v_mac.array(mfi);
            Array4<Real      > const& a_wmac      = w_mac.array(mfi);

            Array4<Real const> const& gmacphi_x_arr = gmacphi_x.const_array(mfi);
            Array4<Real const> const& gmacphi_y_arr = gmacphi_y.const_array(mfi);
            Array4<Real const> const& gmacphi_z_arr = gmacphi_z.const_array(mfi);

            Array4<Real const> const& a_vel       = vel.const_array(mfi);
            Array4<Real const> const& a_f         = vel_forces.const_array(mfi);

            // Not all the arrays have exactly the size of bxg3 but none are bigger
            // In 2-d:
            //  8*ncomp are:  Imx, Ipx, Imy, Ipy, xlo/xhi, ylo/yhi
            //  2       are:  u_ad, v_ad
            // In 3-d:
            // 12*ncomp are:  Imx, Ipx, Imy, Ipy, Imz, Ipz, xlo/xhi, ylo/yhi, zlo/zhi
            //  3       are:  u_ad, v_ad, w_ad
            Box const& bxg3 = amrex::grow(bx,3);
            scratch.resize(bxg3, 4*AMREX_SPACEDIM*ncomp + AMREX_SPACEDIM);
            Real* p  = scratch.dataPtr();

//            Elixir eli = scratch.elixir(); // not needed because of streamSynchronize later

            Box const& xbx = mfi.nodaltilebox(0);
            Box const& ybx = mfi.nodaltilebox(1);
            Box const& zbx = mfi.nodaltilebox(2);

            Box xebx(Box(bx).grow(1).surroundingNodes(0));
            Box yebx(Box(bx).grow(1).surroundingNodes(1));
            Box zebx(Box(bx).grow(1).surroundingNodes(2));

            Box xebx_g2(Box(bx).grow(1).grow(1,1).grow(2,1).surroundingNodes(0));
            Box yebx_g2(Box(bx).grow(1).grow(0,1).grow(2,1).surroundingNodes(1));
            Box zebx_g2(Box(bx).grow(1).grow(0,1).grow(1,1).surroundingNodes(2));

            Array4<Real> Imx = makeArray4(p,bxg2,ncomp);
            p +=         Imx.size();
            Array4<Real> Ipx = makeArray4(p,bxg2,ncomp);
            p +=         Ipx.size();
            Array4<Real> Imy = makeArray4(p,bxg2,ncomp);
            p +=         Imy.size();
            Array4<Real> Ipy = makeArray4(p,bxg2,ncomp);
            p +=         Ipy.size();

            Array4<Real> u_ad = makeArray4(p,xebx_g2,1);
            p +=         u_ad.size();
            Array4<Real> v_ad = makeArray4(p,yebx_g2,1);
            p +=         v_ad.size();

            Array4<Real> Imz = makeArray4(p,bxg2,ncomp);
            p +=         Imz.size();
            Array4<Real> Ipz = makeArray4(p,bxg2,ncomp);
            p +=         Ipz.size();

            Array4<Real> w_ad = makeArray4(p,zebx_g2,1);
            p +=         w_ad.size();

            // This tests on covered cells just in the box itself
            if (flagfab.getType(bx) == FabType::covered)
            {

            // This tests on only regular cells including two rows of ghost cells
            } else if (flagfab.getType(amrex::grow(bx,2)) == FabType::regular) {

                godunov::predict_plm_x (xebx_g2, Imx, Ipx, a_vel, a_vel,
                                        geom, l_dt, h_bcrec, d_bcrec);
                godunov::predict_plm_y (yebx_g2, Imy, Ipy, a_vel, a_vel,
                                        geom, l_dt, h_bcrec, d_bcrec);
                godunov::predict_plm_z (zebx_g2, Imz, Ipz, a_vel, a_vel,
                                        geom, l_dt, h_bcrec, d_bcrec);

                bool use_forces_in_trans = false;
                godunov::make_trans_velocities(AMREX_D_DECL(Box(u_ad), Box(v_ad), Box(w_ad)),
                                               AMREX_D_DECL(u_ad, v_ad, w_ad),
                                               AMREX_D_DECL(Imx, Imy, Imz),
                                               AMREX_D_DECL(Ipx, Ipy, Ipz), a_vel, a_f,
                                               domain, l_dt, d_bcrec, use_forces_in_trans);

                godunov::predict_godunov_on_box(bx, ncomp,
                                       AMREX_D_DECL(xbx, ybx, zbx),
                                       AMREX_D_DECL(a_umac, a_vmac, a_wmac),
                                       a_vel, AMREX_D_DECL(u_ad, v_ad, w_ad),
                                       AMREX_D_DECL(Imx, Imy, Imz),
                                       AMREX_D_DECL(Ipx, Ipy, Ipz), a_f,
                                       domain, dx, l_dt, d_bcrec, use_forces_in_trans,
                                       AMREX_D_DECL(gmacphi_x_arr, gmacphi_y_arr, gmacphi_z_arr),
                                       use_mac_phi_in_godunov, p);

            } else {

                Array4<Real const> const& fcx = fcent[0]->const_array(mfi);
                Array4<Real const> const& fcy = fcent[1]->const_array(mfi);
                Array4<Real const> const& fcz = fcent[2]->const_array(mfi);

                Array4<Real const> const& ccent_arr = ccent.const_array(mfi);
                Array4<Real const> const& vfrac_arr = vfrac.const_array(mfi);

                ebgodunov::predict_plm_x(xebx_g2, Imx, Ipx, a_vel, a_vel,
                                         flagarr, vfrac_arr,
                                         AMREX_D_DECL(fcx,fcy,fcz),ccent_arr,
                                         geom, l_dt, h_bcrec, d_bcrec);
                ebgodunov::predict_plm_y(yebx_g2, Imy, Ipy, a_vel, a_vel,
                                         flagarr, vfrac_arr,
                                         AMREX_D_DECL(fcx,fcy,fcz),ccent_arr,
                                         geom, l_dt, h_bcrec, d_bcrec);
                ebgodunov::predict_plm_z(zebx_g2, Imz, Ipz, a_vel, a_vel,
                                         flagarr, vfrac_arr,
                                         AMREX_D_DECL(fcx,fcy,fcz),ccent_arr,
                                         geom, l_dt, h_bcrec, d_bcrec);

                ebgodunov::make_trans_velocities(AMREX_D_DECL(xebx_g2, yebx_g2, zebx_g2),
                                                 AMREX_D_DECL(u_ad, v_ad, w_ad),
                                                 AMREX_D_DECL(Imx, Imy, Imz),
                                                 AMREX_D_DECL(Ipx, Ipy, Ipz), a_vel,
                                                 flagarr, domain, d_bcrec);

                Array4<Real const> const& apx = areafrac[0]->const_array(mfi);
                Array4<Real const> const& apy = areafrac[1]->const_array(mfi);
                Array4<Real const> const& apz = areafrac[2]->const_array(mfi);

                ebgodunov::predict_godunov_on_box(bx, ncomp,
                                       AMREX_D_DECL(xbx,ybx,zbx),
                                       AMREX_D_DECL(xebx_g2,yebx_g2,zebx_g2),
                                       AMREX_D_DECL(a_umac, a_vmac, a_wmac),
                                       a_vel,
                                       AMREX_D_DECL(u_ad, v_ad, w_ad),
                                       AMREX_D_DECL(Imx, Imy, Imz),
                                       AMREX_D_DECL(Ipx, Ipy, Ipz),
                                       a_f,
                                       domain, dx, l_dt, d_bcrec,
                                       flagarr,
                                       AMREX_D_DECL(apx, apy, apz),
                                       vfrac_arr,
                                       AMREX_D_DECL(fcx, fcy, fcz),
                                       AMREX_D_DECL(gmacphi_x_arr, gmacphi_y_arr, gmacphi_z_arr),
                                       use_mac_phi_in_godunov, p);
            }

            Gpu::streamSynchronize();  // otherwise we might be using too much memory
        }
    }
}

void ebgodunov::predict_godunov_on_box (Box const& bx, int ncomp,
                                        Box const& xbx, Box const& ybx, Box const& zbx,
                                        Box const& xebx, Box const& yebx, Box const& zebx,
                                        Array4<Real> const& qx,
                                        Array4<Real> const& qy,
                                        Array4<Real> const& qz,
                                        Array4<Real const> const& q,
                                        Array4<Real const> const& u_ad,
                                        Array4<Real const> const& v_ad,
                                        Array4<Real const> const& w_ad,
                                        Array4<Real> const& Imx,
                                        Array4<Real> const& Imy,
                                        Array4<Real> const& Imz,
                                        Array4<Real> const& Ipx,
                                        Array4<Real> const& Ipy,
                                        Array4<Real> const& Ipz,
                                        Array4<Real const> const& f,
                                        const Box& domain,
                                        const Real* dx_arr,
                                        Real l_dt,
                                        BCRec  const* pbc,
                                        Array4<EBCellFlag const> const& flag,
                                        Array4<Real const> const& apx,
                                        Array4<Real const> const& apy,
                                        Array4<Real const> const& apz,
                                        Array4<Real const> const& vfrac_arr,
                                        Array4<Real const> const& fcx,
                                        Array4<Real const> const& fcy,
                                        Array4<Real const> const& fcz,
                                        Array4<Real const> const& gmacphi_x,
                                        Array4<Real const> const& gmacphi_y,
                                        Array4<Real const> const& gmacphi_z,
                                        bool l_use_mac_phi_in_godunov,
                                        Real* p)
{
    const Dim3 dlo = amrex::lbound(domain);
    const Dim3 dhi = amrex::ubound(domain);
    Real dx = dx_arr[0];
    Real dy = dx_arr[1];
    Real dz = dx_arr[2];

    Box xebx_g2 = Box(bx).grow(2).surroundingNodes(0);
    Box yebx_g2 = Box(bx).grow(2).surroundingNodes(1);
    Box zebx_g2 = Box(bx).grow(2).surroundingNodes(2);

    Array4<Real> xlo = makeArray4(p, Box(xebx_g2), ncomp);
    p += xlo.size();
    Array4<Real> xhi = makeArray4(p, Box(xebx_g2), ncomp);
    p += xhi.size();
    Array4<Real> ylo = makeArray4(p, Box(yebx_g2), ncomp);
    p += ylo.size();
    Array4<Real> yhi = makeArray4(p, Box(yebx_g2), ncomp);
    p += yhi.size();
    Array4<Real> zlo = makeArray4(p, Box(zebx_g2), ncomp);
    p += zlo.size();
    Array4<Real> zhi = makeArray4(p, Box(zebx_g2), ncomp);
    p += zhi.size();

    amrex::ParallelFor(
        xebx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real lo = Ipx(i-1,j,k,n);
            Real hi = Imx(i  ,j,k,n);

            Real uad = u_ad(i,j,k);
            auto bc = pbc[n];

            Godunov_trans_xbc(i, j, k, n, q, lo, hi, bc.lo(0), bc.hi(0), dlo.x, dhi.x, true);

            xlo(i,j,k,n) = lo;
            xhi(i,j,k,n) = hi;

            Real st = (uad >= 0.) ? lo : hi;
            Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
            Imx(i, j, k, n) = fu*st + (1.0 - fu) *0.5 * (hi + lo); // store xedge
        },
        yebx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real lo = Ipy(i,j-1,k,n);
            Real hi = Imy(i,j  ,k,n);

            Real vad = v_ad(i,j,k);
            auto bc = pbc[n];

            Godunov_trans_ybc(i, j, k, n, q, lo, hi, bc.lo(1), bc.hi(1), dlo.y, dhi.y, true);

            ylo(i,j,k,n) = lo;
            yhi(i,j,k,n) = hi;

            Real st = (vad >= 0.) ? lo : hi;
            Real fu = (amrex::Math::abs(vad) < small_vel) ? 0.0 : 1.0;
            Imy(i, j, k, n) = fu*st + (1.0 - fu)*0.5*(hi + lo); // store yedge
        },
        zebx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real lo = Ipz(i,j,k-1,n);
            Real hi = Imz(i,j,k  ,n);

            Real wad = w_ad(i,j,k);
            auto bc = pbc[n];

            Godunov_trans_zbc(i, j, k, n, q, lo, hi, bc.lo(2), bc.hi(2), dlo.z, dhi.z, true);

            zlo(i,j,k,n) = lo;
            zhi(i,j,k,n) = hi;

            Real st = (wad >= 0.) ? lo : hi;
            Real fu = (amrex::Math::abs(wad) < small_vel) ? 0.0 : 1.0;
            Imz(i, j, k, n) = fu*st + (1.0 - fu)*0.5*(hi + lo); // store zedge
        });

    Array4<Real> divu = makeArray4(Ipx.dataPtr(), grow(bx,1), 1);
    amrex::ParallelFor(Box(divu), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        divu(i,j,k) = 0.0;
    });

    // We can reuse the space in Ipy and Ipz.
    Array4<Real> xedge = Imx;
    Array4<Real> yedge = Imy;
    Array4<Real> zedge = Imz;

    //
    // X-Flux
    //

    // Grow in x-direction
    Box const xbxtmp = Box(xbx).enclosedCells().grow(0,1);

    // Need to grow these in y,z directions 1) because these are on y/z faces
    //                                      2) because we will do tangential interpolation
    Array4<Real> yzlo = makeArray4(Ipy.dataPtr(), amrex::surroundingNodes(xbxtmp,1).grow(2,1).grow(0,1), 1);
    Array4<Real> zylo = makeArray4(Ipz.dataPtr(), amrex::surroundingNodes(xbxtmp,2).grow(1,1).grow(0,1), 1);

    // Add d/dy term to z-faces
    // Start with {zlo,zhi} --> {zylo, zyhi} and upwind using w_ad to {zylo}
    // Add d/dz to y-faces
    // Start with {ylo,yhi} --> {yzlo, yzhi} and upwind using v_ad to {yzlo}

    amrex::ParallelFor(Box(zylo), Box(yzlo),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 0;
        auto bc = pbc[n];
        Real l_zylo, l_zyhi;
        EBGodunov_corner_couple_zy(l_zylo, l_zyhi,
                                   i, j, k, n, l_dt, dy, false,
                                   zlo(i,j,k,n), zhi(i,j,k,n),
                                   q, divu, apx, apy, apz, vfrac_arr, v_ad, yedge);

        Real wad = w_ad(i,j,k);
        Godunov_trans_zbc(i, j, k, n, q, l_zylo, l_zyhi, bc.lo(2), bc.hi(2), dlo.z, dhi.z, true);


        Real st = (wad >= 0.) ? l_zylo : l_zyhi;
        Real fu = (amrex::Math::abs(wad) < small_vel) ? 0.0 : 1.0;
        zylo(i,j,k) = fu*st + (1.0 - fu) * 0.5 * (l_zyhi + l_zylo);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 0;
        auto bc = pbc[n];
        Real l_yzlo, l_yzhi;
        EBGodunov_corner_couple_yz(l_yzlo, l_yzhi,
                                   i, j, k, n, l_dt, dz, false,
                                   ylo(i,j,k,n), yhi(i,j,k,n),
                                   q, divu, apx, apy, apz, vfrac_arr, w_ad, zedge);

        Real vad = v_ad(i,j,k);
        Godunov_trans_ybc(i, j, k, n, q, l_yzlo, l_yzhi, bc.lo(1), bc.hi(1), dlo.y, dhi.y, true);

        Real st = (vad >= 0.) ? l_yzlo : l_yzhi;
        Real fu = (amrex::Math::abs(vad) < small_vel) ? 0.0 : 1.0;
        yzlo(i,j,k) = fu*st + (1.0 - fu) * 0.5 * (l_yzhi + l_yzlo);
    });
    //
    amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (flag(i,j,k).isConnected(-1,0,0))
        {
        constexpr int n = 0;
        auto bc = pbc[n];

        // stl is on the lo side of the lo-x side of cell (i,j,k)
        // sth is on the hi side of the lo-x side of cell (i,j,k)
        Real stl = xlo(i,j,k,n);
        Real sth = xhi(i,j,k,n);

        Real trans_y, trans_z;

        //
        // Left side of interface
        //
        {
        int ic = i-1;
        if (flag(ic,j,k).isRegular())
        {
            stl += - (0.25*l_dt/dy)*(v_ad(ic,j+1,k  )+v_ad(ic,j,k))*
                                    (yzlo(ic,j+1,k  )-yzlo(ic,j,k))
                   - (0.25*l_dt/dz)*(w_ad(ic,j  ,k+1)+w_ad(ic,j,k))*
                                    (zylo(ic,j  ,k+1)-zylo(ic,j,k))
                   + 0.5 * l_dt * f(ic,j,k,n);

        // Only add dt-based terms if we can construct all transverse terms
        //    using non-covered faces
        } else if (apy(ic,j+1,k) > 0. && apy(ic,j,k) > 0. &&
                   apz(ic,j,k+1) > 0. && apz(ic,j,k) > 0.)
        {
            create_transverse_terms_for_xface(ic, j, k, v_ad, w_ad, yzlo, zylo,
                                              apy, apz, fcy, fcz, trans_y, trans_z,
                                              dy, dz);

            stl += -0.5 * l_dt * (trans_y + trans_z)
                   +0.5 * l_dt * f(ic,j,k,n);
        }
        }

        //
        // Right side of interface
        //
        {
        int ic = i;
        if (flag(ic,j,k).isRegular())
        {
             sth += - (0.25*l_dt/dy)*(v_ad(ic,j+1,k  )+v_ad(ic,j,k))*
                                     (yzlo(ic,j+1,k  )-yzlo(ic,j,k))
                    - (0.25*l_dt/dz)*(w_ad(ic,j  ,k+1)+w_ad(ic,j,k))*
                                     (zylo(ic,j  ,k+1)-zylo(ic,j,k))
                    + 0.5 * l_dt * f(ic,j,k,n);

        // Only add dt-based terms if we can construct all transverse terms
        //    using non-covered faces
        } else if (apy(ic,j+1,k) > 0. && apy(ic,j,k) > 0. &&
                   apz(ic,j,k+1) > 0. && apz(ic,j,k) > 0.)
        {
            create_transverse_terms_for_xface(ic, j, k, v_ad, w_ad, yzlo, zylo,
                                              apy, apz, fcy, fcz, trans_y, trans_z,
                                              dy, dz);

            sth += -0.5 * l_dt * (trans_y + trans_z)
                   +0.5 * l_dt * f(ic,j,k,n);
        }
        }

        Real gphi_x;

        if (l_use_mac_phi_in_godunov)
        {
            // Note that the getFluxes call returns (-1/rho G^Mac phi) so we use the negative here
            gphi_x = -gmacphi_x(i,j,k);
            stl -= 0.5 * l_dt * gphi_x;
            sth -= 0.5 * l_dt * gphi_x;
        }

        Godunov_cc_xbc_lo(i, j, k, n, q, stl, sth, bc.lo(0), dlo.x, true);
        Godunov_cc_xbc_hi(i, j, k, n, q, stl, sth, bc.hi(0), dhi.x, true);

        // Prevent backflow
        if ( (i==dlo.x) and (bc.lo(0) == BCType::foextrap || bc.lo(0) == BCType::hoextrap) )
        {
            sth = amrex::min(sth,0.);
            stl = sth;
        }
        if ( (i==dhi.x+1) and (bc.hi(0) == BCType::foextrap || bc.hi(0) == BCType::hoextrap) )
        {
             stl = amrex::max(stl,0.);
             sth = stl;
        }

        Real st = ( (stl+sth) >= 0.) ? stl : sth;
        bool ltm = ( (stl <= 0. && sth >= 0.) || (amrex::Math::abs(stl+sth) < small_vel) );
        qx(i,j,k) = ltm ? 0. : st;

        if (l_use_mac_phi_in_godunov)
            qx(i,j,k) += 0.5 * l_dt * gphi_x;

        } else {
            qx(i,j,k) = 0.;
        }
    });

    //
    // Y-Flux
    //
    // Grow in y-direction
    Box const ybxtmp = Box(ybx).enclosedCells().grow(1,1);

    // Need to grow these in x,z directions 1) because these are on x/z faces
    //                                      2) because we will do tangential interpolation
    Array4<Real> xzlo = makeArray4(Ipy.dataPtr(), amrex::surroundingNodes(ybxtmp,0).grow(2,1).grow(1,1), 1);
    Array4<Real> zxlo = makeArray4(Ipz.dataPtr(), amrex::surroundingNodes(ybxtmp,2).grow(0,1).grow(1,1), 1);

    // Add d/dz to x-faces
    // Start with {xlo,xhi} --> {xzlo, xzhi} and upwind using u_ad to {xzlo}
    // Add d/dx term to z-faces
    // Start with {zlo,zhi} --> {zxlo, zxhi} and upwind using w_ad to {zxlo}
    amrex::ParallelFor(Box(xzlo), Box(zxlo),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 1;
        const auto bc = pbc[n];
        Real l_xzlo, l_xzhi;
        EBGodunov_corner_couple_xz(l_xzlo, l_xzhi,
                                   i, j, k, n, l_dt, dz, false,
                                   xlo(i,j,k,n),  xhi(i,j,k,n),
                                   q, divu, apx, apy, apz, vfrac_arr, w_ad, zedge);

        Real uad = u_ad(i,j,k);
        Godunov_trans_xbc(i, j, k, n, q, l_xzlo, l_xzhi, bc.lo(0), bc.hi(0), dlo.x, dhi.x, true);


        Real st = (uad >= 0.) ? l_xzlo : l_xzhi;
        Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
        xzlo(i,j,k) = fu*st + (1.0 - fu) * 0.5 * (l_xzhi + l_xzlo);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 1;
        const auto bc = pbc[n];
        Real l_zxlo, l_zxhi;
        EBGodunov_corner_couple_zx(l_zxlo, l_zxhi,
                                   i, j, k, n, l_dt, dx, false,
                                   zlo(i,j,k,n), zhi(i,j,k,n),
                                   q, divu, apx, apy, apz, vfrac_arr, u_ad, xedge);

        Real wad = w_ad(i,j,k);
        Godunov_trans_zbc(i, j, k, n, q, l_zxlo, l_zxhi, bc.lo(2), bc.hi(2), dlo.z, dhi.z, true);


        Real st = (wad >= 0.) ? l_zxlo : l_zxhi;
        Real fu = (amrex::Math::abs(wad) < small_vel) ? 0.0 : 1.0;
        zxlo(i,j,k) = fu*st + (1.0 - fu) * 0.5 * (l_zxhi + l_zxlo);
    });
    //
    amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (flag(i,j,k).isConnected(0,-1,0))
        {
        constexpr int n = 1;
        auto bc = pbc[n];

        // stl is on the lo side of the lo-y side of cell (i,j,k)
        // sth is on the hi side of the lo-y side of cell (i,j,k)
        Real stl = ylo(i,j,k,n);
        Real sth = yhi(i,j,k,n);

        Real trans_x, trans_z;

        //
        // Left side of interface
        //
        {
        int jc = j-1;
        if (flag(i,jc,k).isRegular())
        {
            stl += - (0.25*l_dt/dx)*(u_ad(i+1,jc,k  )+u_ad(i,jc,k))*
                                    (xzlo(i+1,jc,k  )-xzlo(i,jc,k))
                   - (0.25*l_dt/dz)*(w_ad(i  ,jc,k+1)+w_ad(i,jc,k))*
                                    (zxlo(i  ,jc,k+1)-zxlo(i,jc,k))
                   + 0.5 * l_dt * f(i,jc,k,n);

        // Only add dt-based terms if we can construct all transverse terms
        //    using non-covered faces
        } else if (apx(i+1,jc,k  ) > 0. && apx(i,jc,k) > 0. &&
                   apz(i  ,jc,k+1) > 0. && apz(i,jc,k) > 0.)
        {
            create_transverse_terms_for_yface(i, jc, k, u_ad, w_ad, xzlo, zxlo,
                                              apx, apz, fcx, fcz, trans_x, trans_z,
                                              dx, dz);

            stl += -0.5 * l_dt * (trans_x + trans_z)
                   +0.5 * l_dt * f(i,jc,k,n);
        }
        }

        //
        // Right side of interface
        //
        {
        int jc = j;
        if (flag(i,jc,k).isRegular())
        {
            sth += - (0.25*l_dt/dx)*(u_ad(i+1,jc,k  )+u_ad(i,jc,k))*
                                    (xzlo(i+1,jc,k  )-xzlo(i,jc,k))
                   - (0.25*l_dt/dz)*(w_ad(i  ,jc,k+1)+w_ad(i,jc,k))*
                                    (zxlo(i  ,jc,k+1)-zxlo(i,jc,k))
                   +0.5 * l_dt * f(i,jc,k,n);

        // Only add dt-based terms if we can construct all transverse terms
        //    using non-covered faces
        } else if (apx(i+1,jc,k  ) > 0. && apx(i,jc,k) > 0. &&
                   apz(i  ,jc,k+1) > 0. && apz(i,jc,k) > 0.)
        {
            create_transverse_terms_for_yface(i, jc, k, u_ad, w_ad, xzlo, zxlo,
                                              apx, apz, fcx, fcz, trans_x, trans_z,
                                              dx, dz);

            sth += -0.5 * l_dt * (trans_x + trans_z)
                   +0.5 * l_dt * f(i,jc,k,n);
        }
        }

        Real gphi_y;

        if (l_use_mac_phi_in_godunov)
        {
            // Note that the getFluxes call returns (-1/rho G^Mac phi) so we use the negative here
            gphi_y = -gmacphi_y(i,j,k);
            stl -= 0.5 * l_dt * gphi_y;
            sth -= 0.5 * l_dt * gphi_y;
        }

        Godunov_cc_ybc_lo(i, j, k, n, q, stl, sth, bc.lo(1), dlo.y, true);
        Godunov_cc_ybc_hi(i, j, k, n, q, stl, sth, bc.hi(1), dhi.y, true);

        // Prevent backflow
        if ( (j==dlo.y) and (bc.lo(1) == BCType::foextrap || bc.lo(1) == BCType::hoextrap) )
        {
            sth = amrex::min(sth,0.);
            stl = sth;
        }
        if ( (j==dhi.y+1) and (bc.hi(1) == BCType::foextrap || bc.hi(1) == BCType::hoextrap) )
        {
            stl = amrex::max(stl,0.);
            sth = stl;
        }

        Real st = ( (stl+sth) >= 0.) ? stl : sth;
        bool ltm = ( (stl <= 0. && sth >= 0.) || (amrex::Math::abs(stl+sth) < small_vel) );
        qy(i,j,k) = ltm ? 0. : st;

        if (l_use_mac_phi_in_godunov)
            qy(i,j,k) += 0.5 * l_dt * gphi_y;

        } else {
            qy(i,j,k) = 0.;
        }
    });

    //
    // Z-Flux
    //
    // Grow in z-direction
    Box const zbxtmp = Box(zbx).enclosedCells().grow(2,1);

    // Need to grow these in x,y directions 1) because these are on x/y faces
    //                                      2) because we will do tangential interpolation
    Array4<Real> xylo = makeArray4(Ipy.dataPtr(), amrex::surroundingNodes(zbxtmp,0).grow(1,1).grow(2,1), 1);
    Array4<Real> yxlo = makeArray4(Ipz.dataPtr(), amrex::surroundingNodes(zbxtmp,1).grow(0,1).grow(2,1), 1);

    amrex::ParallelFor(Box(xylo), Box(yxlo),
    //
    // Add d/dy term to x-faces
    // Start with {xlo,xhi} --> {xylo, xyhi} and upwind using u_ad to {xylo}
    //
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 2;
        const auto bc = pbc[n];
        Real l_xylo, l_xyhi;
        EBGodunov_corner_couple_xy(l_xylo, l_xyhi,
                                   i, j, k, n, l_dt, dy, false,
                                   xlo(i,j,k,n), xhi(i,j,k,n),
                                   q, divu, apx, apy, apz, vfrac_arr, v_ad, yedge);

        Real uad = u_ad(i,j,k);
        Godunov_trans_xbc(i, j, k, n, q, l_xylo, l_xyhi, bc.lo(0), bc.hi(0), dlo.x, dhi.x, true);


        Real st = (uad >= 0.) ? l_xylo : l_xyhi;
        Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
        xylo(i,j,k) = fu*st + (1.0 - fu) * 0.5 * (l_xyhi + l_xylo);
    },
    //
    // Add d/dx term to y-faces
    // Start with {ylo,yhi} --> {yxlo, yxhi} and upwind using v_ad to {yxlo}
    //
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 2;
        const auto bc = pbc[n];
        Real l_yxlo, l_yxhi;
        EBGodunov_corner_couple_yx(l_yxlo, l_yxhi,
                                   i, j, k, n, l_dt, dx, false,
                                   ylo(i,j,k,n), yhi(i,j,k,n),
                                   q, divu, apx, apy, apz, vfrac_arr, u_ad, xedge);

        Real vad = v_ad(i,j,k);
        Godunov_trans_ybc(i, j, k, n, q, l_yxlo, l_yxhi, bc.lo(1), bc.hi(1), dlo.y, dhi.y, true);


        Real st = (vad >= 0.) ? l_yxlo : l_yxhi;
        Real fu = (amrex::Math::abs(vad) < small_vel) ? 0.0 : 1.0;
        yxlo(i,j,k) = fu*st + (1.0 - fu) * 0.5 * (l_yxhi + l_yxlo);
    });
    //
    amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (flag(i,j,k).isConnected(0,0,-1))
        {
        constexpr int n = 2;
        auto bc = pbc[n];

        // stl is on the lo side of the lo-z side of cell (i,j,k)
        // sth is on the hi side of the lo-z side of cell (i,j,k)
        Real stl = zlo(i,j,k,n);
        Real sth = zhi(i,j,k,n);

        Real trans_x, trans_y;

        //
        // Lo side of interface
        //
        {
        int kc = k-1;
        if (flag(i,j,kc).isRegular())
        {
            stl += - (0.25*l_dt/dx)*(u_ad(i+1,j  ,kc)+u_ad(i,j,kc))*
                                    (xylo(i+1,j  ,kc)-xylo(i,j,kc))
                   - (0.25*l_dt/dy)*(v_ad(i  ,j+1,kc)+v_ad(i,j,kc))*
                                    (yxlo(i  ,j+1,kc)-yxlo(i,j,kc))
                   + 0.5 * l_dt * f(i,j,kc,n);

        // Only add dt-based terms if we can construct all transverse terms
        //    using non-covered faces
        } else if (apx(i+1,j  ,kc) > 0. && apx(i,j,kc) > 0. &&
                   apy(i  ,j+1,kc) > 0. && apy(i,j,kc) > 0.)
        {
            create_transverse_terms_for_zface(i, j, kc, u_ad, v_ad, xylo, yxlo,
                                              apx, apy, fcx, fcy, trans_x, trans_y,
                                              dx, dy);

            stl += -0.5 * l_dt * (trans_x + trans_y)
                   +0.5 * l_dt * f(i,j,kc,n);
        }
        }

        //
        // Right side of interface
        //
        {
        int kc = k;
        if (flag(i,j,kc).isRegular())
        {
            sth += - (0.25*l_dt/dx)*(u_ad(i+1,j  ,kc)+u_ad(i,j,kc))*
                                    (xylo(i+1,j  ,kc)-xylo(i,j,kc))
                   - (0.25*l_dt/dy)*(v_ad(i  ,j+1,kc)+v_ad(i,j,kc))*
                                    (yxlo(i  ,j+1,kc)-yxlo(i,j,kc))
                   + 0.5 * l_dt * f(i,j,kc,n);

        // Only add dt-based terms if we can construct all transverse terms
        //    using non-covered faces
        } else if (apx(i+1,j  ,kc) > 0. && apx(i,j,kc) > 0. &&
                   apy(i  ,j+1,kc) > 0. && apy(i,j,kc) > 0.)
        {
            create_transverse_terms_for_zface(i, j, kc, u_ad, v_ad, xylo, yxlo,
                                              apx, apy, fcx, fcy, trans_x, trans_y,
                                              dx, dy);

            sth += -0.5 * l_dt * (trans_x + trans_y)
                   +0.5 * l_dt * f(i,j,kc,n);
        }
        }

        Real gphi_z;

        if (l_use_mac_phi_in_godunov)
        {
            // Note that the getFluxes call returns (-1/rho G^Mac phi) so we use the negative here
            gphi_z = -gmacphi_z(i,j,k);
            stl -= 0.5 * l_dt * gphi_z;
            sth -= 0.5 * l_dt * gphi_z;
        }

        Godunov_cc_zbc_lo(i, j, k, n, q, stl, sth, bc.lo(2), dlo.z, true);
        Godunov_cc_zbc_hi(i, j, k, n, q, stl, sth, bc.hi(2), dhi.z, true);

        // Prevent backflow
        if ( (k==dlo.z) and (bc.lo(2) == BCType::foextrap || bc.lo(2) == BCType::hoextrap) )
        {
            sth = amrex::min(sth,0.);
            stl = sth;
        }
        if ( (k==dhi.z+1) and (bc.hi(2) == BCType::foextrap || bc.hi(2) == BCType::hoextrap) )
        {
            stl = amrex::max(stl,0.);
            sth = stl;
        }

        Real st = ( (stl+sth) >= 0.) ? stl : sth;
        bool ltm = ( (stl <= 0. && sth >= 0.) || (amrex::Math::abs(stl+sth) < small_vel) );
        qz(i,j,k) = ltm ? 0. : st;

        if (l_use_mac_phi_in_godunov)
            qz(i,j,k) += 0.5 * l_dt * gphi_z;

        } else {
            qz(i,j,k) = 0.;
        }
    });
}
