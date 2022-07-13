#include <mfix.H>

#include <mfix_fluid.H>
#include <mfix_bc.H>
#include <mfix_ic.H>

#include <AMReX_GpuDevice.H>

#include <string>
#include <climits>

namespace set_p0_aux {

void compute_p0_bcs (const Box& sbx,
                     const Box& domain,
                     Array4<Real> const& p0_g,
                     Gpu::DeviceVector<Real>& m_bc_p_g, Real& pj,
                     const RealVect& gravity,
                     const Real dx,
                     const Real dy,
                     const Real dz,
                     Array4<const int> const& bct_ilo,
                     Array4<const int> const& bct_ihi,
                     Array4<const int> const& bct_jlo,
                     Array4<const int> const& bct_jhi,
                     Array4<const int> const& bct_klo,
                     Array4<const int> const& bct_khi,
                     const int nlft,
                     const int nrgt,
                     const int nbot,
                     const int ntop,
                     const int ndwn,
                     const int nup,
                     const int nghost,
                     const Real ro_g0,
                     const MFIXBoundaryConditions& boundary_conditions);

void set_p0_bcs (const Box& sbx,
                 const Box& domain,
                 Array4<Real> const& p0_g,
                 Gpu::DeviceVector<Real>& m_bc_p_g,
                 Array4<const int> const& bct_ilo,
                 Array4<const int> const& bct_ihi,
                 Array4<const int> const& bct_jlo,
                 Array4<const int> const& bct_jhi,
                 Array4<const int> const& bct_klo,
                 Array4<const int> const& bct_khi,
                 const int nlft,
                 const int nrgt,
                 const int nbot,
                 const int ntop,
                 const int ndwn,
                 const int nup,
                 const int nghost);

} // end namespace set_p0_aux


using namespace set_p0_aux;


void
mfix::set_p0 (const Box& bx,
              MFIter* mfi,
              const int lev,
              const Box& domain)
{
  Real dx = geom[lev].CellSize(0);
  Real dy = geom[lev].CellSize(1);
  Real dz = geom[lev].CellSize(2);

  Real xlen = geom[0].ProbHi(0) - geom[0].ProbLo(0);
  Real ylen = geom[0].ProbHi(1) - geom[0].ProbLo(1);
  Real zlen = geom[0].ProbHi(2) - geom[0].ProbLo(2);

  const Real tolerance = std::numeric_limits<Real>::epsilon();
  Real offset(-0.5);

  Array4<double> const& array4_p0_g = m_leveldata[lev]->p0_g->array(*mfi);

  const IntVect sbx_lo((*m_leveldata[lev]->p0_g)[*mfi].loVect());
  const IntVect sbx_hi((*m_leveldata[lev]->p0_g)[*mfi].hiVect());
  const Box sbx(sbx_lo, sbx_hi);

  const IntVect dom_lo(domain.loVect());
  const IntVect dom_hi(domain.hiVect());

  const IntVect bx_lo(bx.loVect());
  const IntVect bx_hi(bx.hiVect());

  auto make_ifab_host = [] (IArrayBox const* dfab) -> IArrayBox
  {
      IArrayBox ifab_host(dfab->box(),dfab->nComp(),The_Pinned_Arena());
#ifdef AMREX_USE_GPU
      Gpu::dtoh_memcpy_async(ifab_host.dataPtr(), dfab->dataPtr(), ifab_host.nBytes());
#else
      std::memcpy(ifab_host.dataPtr(), dfab->dataPtr(), ifab_host.nBytes());
#endif
      return ifab_host;
  };

  IArrayBox const& bc_ilo_host = make_ifab_host(bc_list.bc_ilo[lev]);
  IArrayBox const& bc_ihi_host = make_ifab_host(bc_list.bc_ihi[lev]);
  IArrayBox const& bc_jlo_host = make_ifab_host(bc_list.bc_jlo[lev]);
  IArrayBox const& bc_jhi_host = make_ifab_host(bc_list.bc_jhi[lev]);
  IArrayBox const& bc_klo_host = make_ifab_host(bc_list.bc_klo[lev]);
  IArrayBox const& bc_khi_host = make_ifab_host(bc_list.bc_khi[lev]);
  Gpu::synchronize();

  Array4<const int> const& bct_ilo = bc_ilo_host.const_array();
  Array4<const int> const& bct_ihi = bc_ihi_host.const_array();
  Array4<const int> const& bct_jlo = bc_jlo_host.const_array();
  Array4<const int> const& bct_jhi = bc_jhi_host.const_array();
  Array4<const int> const& bct_klo = bc_klo_host.const_array();
  Array4<const int> const& bct_khi = bc_khi_host.const_array();

  const int nlft = amrex::max(0, dom_lo[0]-sbx_lo[0]+1);
  const int nbot = amrex::max(0, dom_lo[1]-sbx_lo[1]+1);
  const int ndwn = amrex::max(0, dom_lo[2]-sbx_lo[2]+1);

  const int nrgt = amrex::max(0, sbx_hi[0]-dom_hi[0]);
  const int ntop = amrex::max(0, sbx_hi[1]-dom_hi[1]);
  const int nup  = amrex::max(0, sbx_hi[2]-dom_hi[2]);

  int delp_dir_loc(m_boundary_conditions.delp_dir());

  Real pj(0);

  Real delp_x = m_boundary_conditions.delp(0);
  Real delp_y = m_boundary_conditions.delp(1);
  Real delp_z = m_boundary_conditions.delp(2);

  // ---------------------------------------------------------------->>>
  //     If the bc's are pressure inflow/outflow then be sure to capture that in p0andgp0
  // ---------------------------------------------------------------->>>

  if(((bct_ilo(dom_lo[0]-1,dom_lo[1],dom_lo[2],0) == BCList::pinf) &&
     (bct_ihi(dom_hi[0]+1,dom_lo[1],dom_lo[2],0) == BCList::pout))
    ||
     ((bct_ihi(dom_hi[0]+1,dom_lo[1],dom_lo[2],0) == BCList::pinf) &&
     (bct_ilo(dom_lo[0]-1,dom_lo[1],dom_lo[2],0) == BCList::pout)))
  {
    delp_dir_loc = 0;

    const int bcv_lo = bct_ilo(dom_lo[0]-1,dom_lo[1],dom_lo[2],1);
    const Real p_lo  = m_h_bc_p_g[bcv_lo];

    const int bcv_hi = bct_ihi(dom_hi[0]+1,dom_lo[1],dom_lo[2],1);
    const Real p_hi  = m_h_bc_p_g[bcv_hi];

    delp_x = p_lo - p_hi;
    m_boundary_conditions.set_delp(0, delp_x);

    pj = p_hi;
  }
  else if(((bct_jlo(dom_lo[0],dom_lo[1]-1,dom_lo[2],0) == BCList::pinf) &&
          (bct_jhi(dom_lo[0],dom_hi[1]+1,dom_lo[2],0) == BCList::pout))
         ||
          ((bct_jhi(dom_lo[0],dom_hi[1]+1,dom_lo[2],0) == BCList::pinf) &&
          (bct_jlo(dom_lo[0],dom_lo[1]-1,dom_lo[2],0) == BCList::pout)))
  {
    delp_dir_loc = 1;

    const int bcv_lo = bct_jlo(dom_lo[0],dom_lo[1]-1,dom_lo[2],1);
    const Real p_lo  = m_h_bc_p_g[bcv_lo];

    const int bcv_hi = bct_jhi(dom_lo[0],dom_hi[1]+1,dom_lo[2],1);
    const Real p_hi  = m_h_bc_p_g[bcv_hi];

    delp_y = p_lo - p_hi;
    m_boundary_conditions.set_delp(1, delp_y);

    pj = p_hi;
  }
  else if(((bct_klo(dom_lo[0],dom_lo[1],dom_lo[2]-1,0) == BCList::pinf) &&
          (bct_khi(dom_lo[0],dom_lo[1],dom_hi[2]+1,0) == BCList::pout))
         ||
          ((bct_khi(dom_lo[0],dom_lo[1],dom_hi[2]+1,0) == BCList::pinf) &&
          (bct_klo(dom_lo[0],dom_lo[1],dom_lo[2]-1,0) == BCList::pout)))
  {
    delp_dir_loc = 2;

    const int bcv_lo = bct_klo(dom_lo[0],dom_lo[1],dom_lo[2]-1,1);
    const Real p_lo  = m_h_bc_p_g[bcv_lo];

    const int bcv_hi = bct_khi(dom_lo[0],dom_lo[1],dom_hi[2]+1,1);
    const Real p_hi  = m_h_bc_p_g[bcv_hi];

    delp_z = p_lo - p_hi;
    m_boundary_conditions.set_delp(2, delp_z);

    pj = p_hi;
  }

  // ---------------------------------------------------------------->>>
  //  Set default value of pj to zero in case no initial conditions are set
  pj = 0;

  //  Make sure that ic_p_g is set if using delp pressure conditions
  for(int icv(0); icv < m_initial_conditions.ic().size(); ++icv)
  {
    if((delp_dir_loc >= 0) && (delp_dir_loc == m_boundary_conditions.delp_dir()))
    {
      if (!m_initial_conditions.ic(icv).fluid.pressure_defined)
      {
        std::cout << "MUST DEFINE ic_p_g if using the DELP pressure condition" << std::endl;
        exit(0);
      }
      pj = m_initial_conditions.ic(icv).fluid.pressure;
    }
    else if((delp_dir_loc >= 0) && (delp_dir_loc != m_boundary_conditions.delp_dir()))
    {
      if( m_initial_conditions.ic(icv).fluid.pressure_defined)
      {
        std::cout << "MUST not define ic_p_g if setting p_inflowandp_outflow" << std::endl;
        exit(0);
      }
    }
    else
    {
      const Real gravity_square_module =
        gravity[0]*gravity[0] + gravity[1]*gravity[1] + gravity[2]*gravity[2];

      if( !m_initial_conditions.ic(icv).fluid.pressure || gravity_square_module > tolerance)
      {
        const Real ro_g0 = m_initial_conditions.ic(icv).fluid.density;

        compute_p0_bcs(sbx, domain, array4_p0_g, m_bc_p_g, pj, gravity, dx, dy,
                       dz, bct_ilo, bct_ihi, bct_jlo, bct_jhi, bct_klo, bct_khi,
                       nlft, nrgt, nbot, ntop, ndwn, nup, nghost_state(), ro_g0,
                       m_boundary_conditions);

        return;
      }

      const Real ic_p_g = m_initial_conditions.ic(icv).fluid.pressure;

      amrex::ParallelFor(sbx, [array4_p0_g,ic_p_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          { array4_p0_g(i,j,k) = ic_p_g; });
    }
  }

  // ---------------------------------------------------------------->>>

  // Here the pressure in each cell is determined from a specified pressure
  // drop across the domain le>. This section requires that the pressure
  // is already defined in all initial condition regions (otherwise this
  // section would be skipped)

  //  This hack allows to set the IC pressure  at L-dx/2 orboth
  //  nodalandCC pressure -> reference value orpressure, AKA IC_P_G,
  //  is set at the last cell center location.
  if(delp_dir_loc != m_boundary_conditions.delp_dir())
    offset = -1.;

  if(amrex::Math::abs(delp_x) > tolerance)
  {
    const Real dpodx = delp_x / xlen;
    pj -= dpodx * dx * (bx_hi[0] - dom_hi[0] + nghost_state() + 2 + offset);

    amrex::ParallelFor(sbx, [pj,dpodx,dx,sbx_hi,array4_p0_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      array4_p0_g(i,j,k) = pj + dpodx*dx * (sbx_hi[0] - i + 1);
    });

    pj += dpodx * dx * (sbx_hi[0] - sbx_lo[0] + 1);
  }

  if(amrex::Math::abs(delp_y) > tolerance)
  {
    const Real dpody = delp_y / ylen;
    pj -= dpody * dy * (bx_hi[1] - dom_hi[1] + nghost_state() + 2 + offset);

    amrex::ParallelFor(sbx, [pj,dpody,dy,sbx_hi,array4_p0_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          array4_p0_g(i,j,k) = pj + dpody*dy * (sbx_hi[1] - j + 1);
        });

    pj += dpody * dy * (sbx_hi[1] - sbx_lo[1] + 1);
  }

  if(amrex::Math::abs(delp_z) > tolerance)
  {
    const Real dpodz = delp_z / zlen;
    pj -= dpodz * dz * (bx_hi[2] - dom_hi[2] + nghost_state() + 2 + offset);

    amrex::ParallelFor(sbx, [pj,dpodz,dz,sbx_hi,array4_p0_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          array4_p0_g(i,j,k) = pj + (dpodz*dz * (sbx_hi[2] - k + 1));
        });

    pj += dpodz * dz * (sbx_hi[2] - sbx_lo[2] + 1);
  }

  // pressure in all initial condition region cells was defined
  set_p0_bcs(sbx, domain, array4_p0_g, m_bc_p_g, bct_ilo, bct_ihi, bct_jlo,
             bct_jhi, bct_klo, bct_khi, nlft, nrgt, nbot, ntop, ndwn, nup, nghost_state());

  return;
}


namespace set_p0_aux {

void compute_p0_bcs (const Box& sbx,
                     const Box& domain,
                     Array4<Real> const& p0_g,
                     Gpu::DeviceVector<Real>& m_bc_p_g,
                     Real& pj,
                     const RealVect& gravity,
                     const Real dx,
                     const Real dy,
                     const Real dz,
                     Array4<const int> const& bct_ilo,
                     Array4<const int> const& bct_ihi,
                     Array4<const int> const& bct_jlo,
                     Array4<const int> const& bct_jhi,
                     Array4<const int> const& bct_klo,
                     Array4<const int> const& bct_khi,
                     const int nlft,
                     const int nrgt,
                     const int nbot,
                     const int ntop,
                     const int ndwn,
                     const int nup,
                     const int nghost,
                     const Real ro_g0,
                     const MFIXBoundaryConditions& boundary_conditions)
{
  const Real tolerance = std::numeric_limits<Real>::epsilon();

  const IntVect dom_lo(domain.loVect());
  const IntVect dom_hi(domain.hiVect());

  const IntVect sbx_lo(sbx.loVect());
  const IntVect sbx_hi(sbx.hiVect());

  // Search oran outflow boundary condition where pressure is specified
  const Real undefined = 9.87654321e31;
  pj = undefined;

  for(int bcv(0); bcv < boundary_conditions.bc().size(); ++bcv)
  {
    const BC_t& bc = boundary_conditions.bc(bcv);

    if (bc.type == BCList::pout)
      pj = bc.fluid.pressure;
  }

  // Either a PO was not specified or PO was specified but not the
  // pressure at the outlet
  if ( std::abs(pj - undefined) < tolerance)
  {
    amrex::ParallelFor(sbx, [p0_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { p0_g(i,j,k) = 0; });

    pj = 0;
  }

  // ----------------------------------------------------------------<<<

  // Set an approximate pressure field assuming that the pressure drop
  // balances the weight of the bed, if the initial pressure-field is not
  // specified

  if(amrex::Math::abs(gravity[0]) > tolerance)
  {
    // Find the average weight per unit area over an x-z slice
    const Real dpodx = -gravity[0]*ro_g0;

    const int bx_lo_x = amrex::max(dom_lo[0], sbx_lo[0]);
    const int bx_hi_x = amrex::min(dom_hi[0]+1, sbx_hi[0]);

    const Box bx({bx_lo_x, sbx_lo[1], sbx_lo[2]},
                 {bx_hi_x, sbx_hi[1], sbx_hi[2]});

    const int upper_stride = (bx_hi_x == dom_hi[0]+1) ? 0 : (dom_hi[0]+1 - sbx_hi[0]);
    const int lower_stride = (bx_lo_x == dom_lo[0]) ? 0 : (sbx_lo[0] - dom_lo[0]);

    const int bx_delta_x = bx_hi_x - bx_lo_x + 1;

    if(gravity[0] < 0)
    {
      pj += upper_stride * dpodx * dx;

      amrex::ParallelFor(bx, [pj,dpodx,dx,bx_hi_x,p0_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            p0_g(i,j,k) = pj + dpodx*dx * (bx_hi_x - i);
          });

      pj += (bx_delta_x + lower_stride) * dpodx*dx;
    }
    else
    {
      pj -= lower_stride * dpodx * dx;

      amrex::ParallelFor(bx, [pj,dpodx,dx,bx_lo_x,p0_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            p0_g(i,j,k) = pj - dpodx*dx * (i - bx_lo_x);
          });

      pj -= (bx_delta_x + upper_stride) * dpodx*dx;
    }
  }
  else if (amrex::Math::abs(gravity[1]) > tolerance)
  {
    const Real dpody = -gravity[1]*ro_g0;

    const int bx_lo_y = amrex::max(dom_lo[1], sbx_lo[1]);
    const int bx_hi_y = amrex::min(dom_hi[1]+1, sbx_hi[1]);

    const Box bx({sbx_lo[0], bx_lo_y, sbx_lo[2]},
                 {sbx_hi[0], bx_hi_y, sbx_hi[2]});

    const int upper_stride = (bx_hi_y == dom_hi[1]+1) ? 0 : (dom_hi[1]+1 - sbx_hi[1]);
    const int lower_stride = (bx_lo_y == dom_lo[1]) ? 0 : (sbx_lo[1] - dom_lo[1]);

    const int bx_delta_y = bx_hi_y - bx_lo_y + 1;

    if(gravity[1] < 0)
    {
      pj += upper_stride * dpody * dy;

      amrex::ParallelFor(bx, [p0_g,dpody,dy,bx_hi_y, pj]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            p0_g(i,j,k) = pj + dpody*dy * (bx_hi_y - j);
          });

      pj += (bx_delta_y + lower_stride) * dpody*dy;
    }
    else
    {
      pj -= lower_stride * dpody * dy;

      amrex::ParallelFor(bx, [p0_g,dpody,dy,bx_lo_y,pj]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            p0_g(i,j,k) = pj - dpody*dy * (j - bx_lo_y);
          });

      pj -= (bx_delta_y + upper_stride) * dpody*dy;
    }
  }
  else if(amrex::Math::abs(gravity[2]) > tolerance)
  {
    const Real dpodz = -gravity[2]*ro_g0;

    const int bx_lo_z = amrex::max(dom_lo[2], sbx_lo[2]);
    const int bx_hi_z = amrex::min(dom_hi[2]+1, sbx_hi[2]);

    const Box bx({sbx_lo[0], sbx_lo[1], bx_lo_z},
                 {sbx_hi[0], sbx_hi[1], bx_hi_z});

    const int upper_stride = (bx_hi_z == dom_hi[2]+1) ? 0 : (dom_hi[2]+1 - sbx_hi[2]);
    const int lower_stride = (bx_lo_z == dom_lo[2]) ? 0 : (sbx_lo[2] - dom_lo[2]);

    const int bx_delta_z = bx_hi_z - bx_lo_z + 1;

    if(gravity[2] < 0)
    {
      pj += upper_stride * dpodz * dz;

      amrex::ParallelFor(bx, [p0_g,dpodz,dz,bx_hi_z,pj]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            p0_g(i,j,k) = pj + dpodz*dz * (bx_hi_z - k);
          });

      pj += (bx_delta_z + lower_stride) * dpodz*dz;
    }
    else
    {
      pj -= lower_stride * dpodz * dz;

      amrex::ParallelFor(bx, [p0_g,dpodz,dz,bx_lo_z,pj]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            p0_g(i,j,k) = pj - dpodz*dz * (k - bx_lo_z);
          });

      pj -= (bx_delta_z + upper_stride) * dpodz*dz;
    }
  }

  set_p0_bcs(sbx, domain, p0_g, m_bc_p_g, bct_ilo, bct_ihi, bct_jlo, bct_jhi,
             bct_klo, bct_khi, nlft, nrgt, nbot, ntop, ndwn, nup, nghost);

  return;
}

void set_p0_bcs (const Box& sbx,
                 const Box& domain,
                 Array4<Real> const& p0_g,
                 Gpu::DeviceVector<Real>& m_bc_p_g,
                 Array4<const int> const& bct_ilo,
                 Array4<const int> const& bct_ihi,
                 Array4<const int> const& bct_jlo,
                 Array4<const int> const& bct_jhi,
                 Array4<const int> const& bct_klo,
                 Array4<const int> const& bct_khi,
                 const int nlft,
                 const int nrgt,
                 const int nbot,
                 const int ntop,
                 const int ndwn,
                 const int nup,
                 const int nghost)
{
  const IntVect sbx_lo(sbx.loVect());
  const IntVect sbx_hi(sbx.hiVect());

  const IntVect dom_lo(domain.loVect());
  const IntVect dom_hi(domain.hiVect());

  Real* m_bc_p_g_ptr = m_bc_p_g.data();

  if (nlft > 0)
  {
    const Box sbx_lo_x(sbx_lo, {dom_lo[0], sbx_hi[1], sbx_hi[2]});

    amrex::ParallelFor(sbx_lo_x,
        [dom_hi,nghost,bct_ilo,dom_lo,m_bc_p_g_ptr,p0_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int jbc = (j > dom_hi[1]+nghost) ? j-1 : j;
          const int kbc = (k > dom_hi[2]+nghost) ? k-1 : k;

          const int bct = bct_ilo(dom_lo[0]-1, jbc, kbc, 0);

          if(bct == BCList::pinf || bct == BCList::pout)
          {
            const int bcv = bct_ilo(dom_lo[0]-1, jbc, kbc, 1);
            p0_g(i,j,k) = m_bc_p_g_ptr[bcv];
          }
        });
  }

  if (nrgt > 0)
  {
    const Box sbx_hi_x({dom_hi[0]+1, sbx_lo[1], sbx_lo[2]}, sbx_hi);

    amrex::ParallelFor(sbx_hi_x,
        [dom_hi,nghost,bct_ihi,m_bc_p_g_ptr,p0_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int jbc = (j > dom_hi[1]+nghost) ? j-1 : j;
          const int kbc = (k > dom_hi[2]+nghost) ? k-1 : k;

          const int bct = bct_ihi(dom_hi[0]+1, jbc, kbc, 0);

          if(bct == BCList::pinf || bct == BCList::pout)
          {
            const int bcv = bct_ihi(dom_hi[0]+1, jbc, kbc, 1);
            p0_g(i,j,k) = m_bc_p_g_ptr[bcv];
          }
        });
  }

  if (nbot > 0)
  {
    const Box sbx_lo_y(sbx_lo, {sbx_hi[0], dom_lo[1], sbx_hi[2]});

    amrex::ParallelFor(sbx_lo_y,
        [dom_hi,nghost,bct_jlo,dom_lo,m_bc_p_g_ptr,p0_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int ibc = (i > dom_hi[0]+nghost) ? i-1 : i;
          const int kbc = (k > dom_hi[2]+nghost) ? k-1 : k;

          const int bct = bct_jlo(ibc, dom_lo[1]-1, kbc, 0);

          if(bct == BCList::pinf || bct == BCList::pout)
          {
            const int bcv = bct_jlo(ibc, dom_lo[1]-1, kbc, 1);
            p0_g(i,j,k) = m_bc_p_g_ptr[bcv];
          }
        });
  }

  if (ntop > 0)
  {
    const Box sbx_hi_y({sbx_lo[0], dom_hi[1]+1, sbx_lo[2]}, sbx_hi);

    amrex::ParallelFor(sbx_hi_y,
        [dom_hi,nghost,bct_jhi,m_bc_p_g_ptr,p0_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int ibc = (i > dom_hi[0]+nghost) ? i-1 : i;
          const int kbc = (k > dom_hi[2]+nghost) ? k-1 : k;

          const int bct = bct_jhi(ibc, dom_hi[1]+1, kbc, 0);

          if(bct == BCList::pinf || bct == BCList::pout)
          {
            const int bcv = bct_jhi(ibc, dom_hi[1]+1, kbc, 1);
            p0_g(i,j,k) = m_bc_p_g_ptr[bcv];
          }
        });
  }

  if (ndwn > 0)
  {
    const Box sbx_lo_z(sbx_lo, {sbx_hi[0], sbx_hi[1], dom_lo[2]});

    amrex::ParallelFor(sbx_lo_z,
        [dom_hi,nghost,bct_klo,dom_lo,m_bc_p_g_ptr,p0_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int ibc = (i > dom_hi[0]+nghost) ? i-1 : i;
          const int jbc = (j > dom_hi[1]+nghost) ? j-1 : j;

          const int bct = bct_klo(ibc, jbc, dom_lo[2]-1, 0);

          if(bct == BCList::pinf || bct == BCList::pout)
          {
            const int bcv = bct_klo(ibc, jbc, dom_lo[2]-1, 1);
            p0_g(i,j,k) = m_bc_p_g_ptr[bcv];
          }
        });
  }

  if (nup > 0)
  {
    const Box sbx_hi_z({sbx_lo[0], sbx_lo[1], dom_hi[2]+1}, sbx_hi);

    amrex::ParallelFor(sbx_hi_z,
        [dom_hi,nghost,bct_khi,m_bc_p_g_ptr,p0_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int ibc = (i > dom_hi[0]+nghost) ? i-1 : i;
          const int jbc = (j > dom_hi[1]+nghost) ? j-1 : j;

          const int bct = bct_khi(ibc, jbc, dom_hi[2]+1, 0);

          if(bct == BCList::pinf || bct == BCList::pout)
          {
            const int bcv = bct_khi(ibc, jbc, dom_hi[2]+1, 1);
            p0_g(i,j,k) = m_bc_p_g_ptr[bcv];
          }
        });
  }

  return;
}

} // end namespace set_p0_aux
