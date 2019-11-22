#include <mfix.H>
#include <bc_mod_F.H>
#include <ic_mod_F.H>
#include <climits>
#include <scales_mod_F.H>
#include <constant_mod_F.H>
#include <param_mod_F.H>
#include <constant_mod_F.H>
#include <mfix_des_F.H>
#include <fld_constants_mod_F.H>

#include <string>

namespace set_p0_aux {

void goto_60(const Box& sbx, const Box& domain, const BcList& bc_list,
             Array4<Real> const& p0_g, Real* m_bc_p_g, Real& pj,
             const Real* gravity, const Real dx, const Real dy, const Real dz,
             Array4<const int> const& bct_ilo,
             Array4<const int> const& bct_ihi,
             Array4<const int> const& bct_jlo,
             Array4<const int> const& bct_jhi,
             Array4<const int> const& bct_klo,
             Array4<const int> const& bct_khi,
             const int nlft, const int nrgt, const int nbot,
             const int ntop, const int ndwn, const int nup,
             const int nghost);

void goto_100(const Box& sbx, const Box& domain, const BcList& bc_list,
              Array4<Real> const& p0_g, Real* m_bc_p_g,
              Array4<const int> const& bct_ilo,
              Array4<const int> const& bct_ihi,
              Array4<const int> const& bct_jlo,
              Array4<const int> const& bct_jhi,
              Array4<const int> const& bct_klo,
              Array4<const int> const& bct_khi,
              const int nlft, const int nrgt, const int nbot,
              const int ntop, const int ndwn, const int nup,
              const int nghost);

} // end namespace set_p0_aux


using namespace set_p0_aux;


void
mfix::set_p0(const Box& bx,
             MFIter* mfi,
             const int lev,
             const Box& domain,
             const Real xlen,
             const Real ylen,
             const Real zlen,
             const int delp_dir_in)
{
  Real dx = geom[lev].CellSize(0);
  Real dy = geom[lev].CellSize(1);
  Real dz = geom[lev].CellSize(2);

  const Real tolerance = std::numeric_limits<Real>::epsilon();
  Real offset(-0.5);

  Array4<double> const& array4_p0_g = p0_g[lev]->array(*mfi);

  const IntVect sbx_lo((*p0_g[lev])[*mfi].loVect());
  const IntVect sbx_hi((*p0_g[lev])[*mfi].hiVect());
  const Box sbx(sbx_lo, sbx_hi);

  const IntVect dom_lo(domain.loVect());
  const IntVect dom_hi(domain.hiVect());

  const IntVect bx_lo(bx.loVect());
  const IntVect bx_hi(bx.hiVect());

  Array4<const int> const& bct_ilo = bc_ilo[lev]->array();
  Array4<const int> const& bct_ihi = bc_ihi[lev]->array();
  Array4<const int> const& bct_jlo = bc_jlo[lev]->array();
  Array4<const int> const& bct_jhi = bc_jhi[lev]->array();
  Array4<const int> const& bct_klo = bc_klo[lev]->array();
  Array4<const int> const& bct_khi = bc_khi[lev]->array();

  const int nlft = std::max(0, dom_lo[0]-sbx_lo[0]+1);
  const int nbot = std::max(0, dom_lo[1]-sbx_lo[1]+1);
  const int ndwn = std::max(0, dom_lo[2]-sbx_lo[2]+1);

  const int nrgt = std::max(0, sbx_hi[0]-dom_hi[0]);
  const int ntop = std::max(0, sbx_hi[1]-dom_hi[1]);
  const int nup  = std::max(0, sbx_hi[2]-dom_hi[2]);

  int delp_dir(delp_dir_in);

  Real pj(0);

  Real delp_x = get_delp_x();
  Real delp_y = get_delp_y();
  Real delp_z = get_delp_z();

  // ---------------------------------------------------------------->>>
  //     If the bc's are pressure inflow/outflow then be sure to capture that in p0andgp0
  // ---------------------------------------------------------------->>>

  const amrex::Real P_ref = get_p_ref();
  const amrex::Real P_scale = get_p_scale();

  if(((bct_ilo(dom_lo[0]-1,dom_lo[1],dom_lo[2],0) == bc_list.get_pinf()) and
     (bct_ihi(dom_hi[0]+1,dom_lo[1],dom_lo[2],0) == bc_list.get_pout()))
    or
     ((bct_ihi(dom_hi[0]+1,dom_lo[1],dom_lo[2],0) == bc_list.get_pinf()) and
     (bct_ilo(dom_lo[0]-1,dom_lo[1],dom_lo[2],0) == bc_list.get_pout())))
  {
    delp_dir = 0;

    const int bcv_lo = bct_ilo(dom_lo[0]-1,dom_lo[1],dom_lo[2],1);
    const Real p_lo  = scale_pressure_cpp(m_bc_p_g[bcv_lo], P_ref, P_scale);

    const int bcv_hi = bct_ihi(dom_hi[0]+1,dom_lo[1],dom_lo[2],1);
    const Real p_hi  = scale_pressure_cpp(m_bc_p_g[bcv_hi], P_ref, P_scale);

    delp_x = p_lo - p_hi;
    set_delp_x(delp_x);

    pj = p_hi;
  }
  else if(((bct_jlo(dom_lo[0],dom_lo[1]-1,dom_lo[2],0) == bc_list.get_pinf()) and
          (bct_jhi(dom_lo[0],dom_hi[1]+1,dom_lo[2],0) == bc_list.get_pout()))
         or
          ((bct_jhi(dom_lo[0],dom_hi[1]+1,dom_lo[2],0) == bc_list.get_pinf()) and
          (bct_jlo(dom_lo[0],dom_lo[1]-1,dom_lo[2],0) == bc_list.get_pout())))
  {
    delp_dir = 1;

    const int bcv_lo = bct_jlo(dom_lo[0],dom_lo[1]-1,dom_lo[2],1);
    const Real p_lo  = scale_pressure_cpp(m_bc_p_g[bcv_lo], P_ref, P_scale);

    const int bcv_hi = bct_jhi(dom_lo[0],dom_hi[1]+1,dom_lo[2],1);
    const Real p_hi  = scale_pressure_cpp(m_bc_p_g[bcv_hi], P_ref, P_scale);

    delp_y = p_lo - p_hi;
    set_delp_y(delp_y);

    pj = p_hi;
  }
  else if(((bct_klo(dom_lo[0],dom_lo[1],dom_lo[2]-1,0) == bc_list.get_pinf()) and
          (bct_khi(dom_lo[0],dom_lo[1],dom_hi[2]+1,0) == bc_list.get_pout()))
         or
          ((bct_khi(dom_lo[0],dom_lo[1],dom_hi[2]+1,0) == bc_list.get_pinf()) and
          (bct_klo(dom_lo[0],dom_lo[1],dom_lo[2]-1,0) == bc_list.get_pout())))
  {
    delp_dir = 2;

    const int bcv_lo = bct_klo(dom_lo[0],dom_lo[1],dom_lo[2]-1,1);
    const Real p_lo  = scale_pressure_cpp(m_bc_p_g[bcv_lo], P_ref, P_scale);

    const int bcv_hi = bct_khi(dom_lo[0],dom_lo[1],dom_hi[2]+1,1);
    const Real p_hi  = scale_pressure_cpp(m_bc_p_g[bcv_hi], P_ref, P_scale);

    delp_z = p_lo - p_hi;
    set_delp_z(delp_z);

    pj = p_hi;
  }

  // ---------------------------------------------------------------->>>
  //  Set default value of pj to zero in case no initial conditions are set
  pj = 0;

  //  Make sure that ic_p_g is set if using delp pressure conditions
  for(int icv(1); icv <= get_dim_ic(); ++icv)
  {
    if(ic_defined_cpp(icv))
    {
      if((delp_dir >= 0) and (delp_dir == delp_dir_in))
      {
        if (not is_defined_db_cpp(get_ic_p_g(icv)))
        {
          std::cout << "MUST DEFINE ic_p_g if using the DELP pressure condition" << std::endl;
          exit(0);
        }
        pj = get_ic_p_g(icv);
      }
      else if((delp_dir >= 0) and (delp_dir != delp_dir_in))
      {
        if(is_defined_db_cpp(get_ic_p_g(icv)))
        {
          std::cout << "MUST not define ic_p_g if setting p_inflowandp_outflow" << std::endl;
          exit(0);
        }
      }
      else
      {
        const amrex::Real gravity_square_module =
          gravity[0]*gravity[0] + gravity[1]*gravity[1] + gravity[2]*gravity[2];

        if(is_undefined_db_cpp(get_ic_p_g(icv)) or gravity_square_module > tolerance)
        {
          goto_60(sbx, domain, bc_list, array4_p0_g, m_bc_p_g.data(), pj,
              gravity, dx, dy, dz, bct_ilo, bct_ihi, bct_jlo, bct_jhi, bct_klo,
              bct_khi, nlft, nrgt, nbot, ntop, ndwn, nup, nghost);
          return;
        }

        const Real ic_p_g = get_ic_p_g(icv);

        amrex::ParallelFor(sbx, [array4_p0_g,ic_p_g]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            { array4_p0_g(i,j,k) = ic_p_g; });
      }
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
  if(delp_dir != delp_dir_in)
    offset = -1.;

  if(std::abs(delp_x) > tolerance)
  {
    const amrex::Real dpodx = delp_x / xlen;
    pj -= dpodx * dx * (bx_hi[0] - dom_hi[0] + nghost + 2 + offset);

    amrex::ParallelFor(sbx, [pj,dpodx,dx,sbx_hi,array4_p0_g,P_ref,P_scale]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const amrex::Real local_pj = pj + dpodx*dx * (sbx_hi[0] - i + 1);
      array4_p0_g(i,j,k) = scale_pressure_cpp(local_pj, P_ref, P_scale);
    });

    pj += dpodx * dx * (sbx_hi[0] - sbx_lo[0] + 1);
  }

  if(std::abs(delp_y) > tolerance)
  {
    const Real dpody = delp_y / ylen;
    pj -= dpody * dy * (bx_hi[1] - dom_hi[1] + nghost + 2 + offset);

    amrex::ParallelFor(sbx, [pj,dpody,dy,sbx_hi,array4_p0_g,P_ref,P_scale]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const amrex::Real local_pj = pj + dpody*dy * (sbx_hi[1] - j + 1);
          array4_p0_g(i,j,k) = scale_pressure_cpp(local_pj, P_ref, P_scale);
        });

    pj += dpody * dy * (sbx_hi[1] - sbx_lo[1] + 1);
  }

  if(std::abs(delp_z) > tolerance)
  {
    const Real dpodz = delp_z / zlen;
    pj -= dpodz * dz * (bx_hi[2] - dom_hi[2] + nghost + 2 + offset);

    amrex::ParallelFor(sbx, [pj,dpodz,dz,sbx_hi,array4_p0_g,P_ref,P_scale]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const amrex::Real local_pj = pj + (dpodz*dz * (sbx_hi[2] - k + 1));
          array4_p0_g(i,j,k) = scale_pressure_cpp(local_pj, P_ref, P_scale);
        });

    pj += dpodz * dz * (sbx_hi[2] - sbx_lo[2] + 1);
  }

  // pressure in all intial condition region cells was defined
  goto_100(sbx, domain, bc_list, array4_p0_g, m_bc_p_g.data(), bct_ilo, bct_ihi,
           bct_jlo, bct_jhi, bct_klo, bct_khi, nlft, nrgt, nbot, ntop, ndwn, nup,
           nghost);

  return;
}


namespace set_p0_aux {

void goto_60(const Box& sbx,
             const Box& domain,
             const BcList& bc_list,
             Array4<Real> const& p0_g,
             Real* m_bc_p_g,
             Real& pj,
             const Real* gravity,
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
             const int nghost)
{
  const Real tolerance = std::numeric_limits<Real>::epsilon();

  const IntVect dom_lo(domain.loVect());
  const IntVect dom_hi(domain.hiVect());

  const IntVect sbx_lo(sbx.loVect());
  const IntVect sbx_hi(sbx.hiVect());

  // Search oran outflow boundary condition where pressure is specified
  pj = get_undefined();

  for(int icv(1); icv <= get_dim_bc(); ++icv)
  {
    if(get_bc_defined(icv))
    {
      char bc_type[16];
      get_bc_type(icv, bc_type);

      std::string bc_type_str(bc_type);

      if(bc_type_str.compare("P_OUTFLOW") == 0 or bc_type_str.compare("PO") == 0)
      {
        pj = m_bc_p_g[icv];  
      }
    }
  }

  // Either a PO was not specified or PO was specified but not the
  // pressure at the outlet
  if (is_undefined_db_cpp(pj))
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

  const amrex::Real P_ref = get_p_ref();
  const amrex::Real P_scale = get_p_scale();

  const Real ro_g0 = get_ro_g0();

  if(std::abs(gravity[0]) > tolerance)
  {
    // Find the average weight per unit area over an x-z slice
    const Real dpodx = -gravity[0]*ro_g0;

    const int bx_lo_x = std::max(dom_lo[0], sbx_lo[0]);
    const int bx_hi_x = std::min(dom_hi[0]+1, sbx_hi[0]);

    const Box bx({bx_lo_x, sbx_lo[1], sbx_lo[2]},
                 {bx_hi_x, sbx_hi[1], sbx_hi[2]});

    const int upper_stride = (bx_hi_x == dom_hi[0]+1) ? 0 : (dom_hi[0]+1 - sbx_hi[0]);
    const int lower_stride = (bx_lo_x == dom_lo[0]) ? 0 : (sbx_lo[0] - dom_lo[0]);

    const int bx_delta_x = bx_hi_x - bx_lo_x + 1;

    if(gravity[0] < 0)
    {
      pj += upper_stride * dpodx * dx;

      amrex::ParallelFor(bx, [pj,dpodx,dx,bx_hi_x,p0_g,P_ref,P_scale]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const amrex::Real local_pj = pj + dpodx*dx * (bx_hi_x - i);
            p0_g(i,j,k) = scale_pressure_cpp(local_pj, P_ref, P_scale);
          });

      pj += (bx_delta_x + lower_stride) * dpodx*dx;
    }
    else
    {
      pj -= lower_stride * dpodx * dx;

      amrex::ParallelFor(bx, [pj,dpodx,dx,bx_lo_x,p0_g,P_ref,P_scale]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const amrex::Real local_pj = pj - dpodx*dx * (i - bx_lo_x);
            p0_g(i,j,k) = scale_pressure_cpp(local_pj, P_ref, P_scale);
          });

      pj -= (bx_delta_x + upper_stride) * dpodx*dx;
    }
  }
  else if (abs(gravity[1]) > tolerance)
  {
    const Real dpody = -gravity[1]*ro_g0;

    const int bx_lo_y = std::max(dom_lo[1], sbx_lo[1]);
    const int bx_hi_y = std::min(dom_hi[1]+1, sbx_hi[1]);

    const Box bx({sbx_lo[0], bx_lo_y, sbx_lo[2]},
                 {sbx_hi[0], bx_hi_y, sbx_hi[2]});

    const int upper_stride = (bx_hi_y == dom_hi[1]+1) ? 0 : (dom_hi[1]+1 - sbx_hi[1]);
    const int lower_stride = (bx_lo_y == dom_lo[1]) ? 0 : (sbx_lo[1] - dom_lo[1]);

    const int bx_delta_y = bx_hi_y - bx_lo_y + 1;

    if(gravity[1] < 0)
    {
      pj += upper_stride * dpody * dy;

      amrex::ParallelFor(bx, [p0_g,dpody,dy,bx_hi_y,P_ref,P_scale,pj]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const amrex::Real local_pj = pj + dpody*dy * (bx_hi_y - j);
            p0_g(i,j,k) = scale_pressure_cpp(local_pj, P_ref, P_scale);
          });
      
      pj += (bx_delta_y + lower_stride) * dpody*dy;
    }
    else
    {
      pj -= lower_stride * dpody * dy;

      amrex::ParallelFor(bx, [p0_g,dpody,dy,bx_lo_y,P_ref,P_scale,pj]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const amrex::Real local_pj = pj - dpody*dy * (j - bx_lo_y);
            p0_g(i,j,k) = scale_pressure_cpp(local_pj, P_ref, P_scale);
          });
      
      pj -= (bx_delta_y + upper_stride) * dpody*dy;
    }
  }
  else if(std::abs(gravity[2]) > tolerance)
  {
    const Real dpodz = -gravity[2]*ro_g0;

    const int bx_lo_z = std::max(dom_lo[2], sbx_lo[2]);
    const int bx_hi_z = std::min(dom_hi[2]+1, sbx_hi[2]);

    const Box bx({sbx_lo[0], sbx_lo[1], bx_lo_z},
                 {sbx_hi[0], sbx_hi[1], bx_hi_z});

    const int upper_stride = (bx_hi_z == dom_hi[2]+1) ? 0 : (dom_hi[2]+1 - sbx_hi[2]);
    const int lower_stride = (bx_lo_z == dom_lo[2]) ? 0 : (sbx_lo[2] - dom_lo[2]);

    const int bx_delta_z = bx_hi_z - bx_lo_z + 1;

    if(gravity[2] < 0)
    {
      pj += upper_stride * dpodz * dz;

      amrex::ParallelFor(bx, [p0_g,dpodz,dz,bx_hi_z,P_ref,P_scale,pj]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const amrex::Real local_pj = pj + dpodz*dz * (bx_hi_z - k);
            p0_g(i,j,k) = scale_pressure_cpp(local_pj, P_ref, P_scale);
          });

      pj += (bx_delta_z + lower_stride) * dpodz*dz;
    }
    else
    {
      pj -= lower_stride * dpodz * dz;

      amrex::ParallelFor(bx, [p0_g,dpodz,dz,bx_lo_z,P_ref,P_scale,pj]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const amrex::Real local_pj = pj - dpodz*dz * (k - bx_lo_z);
            p0_g(i,j,k) = scale_pressure_cpp(local_pj, P_ref, P_scale);
          });
      
      pj -= (bx_delta_z + upper_stride) * dpodz*dz;
    }
  }

  goto_100(sbx, domain, bc_list, p0_g, m_bc_p_g, bct_ilo, bct_ihi,
           bct_jlo, bct_jhi, bct_klo, bct_khi, nlft, nrgt, nbot, ntop, ndwn, nup,
           nghost);
  
  return;
}

void goto_100(const Box& sbx,
              const Box& domain,
              const BcList& bc_list,
              Array4<Real> const& p0_g,
              Real* m_bc_p_g,
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

  const amrex::Real P_ref = get_p_ref();
  const amrex::Real P_scale = get_p_scale();

  if (nlft > 0)
  {
    const Box sbx_lo_x(sbx_lo, {dom_lo[0], sbx_hi[1], sbx_hi[2]});

    amrex::ParallelFor(sbx_lo_x,
        [dom_hi,nghost,bct_ilo,bc_list,dom_lo,P_ref,P_scale,m_bc_p_g,p0_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int jbc = (j > dom_hi[1]+nghost) ? j-1 : j;
          const int kbc = (k > dom_hi[2]+nghost) ? k-1 : k;

          const int bct = bct_ilo(dom_lo[0]-1, jbc, kbc, 0);

          if(bct == bc_list.get_pinf() or bct == bc_list.get_pout())
          {
            const int bcv = bct_ilo(dom_lo[0]-1, jbc, kbc, 1);
            p0_g(i,j,k) = scale_pressure_cpp(m_bc_p_g[bcv], P_ref, P_scale);
          }
        });
  }

  if (nrgt > 0)
  {
    const Box sbx_hi_x({dom_hi[0]+1, sbx_lo[1], sbx_lo[2]}, sbx_hi);

    amrex::ParallelFor(sbx_hi_x,
        [dom_hi,nghost,bct_ihi,bc_list,dom_lo,P_ref,P_scale,m_bc_p_g,p0_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int jbc = (j > dom_hi[1]+nghost) ? j-1 : j;
          const int kbc = (k > dom_hi[2]+nghost) ? k-1 : k;

          const int bct = bct_ihi(dom_hi[0]+1, jbc, kbc, 0);

          if(bct == bc_list.get_pinf() or bct == bc_list.get_pout())
          {
            const int bcv = bct_ihi(dom_hi[0]+1, jbc, kbc, 1);
            p0_g(i,j,k) = scale_pressure_cpp(m_bc_p_g[bcv], P_ref, P_scale);
          }
        });
  }

  if (nbot > 0)
  {
    const Box sbx_lo_y(sbx_lo, {sbx_hi[0], dom_lo[1], sbx_hi[2]});

    amrex::ParallelFor(sbx_lo_y,
        [dom_hi,nghost,bct_jlo,bc_list,dom_lo,P_ref,P_scale,m_bc_p_g,p0_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int ibc = (i > dom_hi[0]+nghost) ? i-1 : i;
          const int kbc = (k > dom_hi[2]+nghost) ? k-1 : k;

          const int bct = bct_jlo(ibc, dom_lo[1]-1, kbc, 0);

          if(bct == bc_list.get_pinf() or bct == bc_list.get_pout())
          {
            const int bcv = bct_jlo(ibc, dom_lo[1]-1, kbc, 1);
            p0_g(i,j,k) = scale_pressure_cpp(m_bc_p_g[bcv], P_ref, P_scale);
          }
        });
  }

  if (ntop > 0)
  {
    const Box sbx_hi_y({sbx_lo[0], dom_hi[1]+1, sbx_lo[2]}, sbx_hi);

    amrex::ParallelFor(sbx_hi_y,
        [dom_hi,nghost,bct_jhi,bc_list,dom_lo,P_ref,P_scale,m_bc_p_g,p0_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int ibc = (i > dom_hi[0]+nghost) ? i-1 : i;
          const int kbc = (k > dom_hi[2]+nghost) ? k-1 : k;

          const int bct = bct_jhi(ibc, dom_hi[1]+1, kbc, 0);

          if(bct == bc_list.get_pinf() or bct == bc_list.get_pout())
          {
            const int bcv = bct_jhi(ibc, dom_hi[1]+1, kbc, 1);
            p0_g(i,j,k) = scale_pressure_cpp(m_bc_p_g[bcv], P_ref, P_scale);
          }
        });
  }

  if (ndwn > 0)
  {
    const Box sbx_lo_z(sbx_lo, {sbx_hi[0], sbx_hi[1], dom_lo[2]});

    amrex::ParallelFor(sbx_lo_z,
        [dom_hi,nghost,bct_klo,bc_list,dom_lo,P_ref,P_scale,m_bc_p_g,p0_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int ibc = (i > dom_hi[0]+nghost) ? i-1 : i;
          const int jbc = (j > dom_hi[1]+nghost) ? j-1 : j;

          const int bct = bct_klo(ibc, jbc, dom_lo[2]-1, 0);

          if(bct == bc_list.get_pinf() or bct == bc_list.get_pout())
          {
            const int bcv = bct_klo(ibc, jbc, dom_lo[2]-1, 1);
            p0_g(i,j,k) = scale_pressure_cpp(m_bc_p_g[bcv], P_ref, P_scale);
          }
        });
  }

  if (nup > 0)
  {
    const Box sbx_hi_z({sbx_lo[0], sbx_lo[1], dom_hi[2]+1}, sbx_hi);

    amrex::ParallelFor(sbx_hi_z,
        [dom_hi,nghost,bct_khi,bc_list,dom_lo,P_ref,P_scale,m_bc_p_g,p0_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int ibc = (i > dom_hi[0]+nghost) ? i-1 : i;
          const int jbc = (j > dom_hi[1]+nghost) ? j-1 : j;

          const int bct = bct_khi(ibc, jbc, dom_hi[2]+1, 0);

          if(bct == bc_list.get_pinf() or bct == bc_list.get_pout())
          {
            const int bcv = bct_khi(ibc, jbc, dom_hi[2]+1, 1);
            p0_g(i,j,k) = scale_pressure_cpp(m_bc_p_g[bcv], P_ref, P_scale);
          }
        });
  }

  return;
}

} // end namespace set_p0_aux
