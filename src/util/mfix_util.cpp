#include <mfix.H>

void
mfix::check_for_nans (int lev)
{
  bool ug_has_nans = m_leveldata[lev]->vel_g->contains_nan(0);
  bool vg_has_nans = m_leveldata[lev]->vel_g->contains_nan(1);
  bool wg_has_nans = m_leveldata[lev]->vel_g->contains_nan(2);
  bool pg_has_nans = m_leveldata[lev]->p_g->contains_nan(0);

  if (ug_has_nans)
    amrex::Print() << "WARNING: u_g contains NaNs!!!";

  if (vg_has_nans)
    amrex::Print() << "WARNING: v_g contains NaNs!!!";

  if (wg_has_nans)
    amrex::Print() << "WARNING: w_g contains NaNs!!!";

  if (pg_has_nans)
    amrex::Print() << "WARNING: p_g contains NaNs!!!";
}

//
// Print the maximum values of the velocity components
//
void
mfix::mfix_print_max_vel (int lev)
{
    amrex::Print() << "   max(abs(u/v/w/p))  = "
                   << m_leveldata[lev]->vel_g->norm0(0,0,false,true) << "  "
                   << m_leveldata[lev]->vel_g->norm0(1,0,false,true) << "  "
                   << m_leveldata[lev]->vel_g->norm0(2,0,false,true) << "  "
                   << m_leveldata[lev]->p_g->norm0(0,0,false,true) << std::endl;
}


//
// Print the maximum values of the pressure gradient components
//
void
mfix::mfix_print_max_gp (int lev)
{
    amrex::Print() << "   max(abs(gpx/gpy/gpz))  = "
                   << m_leveldata[lev]->gp->norm0(0,0,false,true) << "  "
                   << m_leveldata[lev]->gp->norm0(1,0,false,true) << "  "
                   << m_leveldata[lev]->gp->norm0(2,0,false,true) <<  std::endl;
}


void
mfix::mfix_compute_vort ()
{
    BL_PROFILE("mfix::mfix_compute_vort");

    for (int lev = 0; lev < nlev; lev++)
    {
      Box domain(geom[lev].Domain());

      m_leveldata[lev]->vort->setVal(0.0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter mfi(*m_leveldata[lev]->vel_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
          // Tilebox
          Box bx = mfi.tilebox ();
          const Real* dx = geom[lev].CellSize();

          Array4<Real> const& vorticity = m_leveldata[lev]->vort->array(mfi);
          Array4<Real> const& velocity_g = m_leveldata[lev]->vel_g->array(mfi);

          const Real odx(1./dx[0]), ody(1./dx[1]), odz(1./dx[2]);

          // This is to check efficiently if this tile contains any eb stuff
          const EBFArrayBox& vel_fab =
            static_cast<EBFArrayBox const&>((*m_leveldata[lev]->vel_g)[mfi]);

          const EBCellFlagFab& flags = vel_fab.getEBCellFlagFab();

          if (flags.getType(amrex::grow(bx,0)) == FabType::regular)
          {
            amrex::ParallelFor(bx, [odx,ody,odz,velocity_g,vorticity]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              Real uy = .5*ody*(velocity_g(i,j+1,k,0) - velocity_g(i,j-1,k,0));
              Real uz = .5*odz*(velocity_g(i,j,k+1,0) - velocity_g(i,j,k-1,0));
              Real vx = .5*odx*(velocity_g(i+1,j,k,1) - velocity_g(i-1,j,k,1));
              Real vz = .5*odz*(velocity_g(i,j,k+1,1) - velocity_g(i,j,k-1,1));
              Real wx = .5*odx*(velocity_g(i+1,j,k,2) - velocity_g(i-1,j,k,2));
              Real wy = .5*ody*(velocity_g(i,j+1,k,2) - velocity_g(i,j-1,k,2));

              vorticity(i,j,k) = std::sqrt((wy-vz)*(wy-vz) +
                                           (uz-wx)*(uz-wx) +
                                           (vx-uy)*(vx-uy));
            });
          }
       }
    }
}

Real
mfix::volWgtSum (int lev, const MultiFab& mf, int comp, bool local)
{
    BL_PROFILE("mfix::volWgtSum()");

    const MultiFab* volfrac =  &(ebfactory[lev]->getVolFrac());

    Real sum = amrex::ReduceSum(mf, *volfrac, 0,
        [comp] AMREX_GPU_HOST_DEVICE (Box const & bx,
                                      Array4<const Real> const rho,
                                      Array4<const Real> const vfrc)
        {
          Real dm = 0.0;

          amrex::Loop(bx, [rho,vfrc,comp,&dm] (int i, int j, int k) noexcept
              { dm += rho(i,j,k,comp) * vfrc(i,j,k); });
          
          return dm;
        });   

    if (!local)
        ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
mfix::volEpsWgtSum (int lev, const MultiFab& mf, int comp, bool local)
{
    BL_PROFILE("mfix::volEpsWgtSum()");

    const MultiFab* volfrac =  &(ebfactory[lev]->getVolFrac());

    Real sum = amrex::ReduceSum(mf, *volfrac, *(m_leveldata[lev]->ep_g), 0,
        [comp] AMREX_GPU_HOST_DEVICE (Box const & bx,
                                      Array4<const Real> const rho,
                                      Array4<const Real> const vfrc,
                                      Array4<const Real> const ep)
        {
          Real dm = 0.0;

          amrex::Loop(bx, [rho,vfrc,ep,comp,&dm] (int i, int j, int k) noexcept
              { dm += rho(i,j,k,comp) * vfrc(i,j,k) * ep(i,j,k); });
          
          return dm;
        });   

    if (!local)
        ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}
