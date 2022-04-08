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
mfix::mfix_print_max_vel (int lev,
                          const Vector<MultiFab*>& vel_g_in,
                          const Vector<MultiFab*>& p_g_in)
{
    amrex::Print() << "   max(abs(u/v/w/p))  = "
                   << vel_g_in[lev]->norm0(0,0,false,true) << "  "
                   << vel_g_in[lev]->norm0(1,0,false,true) << "  "
                   << vel_g_in[lev]->norm0(2,0,false,true) << "  "
                   << p_g_in[lev]->norm0(0,0,false,true) << std::endl;
}


//
// Print the maximum values of the pressure gradient components
//
void
mfix::mfix_print_max_gp (int lev,
                         const Vector<MultiFab*>& gp_g_in)
{
    amrex::Print() << "   max(abs(gpx/gpy/gpz))  = "
                   << gp_g_in[lev]->norm0(0,0,false,true) << "  "
                   << gp_g_in[lev]->norm0(1,0,false,true) << "  "
                   << gp_g_in[lev]->norm0(2,0,false,true) <<  std::endl;
}


void
mfix::mfix_compute_vort ()
{
    BL_PROFILE("mfix::mfix_compute_vort");

    for (int lev = 0; lev < nlev; lev++)
    {
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


namespace MfixIO {

//
//
//
void
MfixRW::ReportGridStats () const
{
  BL_PROFILE("mfix::volEpsWgtSum()");

  std::vector<long> counts(6,0);

  int lev = 0;

  const MultiFab* volfrac =  &(ebfactory[lev]->getVolFrac());

  // Count the number of regular cells
  counts[0] = static_cast<int>(amrex::ReduceSum(*volfrac, *(m_leveldata[lev]->ep_g), 0,
    [=] AMREX_GPU_HOST_DEVICE (Box const & bx,
                               Array4<const Real> const & vfrc,
                               Array4<const Real> const & ep) -> int
    {
      int dm = 0;

      amrex::Loop(bx, [vfrc,ep,&dm] (int i, int j, int k) noexcept
      {if(vfrc(i,j,k)==1.0) dm += 1;});

      return dm;
    }));

  // Count the number of covered cells
  counts[1] = static_cast<int>(amrex::ReduceSum( *volfrac, *(m_leveldata[lev]->ep_g), 0,
    [=] AMREX_GPU_HOST_DEVICE (Box const & bx,
                               Array4<const Real> const & vfrc,
                               Array4<const Real> const & ep) -> int
    {
      int dm = 0;

      amrex::Loop(bx, [vfrc,ep,&dm] (int i, int j, int k) noexcept
      {if(vfrc(i,j,k)==0.0) dm += 1;});

      return dm;
    }));

  // Count the number of cut cells
  counts[2] = static_cast<int>(amrex::ReduceSum( *volfrac, *(m_leveldata[lev]->ep_g), 0,
    [=] AMREX_GPU_HOST_DEVICE (Box const & bx,
                               Array4<const Real> const & vfrc,
                               Array4<const Real> const & ep) -> int
    {
      int dm = 0;

      amrex::Loop(bx, [vfrc,ep,&dm] (int i, int j, int k) noexcept
      {if(0.0 < vfrc(i,j,k) && vfrc(i,j,k) < 1.0) dm += 1;});

      return dm;
    }));

  int regular(0), covered(0), cut(0);

#ifdef _OPENMP
#pragma omp parallel reduction(+:regular, covered, cut) if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(*m_leveldata[lev]->vel_g,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    const auto& vel_fab   =
      static_cast<EBFArrayBox const&>((*m_leveldata[lev]->vel_g)[mfi]);

    const Box& bx     = mfi.tilebox();
    const auto& flags = vel_fab.getEBCellFlagFab();

    // Count number of regular grids
    if (flags.getType(amrex::grow(bx,1)) == FabType::regular ) {
      regular += 1;
    } else if (flags.getType(amrex::grow(bx,1)) == FabType::covered ) {
      covered += 1;
    } else {
      cut += 1;
    }
  }

  counts[3] = regular;
  counts[4] = covered;
  counts[5] = cut;

  ParallelDescriptor::ReduceLongSum(counts.data(), 6);

  if(ParallelDescriptor::IOProcessor()){
    printf("\n\n****************************************\n");
    printf("  Coverage report:  Grids        Cells\n");
    printf("          regular:  %5ld   %10ld\n", counts[3], counts[0]);
    printf("          covered:  %5ld   %10ld\n", counts[4], counts[1]);
    printf("              cut:  %5ld   %10ld\n", counts[5], counts[2]);
    printf("****************************************\n\n");
  }

}

}


//
// Print the minimum volume fraction and cell location.
//
IntVect
mfix::mfix_print_min_epg ()
{

#ifndef AMREX_USE_GPU

  for (int lev = 0; lev <= finest_level; lev++) {

    const Real tolerance = std::numeric_limits<Real>::epsilon();
    auto& ld = *m_leveldata[lev];
    const Real min_epg = ld.ep_g->min(0);

    for (MFIter mfi(*ld.vel_g,false); mfi.isValid(); ++mfi) {
      Box const& bx = mfi.tilebox();
      Array4<Real const> const& epg = ld.ep_g->const_array(mfi);

      IntVect epg_cell = {-100,-100,-100};
      int found(0);

      amrex::ParallelFor(bx, [epg, min_epg, &found, &epg_cell, tolerance]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        if( amrex::Math::abs(epg(i,j,k) - min_epg) < tolerance ){
          epg_cell[0] = i;
          epg_cell[1] = j;
          epg_cell[2] = k;
          found +=1;
        }
      });

      if(found > 0){
        amrex::Print(Print::AllProcs)
          << std::endl << std::endl << "min epg "  << min_epg
          << "  at " << epg_cell[0] << "  " << epg_cell[1] << "  " << epg_cell[2]
          << "   total found " << found << std::endl << std::endl;

        return epg_cell;

      }

      //AMREX_ALWAYS_ASSERT(min_epg > 0.275);

    } // mfi
  } // lev
#endif
  IntVect fake = {0,0,0};
  return fake;

}


IntVect
mfix::mfix_locate_max_eps (Vector< MultiFab* >& ep_s_in, const Real max_eps)
{

  int imax(-100), jmax(-100), kmax(-100);

#ifdef AMREX_USE_GPU
  amrex::Gpu::DeviceScalar<int> imax_gpu(-100);
  amrex::Gpu::DeviceScalar<int> jmax_gpu(-100);
  amrex::Gpu::DeviceScalar<int> kmax_gpu(-100);

  int* p_imax = imax_gpu.dataPtr();
  int* p_jmax = jmax_gpu.dataPtr();
  int* p_kmax = kmax_gpu.dataPtr();
#endif

  for (int lev = 0; lev <= finest_level; lev++) {

    constexpr Real tolerance = std::numeric_limits<Real>::epsilon();

    for (MFIter mfi(*ep_s_in[lev],false); mfi.isValid(); ++mfi) {
      Box const& bx = mfi.tilebox();
      // Array4<Real const> const& epg = ld.ep_g->const_array(mfi);
      Array4<const Real> const& eps = ep_s_in[lev]->const_array(mfi);

#ifdef AMREX_USE_GPU
      amrex::ParallelFor(bx, [eps, max_eps, p_imax, p_jmax, p_kmax]
#else
      amrex::ParallelFor(bx, [eps, max_eps, &imax, &jmax, &kmax]
#endif
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if( amrex::Math::abs(eps(i,j,k) - max_eps) < tolerance ){
#ifdef AMREX_USE_GPU
          *p_imax = i;
          *p_jmax = j;
          *p_kmax = k;
#else
          imax = i;
          jmax = j;
          kmax = k;
#endif
        }
      });
    } // mfi
  } // lev


  Gpu::synchronize();

#ifdef AMREX_USE_GPU
  imax = imax_gpu.dataValue();
  jmax = jmax_gpu.dataValue();
  kmax = kmax_gpu.dataValue();
#endif

  ParallelDescriptor::ReduceIntMax({imax, jmax, kmax});

  return {imax, jmax, kmax};
}
