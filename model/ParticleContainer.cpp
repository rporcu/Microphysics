#include <limits>

#include <ParticleContainer.H>
#include <ParticleIterator.H>

int     MyParticleContainer::do_tiling = 0;
IntVect MyParticleContainer::tile_size   { D_DECL(1024000,8,8) };

MyParticleContainer::MyParticleContainer (AmrCore* amr_core)
    : ParticleContainer<PIdx::nattribs,0,std::vector<Particle<PIdx::nattribs,0> > >
      (amr_core->GetParGDB())
{
    ReadStaticParameters();

    this->SetVerbose(0);

    m_particles.reserve(m_gdb->maxLevel()+1);
    m_particles.resize (m_gdb->finestLevel()+1);

    m_partdata.reserve(m_gdb->maxLevel()+1);
    m_partdata.resize (m_gdb->finestLevel()+1);
}

void
MyParticleContainer::ReadStaticParameters ()
{
    static bool initialized = false;
    if (!initialized)
    {
	ParmParse pp("particles");

	pp.query("do_tiling",  do_tiling);

	Array<int> ts(BL_SPACEDIM);
	if (pp.queryarr("tile_size", ts)) {
	    tile_size = IntVect(ts);
	}

	initialized = true;
    }
}

void
MyParticleContainer::AllocData ()
{
    m_particles.resize(m_gdb->finestLevel()+1);
    m_partdata.resize(m_gdb->finestLevel()+1);

    for (int lev = 0; lev <= m_gdb->finestLevel(); ++ lev)
    {
	auto& partleveldata = m_partdata[lev];
	const BoxArray& ba = m_gdb->ParticleBoxArray(lev);
	const DistributionMapping& dm = m_gdb->ParticleDistributionMap(lev);

	MultiFab foo(ba, 1, 0, dm, Fab_noallocate);
	for (MFIter mfi(foo); mfi.isValid(); ++mfi)
	{
	    int i = mfi.index();
	    partleveldata[i] = Array<std::unique_ptr<Array<Real> > > (PIdx::npartdata);
	    for (auto& d : partleveldata[i])
	    {
		d.reset(new Array<Real>());
	    }
	}
    }
}

void
MyParticleContainer::Evolve (int lev, Real dt)
{
    BL_PROFILE("MyPC::Evolve()");

    const Geometry& gm  = m_gdb->Geom(lev);

#if (BL_SPACEDIM == 3)
    const Real* dx = gm.CellSize();
#elif (BL_SPACEDIM == 2)
    Real dx[3] = { gm.CellSize(0), std::numeric_limits<Real>::quiet_NaN(), gm.CellSize(1) };
#endif
}

