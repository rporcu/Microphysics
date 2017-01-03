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
MyParticleContainer::InitData()
{
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

