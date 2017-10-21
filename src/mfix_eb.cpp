#include <AMReX.H>
#include <AMReX_GeometryShop.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_SphereIF.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBTower.H>
#include <AMReX_EBFArrayBox.H>

#include <mfix_level.H>

void
mfix_level::make_eb_geometry(int lev)
{
    Box domain(geom[lev].Domain());
    Real dx = geom[lev].CellSize()[0];

    Real size = 1.0;

    // set up our sphere
    Real radius = 0.25*size;
    RealVect center = 0.5*size*RealVect::Unit;
    bool insideRegular = true;
    SphereIF sphere(radius, center, insideRegular);
    GeometryShop workshop(sphere);
    EBIndexSpace* ebis = AMReX_EBIS::instance();
    ebis->define(domain, RealVect::Zero, dx, workshop);

    // set up ebfactory
    int m_eb_basic_grow_cells = 5;
    int m_eb_volume_grow_cells = 4;
    int m_eb_full_grow_cells = 2;
    EBSupport m_eb_support_level = EBSupport::full;

    EBTower::Build();

    EBFArrayBoxFactory ebfactory(geom[lev], grids[lev], dmap[lev],
                                 {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                                 m_eb_support_level);

    std::array<const MultiCutFab*, AMREX_SPACEDIM> areafrac;
    const MultiFab* volfrac;
    const MultiCutFab* bndrycent;

    areafrac  =  ebfactory.getAreaFrac();
    volfrac   = &ebfactory.getVolFrac();
    bndrycent = &ebfactory.getBndryCent();
}
