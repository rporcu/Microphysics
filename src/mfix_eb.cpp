#include <AMReX.H>
#include <AMReX_PlaneIF.H>
#include <AMReX_SphereIF.H>

#include <mfix_level.H>

void
mfix_level::make_eb_geometry(int lev)
{
    bool make_sphere = false;
    bool make_planes = true;

    if (make_sphere && make_planes)
       amrex::Abort("Need to choose between planes and sphere");

    Box domain(geom[lev].Domain());
    Real dx = geom[lev].CellSize()[0];

    Real size = 1.0;

    GeometryShop* workshop;  

    // Set up our sphere
    if (make_sphere) 
    {
       Real radius = 0.25*size;
       RealVect center = 0.5*size*RealVect::Unit;
       bool insideRegular = true;
       SphereIF sphere(radius, center, insideRegular);
       workshop = new GeometryShop(sphere);
    }

    // Set up a plane
    if (make_planes) 
    {
       RealVect normal = RealVect(0,0,1);
       RealVect center = RealVect(0,0,0.2);
       PlaneIF plane(normal,center,true);
       workshop = new GeometryShop(plane);
    }

    // This part is generic once you have defined the workshop
    EBIndexSpace* ebis = AMReX_EBIS::instance();
    ebis->define(domain, RealVect::Zero, dx, *workshop);

    // set up ebfactory
    int m_eb_basic_grow_cells = 5;
    int m_eb_volume_grow_cells = 4;
    int m_eb_full_grow_cells = 2;
    EBSupport m_eb_support_level = EBSupport::full;

    EBTower::Build();

    ebfactory = new EBFArrayBoxFactory(geom[lev], grids[lev], dmap[lev],
                 {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                  m_eb_support_level);
}
