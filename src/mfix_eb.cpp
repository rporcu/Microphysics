#include <AMReX.H>
#include <AMReX_PlaneIF.H>
#include <AMReX_SphereIF.H>
#include <AMReX_IntersectionIF.H>

#include <mfix_level.H>
#include <mfix_F.H>

void
mfix_level::make_eb_geometry(int lev)
{
    if (geom[lev].isAllPeriodic()) return;

    Box domain(geom[lev].Domain());
    Real dx = geom[lev].CellSize()[0];

    GeometryShop* workshop;

    int exists;
    RealVect normal, center;
    PlaneIF* plane;
    Vector<BaseIF*> planes;
    planes.resize(0);

    for (int i = 1; i <= 500; i++) {
      mfix_get_walls(&i, &exists, &normal, &center);
      if(exists){
        std::cout << "Normal " << normal << std::endl;
        std::cout << "Center " << center << std::endl;
        plane = new PlaneIF(normal,center,true);
        planes.push_back(plane);
      }
    }
    IntersectionIF all_planes(planes);
    workshop = new GeometryShop(all_planes);

    // This part is generic once you have defined the workshop
    EBIndexSpace* ebis = AMReX_EBIS::instance();
    ebis->define(domain, RealVect::Zero, dx, *workshop);

    // set up ebfactory
    int m_eb_basic_grow_cells = 2;
    int m_eb_volume_grow_cells = 2;
    int m_eb_full_grow_cells = 2;
    EBSupport m_eb_support_level = EBSupport::full;

    EBTower::Build();

    ebfactory = std::unique_ptr<EBFArrayBoxFactory>(new EBFArrayBoxFactory(geom[lev], grids[lev], dmap[lev],
                 {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                  m_eb_support_level));
}
