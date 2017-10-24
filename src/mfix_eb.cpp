#include <AMReX.H>
#include <AMReX_PlaneIF.H>
#include <AMReX_SphereIF.H>
#include <AMReX_IntersectionIF.H>

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

#if 0
    RealVect normal = RealVect(0,1,0);
    RealVect center = RealVect(0,1e-15,0);
    PlaneIF plane(normal,center,true);
    workshop = new GeometryShop(plane);
#else

    // Set up a plane
    if (make_planes) 
    {
       RealVect normal, center;
       PlaneIF* plane;
       Vector<BaseIF*> planes;
       planes.resize(0);

       if (!geom[lev].isPeriodic(0))
       { 
          normal = RealVect(0,0,1);
          center = RealVect(0,0,1e-15);
          plane = new PlaneIF(normal,center,true);
          planes.push_back(plane);

          normal = RealVect(0,0,-1);
          center = RealVect(0,0,0.9999999999999);
          plane = new PlaneIF(normal,center,true);
          planes.push_back(plane);
       } 

       if (!geom[lev].isPeriodic(1))
       { 
          normal = RealVect(0,1,0);
          center = RealVect(0,1e-15,0);
          plane = new PlaneIF(normal,center,true);
          planes.push_back(plane);

          normal = RealVect(0,-1,0);
          center = RealVect(0,0.9999999999999,0);
          plane = new PlaneIF(normal,center,true);
          planes.push_back(plane);
       } 

       if (!geom[lev].isPeriodic(2))
       { 
          normal = RealVect(1,0,0);
          center = RealVect(1e-15,0,0);
          plane = new PlaneIF(normal,center,true);
          planes.push_back(plane);

          normal = RealVect(-1,0,0);
          center = RealVect(0.9999999999999,0,0);
          plane = new PlaneIF(normal,center,true);
          planes.push_back(plane);
       } 

       IntersectionIF all_planes(planes);
       workshop = new GeometryShop(all_planes);
    }
#endif

    // This part is generic once you have defined the workshop
    EBIndexSpace* ebis = AMReX_EBIS::instance();
    ebis->define(domain, RealVect::Zero, dx, *workshop);

    // set up ebfactory
    int m_eb_basic_grow_cells = 2;
    int m_eb_volume_grow_cells = 2;
    int m_eb_full_grow_cells = 2;
    EBSupport m_eb_support_level = EBSupport::full;

    EBTower::Build();

    ebfactory = new EBFArrayBoxFactory(geom[lev], grids[lev], dmap[lev],
                 {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                  m_eb_support_level);
}
