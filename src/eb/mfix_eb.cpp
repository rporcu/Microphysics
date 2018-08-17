#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

#include <AMReX_GeometryShop.H>
#include <AMReX_SphereIF.H>
#include <AMReX_PlaneIF.H>
#include <AMReX_AllRegularService.H>
#include <AMReX_FlatPlateGeom.H>
#include <AMReX_EBISLayout.H>
#include <AMReX_EBGraph.H>
#include <AMReX_EBDebugOut.H>
#include <AMReX_EBCellFAB.H>
#include <AMReX_EBCellFactory.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_UnionIF.H>
#include <AMReX_TransformIF.H>
#include <AMReX_ComplementIF.H>
#include <AMReX_IntersectionIF.H>
#include <AMReX_LatheIF.H>
#include <AMReX_PolynomialIF.H>
#include <AMReX_AnisotropicDxPlaneIF.H>
#include <AMReX_AnisotropicIF.H>

//#include <AMReX_VisMF.H>  // amrex::VisMF::Write(MultiFab)
//#include <sstream>

#include <algorithm>
#include <AMReX_EB_levelset.H>
#include <mfix_level.H>
#include <mfix_eb_F.H>

void
mfix_level::make_eb_geometry (int lev)
{
#if 0
    if (lev == 0)
    {
       bool hourglass    = false;
       bool clr          = false;
       bool clr_riser    = false;

       ParmParse pp("mfix");
       pp.query("hourglass", hourglass);
       pp.query("clr", clr);
       pp.query("clr_riser", clr_riser);

       if (hourglass) {
           make_eb_hourglass(lev);
       } else if(clr) {
           make_eb_clr(lev);
       } else if(clr_riser) {
           make_eb_clr_riser(lev);
       } else {
           make_eb_cylinder(lev);
       }
    }
#endif

    make_eb_cylinder(lev);
}

void
mfix_level::make_eb_cylinder(int lev)
{
    ParmParse pp("mfix");

    bool use_walls   = true;
    bool use_poly2   = false;
    bool use_divider = false;

    pp.query("use_walls",   use_walls);
    pp.query("use_poly2",   use_poly2);
    pp.query("use_divider", use_divider);


   /****************************************************************************
    * IF using divider => Extract parameters                                   *
    ****************************************************************************/

    int div_dir = 0;
    Real div_pos, div_height, div_width;

    if(use_divider) {
        pp.query("divider_position", div_pos);
        pp.query("divider_height",   div_height);
        pp.query("divider_width",    div_width);
        pp.query("div_dir",          div_dir);
    }

   /****************************************************************************
    *                                                                          *
    *  Set some of the parameters to be used in making the EB factories        *
    *                                                                          *
    ****************************************************************************/
       int max_level_here = 0;

   /****************************************************************************
    *                                                                          *
    * Fill level-set:                                                          *
    *                                                                          *
    ****************************************************************************/

#if 0
    if (solve_dem)
    {
   
       amrex::Print() << "Creating the levelset ..." << std::endl;

       bool eb_verbosity = true;
       int grid_size = 16;

       fill_levelset( lev, use_walls, use_poly2,
                      * impfunc_walls_part.get(), * impfunc_poly2.get(),
                      max_level_here, grid_size, eb_verbosity            );

       // store copy of level set (for plotting).
       std::unique_ptr<MultiFab> ls_data = level_set->coarsen_data();
       ls[lev]->copy(* ls_data, 0, 0, 1, 0, 0);
       ls[lev]->FillBoundary(geom[lev].periodicity());

       amrex::Print() << "Done making the levelset ..." << std::endl;
    }
#endif

    /***************************************************************************
     *                                                                         *
     * Build standard EB Factories                                             *
     *                                                                         *
     ***************************************************************************/

    // dx : cell size used by EBIndexSpace
    Box domain(geom[lev].Domain());

    // set up ebfactory
    int m_eb_basic_grow_cells = nghost;
    int m_eb_volume_grow_cells = nghost;
    int m_eb_full_grow_cells = nghost;
    EBSupport m_eb_support_level = EBSupport::full;

    EB2::useEB2(true);

   /****************************************************************************
    *                                                                          *
    * THIS FILLS PARTICLE EBFACTORY
    *                                                                          *
   /****************************************************************************/

    if (solve_dem)
    {
       amrex::Print() << " " << std::endl;
       amrex::Print() << "Now  making the particle ebfactory ..." << std::endl;

       bool inside = true;
       Real radius = 0.0002;
       Real height = 0.008;
       int direction = 0;
       Array<Real,3> center{0.002,0.0005,0.0005};
       EB2::CylinderIF my_cyl(radius, height, direction, center, inside);

       Array<Real,3> point{0.000,0.0005,0.0005};
       Array<Real,3> normal{1.0,0.0,0.0};
       EB2::PlaneIF my_plane(point,normal);

       auto gshop = EB2::makeShop(EB2::makeUnion(my_cyl,my_plane));
       int max_coarsening_level = 100;
       EB2::Build(gshop, geom.back(), max_level_here, max_level_here+max_coarsening_level);
       const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
       const EB2::Level& eb_level = eb_is.getLevel(geom[lev]);

       particle_ebfactory.reset(new EBFArrayBoxFactory(eb_level, 
                   geom[lev], grids[lev], dmap[lev],
                   {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                   m_eb_support_level
               )
            );

//     eb_normals = pc->EBNormals(lev, particle_ebfactory.get(), dummy.get());

       amrex::Print() << "Done making the particle ebfactory ..." << std::endl;
       amrex::Print() << " " << std::endl;
    }

   /****************************************************************************
    *                                                                          *
    * THIS FILLS FLUID EBFACTORY
    *                                                                          *
   /****************************************************************************/

    if (solve_fluid)
    {
       amrex::Print() << "Now  making the fluid ebfactory ..." << std::endl;

       bool inside = true;
       Real radius = 0.0002;
       Real height = 0.008;
       int direction = 0;
       Array<Real,3> center{0.002,0.0005,0.0005};
       EB2::CylinderIF my_cyl(radius, height, direction, center, inside);

       auto gshop = EB2::makeShop(my_cyl);
       int max_coarsening_level = 100;
       EB2::Build(gshop, geom.back(), max_level_here, max_level_here+max_coarsening_level);

       const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
       const EB2::Level& eb_level = eb_is.getLevel(geom[lev]);

       ebfactory[lev].reset(new EBFArrayBoxFactory(eb_level, 
                   geom[lev], grids[lev], dmap[lev],
                   {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                   m_eb_support_level
               )
            );

       amrex::Print() << "Done making the fluid ebfactory ..." << std::endl;
    }

}

std::unique_ptr<BaseIF> mfix_level::get_walls(int lev, bool anisotropic, bool & has_walls) 
{
    // Extracts all walls from the mfix.dat

    has_walls = false;  // will be set to true if there are any walls

    // Walls can be defined per phase => Itterarte over all phases and check
    // each for walls in the mfix.dat
    Vector<std::unique_ptr<BaseIF>> planes;
    for (int i = 1; i <= 500; i++) {
        int exists;
        RealVect normal, center;
        mfix_get_walls(& i, & exists, & normal, & center);
        if(exists) {
            has_walls = true;
            amrex::Print() << "Normal " << normal << std::endl;
            amrex::Print() << "Center " << center << std::endl;
            std::unique_ptr<PlaneIF> plane;
            if(anisotropic) {
                RealVect dxVec;
                for(int idir = 0; idir < 3; idir++)
                    dxVec[idir] = geom[lev].CellSize()[idir];
                plane = std::unique_ptr<PlaneIF>(new AnisotropicDxPlaneIF(normal, center, true, dxVec));
            } else {
                plane = std::unique_ptr<PlaneIF>(new PlaneIF(normal, center, true));
            }
            planes.push_back(std::move(plane));
        }
    }

    // The IntersectionIF constructor requires a vector of pointers to the
    // individual PlaneIFs => Construct a pointer-vector. Note that these
    // pointers become invalid once the `planes` vector goes out of scope.
    Vector<BaseIF *> plane_ptrs;
    for(std::unique_ptr<BaseIF> & pl : planes)
        plane_ptrs.push_back(pl.get());
    IntersectionIF all_planes(plane_ptrs);

    return std::unique_ptr<BaseIF>(all_planes.newImplicitFunction());
}

std::unique_ptr<BaseIF> mfix_level::get_real_walls(int lev, bool anisotropic, bool & has_real_walls) 
{
    // Extracts all walls from the mfix.dat

    has_real_walls = false;  // will be set to true if there are any walls

    // Walls can be defined per phase => Itterarte over all phases and check
    // each for walls in the mfix.dat
    Vector<std::unique_ptr<BaseIF>> planes;
    for (int i = 1; i <= 500; i++) {
        int exists;
        RealVect normal, center;
        mfix_get_real_walls(& i, & exists, & normal, & center);
        if(exists) {
            has_real_walls = true;
            amrex::Print() << "Normal " << normal << std::endl;
            amrex::Print() << "Center " << center << std::endl;
            std::unique_ptr<PlaneIF> plane;
            if(anisotropic) {
                RealVect dxVec;
                for(int idir = 0; idir < 3; idir++)
                    dxVec[idir] = geom[lev].CellSize()[idir];
                plane = std::unique_ptr<PlaneIF>(new AnisotropicDxPlaneIF(normal, center, true, dxVec));
            } else {
                plane = std::unique_ptr<PlaneIF>(new PlaneIF(normal, center, true));
            }
            planes.push_back(std::move(plane));
        }
    }

    // The IntersectionIF constructor requires a vector of pointers to the
    // individual PlaneIFs => Construct a pointer-vector. Note that these
    // pointers become invalid once the `planes` vector goes out of scope.
    Vector<BaseIF *> plane_ptrs;
    for(std::unique_ptr<BaseIF> & pl : planes)
        plane_ptrs.push_back(pl.get());
    IntersectionIF all_planes(plane_ptrs);

    return std::unique_ptr<BaseIF>(all_planes.newImplicitFunction());
}
