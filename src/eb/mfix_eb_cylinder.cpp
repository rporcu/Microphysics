#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Intersection.H>

//#include <AMReX_VisMF.H>  // amrex::VisMF::Write(MultiFab)
//#include <sstream>

#include <algorithm>
#include <AMReX_EB_levelset.H>
#include <mfix_level.H>
#include <mfix_eb_F.H>


/****************************************************************************
 * Function to create a simple cylinder EB.                                 *
 *                                                                          *
 * Comments: The cylinder can either be closed or open at the bottom. At    *
 * this time, there is no option to put a "lid" on the top of the cylinder. *
 *                                                                          *
 * The bottom of the cylinder is only applied to particles. This allows for *
 * the fluid to not have an EB surface generated that would block gas flow  *
 * from a mass inflow boundary. A more general implementation would check   *
 * the mfix.dat file for inflows and outflows and correctly cap the top and *
 * bottom of the cylinder.                                                  *
 *                                                                          *
 ****************************************************************************/
void
mfix_level::make_eb_cylinder(int lev)
{
    ParmParse pp("cylinder");

    int max_level_here = 0;

    /****************************************************************************
     * Get cylinder information from inputs file.                               *
     ****************************************************************************/
    bool inside       = true;
    bool close_bottom = true;
    Real offset    = 1.0e-8;

    Real radius    = 0.0002;
    Real height    = 0.0080;

    int direction  = 0;
    Vector<Real> centervec(3);

    pp.query("internal_flow", inside);
    pp.query("closed_bottom", close_bottom);
    pp.query("bottom_offset", offset);

    pp.query("radius",     radius);
    pp.query("height",     height);
    pp.query("direction", direction);
    pp.getarr("center",   centervec,  0, 3);
    Array<Real,3> center={centervec[0],centervec[1],centervec[2]};

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

    amrex::Print() << " " << std::endl;
    amrex::Print() << " Internal Flow: " << inside << std::endl;
    amrex::Print() << " Radius:    " << radius    << std::endl;
    amrex::Print() << " Height:    " << height    << std::endl;
    amrex::Print() << " Offset:    " << offset    << std::endl;
    amrex::Print() << " Direction: " << direction << std::endl;
    amrex::Print() << " Center:    " << center[0] << ", "
                                     << center[1] << ", "
                                     << center[2] << std::endl;

    // Create the cylinder -- used for both fluid and particles
    EB2::CylinderIF my_cyl(radius, height, direction, center, inside);


   /****************************************************************************
    *                                                                          *
    * THIS FILLS PARTICLE EBFACTORY
    *                                                                          *
    ****************************************************************************/
    if (solve_dem)
    {
       amrex::Print() << " " << std::endl;
       amrex::Print() << "Now  making the particle ebfactory ..." << std::endl;

       if(close_bottom) {

         Array<Real,3> point{0.0, 0.0, 0.0};
         Array<Real,3> normal{0.0, 0.0, 0.0};

         point[direction] = geom[lev].ProbLo(direction) + offset;
         normal[direction] = -1.0;

         amrex::Print() << "Capping bottom: " << std::endl;
         amrex::Print() << "   Point:  " << point[0]  << ", "
                                         << point[1]  << ", "
                                         << point[2]  << std::endl;

         amrex::Print() << "   Normal: " << normal[0] << ", "
                                         << normal[1] << ", "
                                         << normal[2] << std::endl;

         EB2::PlaneIF my_plane(point,normal);

         auto gshop = EB2::makeShop(EB2::makeUnion(my_cyl,my_plane));
         int max_coarsening_level = 100;
         EB2::Build(gshop, geom.back(), max_level_here,
                    max_level_here+max_coarsening_level);
       }
       else {
         auto gshop = EB2::makeShop(my_cyl);
         int max_coarsening_level = 100;
         EB2::Build(gshop, geom.back(), max_level_here,
                    max_level_here+max_coarsening_level);
       }

       const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
       const EB2::Level& eb_level = eb_is.getLevel(geom[lev]);

       particle_ebfactory.reset(new EBFArrayBoxFactory(eb_level,
                   geom[lev], grids[lev], dmap[lev],
                   {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                    m_eb_full_grow_cells}, m_eb_support_level)
            );

//     eb_normals = pc->EBNormals(lev, particle_ebfactory.get(), dummy.get());

       amrex::Print() << "Done making the particle ebfactory ..." << std::endl;
       amrex::Print() << " " << std::endl;
    }

   /****************************************************************************
    *                                                                          *
    * THIS FILLS FLUID EBFACTORY
    *                                                                          *
    ****************************************************************************/

    if (solve_fluid)
    {
       amrex::Print() << "Now  making the fluid ebfactory ..." << std::endl;

       auto gshop = EB2::makeShop(my_cyl);
       int max_coarsening_level = 100;
       EB2::Build(gshop, geom.back(), max_level_here,
                  max_level_here+max_coarsening_level);

       const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
       const EB2::Level& eb_level = eb_is.getLevel(geom[lev]);

       ebfactory[lev].reset(new EBFArrayBoxFactory(eb_level,
                   geom[lev], grids[lev], dmap[lev],
                   {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                    m_eb_full_grow_cells}, m_eb_support_level)
            );

       amrex::Print() << "Done making the fluid ebfactory ..." << std::endl;
    }



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
}
