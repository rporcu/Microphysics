#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

//#include <AMReX_VisMF.H>  // amrex::VisMF::Write(MultiFab)
//#include <sstream>

#include <algorithm>
#include <AMReX_EB_LSCore.H>
#include <AMReX_EB_levelset.H>
#include <mfix.H>
#include <mfix_eb_F.H>


/********************************************************************************
 *                                                                              *
 * Function to create a simple cylinder EB.                                     *
 *                                                                              *
 * Comments: The cylinder can either be closed or open at the bottom. At this   *
 * time, there is no option to put a "lid" on the top of the cylinder.          *
 *                                                                              *
 * The bottom of the cylinder is only applied to particles. This allows for the *
 * fluid to not have an EB surface generated that would block gas flow from a   *
 * mass inflow boundary. A more general implementation would check the mfix.dat *
 * file for inflows and outflows and correctly cap the top and bottom of the    *
 * cylinder.                                                                    *
 *                                                                              *
 ********************************************************************************/
void mfix::make_eb_cylinder()
{
    ParmParse pp("cylinder");

    int max_level_here = 0;

    /****************************************************************************
     * Get cylinder information from inputs file.                               *
     ****************************************************************************/
    bool inside       = true;
    bool close_bottom = true;
    Real offset       = 1.0e-8;

    Real radius    = 0.0002;
    Real height    = -1.;

    int direction  = 0;
    Vector<Real> centervec(3);

    pp.query("internal_flow", inside);
    pp.query("closed_bottom", close_bottom);
    pp.query("bottom_offset", offset);

    pp.query("radius",     radius);
    pp.query("height",     height);
    pp.query("direction",  direction);
    pp.getarr("center",    centervec,  0, 3);
    Array<Real,3> center={centervec[0], centervec[1], centervec[2]};

    /****************************************************************************
     *                                                                          *
     * Build standard EB Factories                                              *
     *                                                                          *
     ****************************************************************************/

    // set up ebfactory
    EBSupport m_eb_support_level = EBSupport::full;

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
    amrex::Print() << "Building the cylinder (side wall) geometry ..." << std::endl;

    // Build the Cylinder geometry first representing the curved walls (this is
    // always present regardless of user input).
    EB2::CylinderIF my_cyl(radius, height, direction, center, inside);

    auto gshop_cyl = EB2::makeShop(my_cyl);
    int max_coarsening_level = 100;
    EB2::Build(gshop_cyl, geom.back(), max_level_here,
               max_level_here + max_coarsening_level);

    const EB2::IndexSpace & ebis_cyl = EB2::IndexSpace::top();

    amrex::Print() << "Done building the cylinder geometry" << std::endl;


    /****************************************************************************
    *                                                                           *
    * THIS FILLS PARTICLE EBFACTORY                                             *
    * NOTE: the "bottom" plane is only applied to particles => fluid EBFactory  *
    *       only needs the cylinder walls.                                      *
    *                                                                           *
    *****************************************************************************/

    if (solve_dem)
    {
       amrex::Print() << " " << std::endl;
       amrex::Print() << "Now making the particle ebfactory ..." << std::endl;

       // Plane IF generator for bottom plane (defined iff close_bottom==true)
       std::unique_ptr<GShopLSFactory<EB2::PlaneIF>> wall_lsfactory;

       if (close_bottom)
       {
           Array<Real,3> point{0.0, 0.0, 0.0};
           Array<Real,3> normal{0.0, 0.0, 0.0};

           point[direction] = geom[0].ProbLo(direction) + offset;
           normal[direction] = -1.0;

           amrex::Print() << "Capping bottom: " << std::endl;
           amrex::Print() << "   Point:  " << point[0]  << ", "
                                           << point[1]  << ", "
                                           << point[2]  << std::endl;

           amrex::Print() << "   Normal: " << normal[0] << ", "
                                           << normal[1] << ", "
                                           << normal[2] << std::endl;

           EB2::PlaneIF my_plane(point, normal);

           auto gshop = EB2::makeShop(EB2::makeUnion(my_cyl, my_plane));
           EB2::Build(gshop, geom.back(), max_level_here,
                      max_level_here + max_coarsening_level);

           EB2::GeometryShop<EB2::PlaneIF> gshop_wall = EB2::makeShop(my_plane);
           wall_lsfactory.reset(new GShopLSFactory<EB2::PlaneIF>(gshop_wall, * level_set));
       }

       const EB2::IndexSpace & eb_is = EB2::IndexSpace::top();

       for (int lev = 0; lev < nlev; lev++)
       {
          // Intercept the cylinder + plane level (if it is built)
          eb_level_particles = & eb_is.getLevel(geom[lev]);

          particle_ebfactory[lev].reset(new EBFArrayBoxFactory(* eb_level_particles,
                                                               geom[lev], grids[lev], dmap[lev],
                                                               {m_eb_basic_grow_cells,
                                                                m_eb_volume_grow_cells,
                                                                m_eb_full_grow_cells},
                                                               m_eb_support_level));

          eb_normals = pc->EBNormals(lev, particle_ebfactory[lev].get(), dummy.get());
       }

       /*************************************************************************
        *                                                                       *
        * Fill level-set:                                                       *
        *                                                                       *
        *************************************************************************/

       if (!levelset__restart)
       {
          for (int lev = 0; lev < nlev; lev++)
          {
           amrex::Print() << "Creating the levelset ..." << std::endl;

           // If there is a bottom plane, fill level set with plane IF first
           if(close_bottom) {
               std::unique_ptr<MultiFab> mf_impfunc_wall = wall_lsfactory->fill_impfunc();
               level_set->intersection_impfunc(* mf_impfunc_wall);
           }

           // Construct EB2 Index space based on the refined geometry
           // (level_set->get_eb_geom()). The IndexSpace's geometry needs to
           // match the one used by the eb_factory later.

           auto gshop_cyl = EB2::makeShop(my_cyl);
           int max_coarsening_level = 100;
           EB2::Build(gshop_cyl, level_set->get_eb_geom(), max_level_here,
                      max_level_here + max_coarsening_level);

           const EB2::IndexSpace & ebis_cyl    = EB2::IndexSpace::top();
           const EB2::Level & ebis_lev_cyl_ref = ebis_cyl.getLevel(level_set->get_eb_geom());

           GShopLSFactory<EB2::CylinderIF> cyl_lsfactory(gshop_cyl, * level_set);
           std::unique_ptr<MultiFab> mf_impfunc_cyl = cyl_lsfactory.fill_impfunc();

           // Construct EBFABFactory based on the refined EB geometry (built above).
           int eb_grow = level_set->get_eb_pad();
           EBFArrayBoxFactory eb_factory_cylinder(ebis_lev_cyl_ref,
                                                  level_set->get_eb_geom(),
                                                  level_set->get_eb_ba(),
                                                  level_set->get_dm(),
                                                  {eb_grow, eb_grow, eb_grow}, EBSupport::full);

           level_set->intersection_ebf(eb_factory_cylinder, * mf_impfunc_cyl );

           amrex::Print() << "Done making the levelset ..." << std::endl;
          }

       } else {
           amrex::Print() << "Loaded level-set is fine => skipping levelset calculation."
                          << std::endl;
       }

       amrex::Print() << "Done making the particle ebfactories ..." << std::endl;
       amrex::Print() << " " << std::endl;
    }

    /****************************************************************************
    *                                                                           *
    * THIS FILLS FLUID EBFACTORY                                                *
    *                                                                           *
    *****************************************************************************/

    if (solve_fluid)
    {
       amrex::Print() << "Now  making the fluid ebfactory ..." << std::endl;

       for (int lev = 0; lev < nlev; lev++)
       {
          const EB2::Level & ebis_lev_cyl  = ebis_cyl.getLevel(geom[lev]);
          eb_level_fluid = & ebis_lev_cyl;

          ebfactory[lev].reset(new EBFArrayBoxFactory(
                      * eb_level_fluid,
                      geom[lev], grids[lev], dmap[lev],
                      {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                       m_eb_full_grow_cells}, m_eb_support_level));

         amrex::Print() << "Done making the fluid ebfactory ..." << std::endl;
       }
    }
}


void mfix::make_amr_cylinder()
{
    ParmParse pp("cylinder");

    int max_level_here = 0;

    /****************************************************************************
     * Get cylinder information from inputs file.                               *
     ****************************************************************************/

    bool inside       = true;
    bool close_bottom = true;
    Real offset       = 1.0e-8;

    Real radius    = 0.0002;
    Real height    = -1.;

    int direction  = 0;
    Vector<Real> centervec(3);

    pp.query("internal_flow", inside);
    pp.query("closed_bottom", close_bottom);
    pp.query("bottom_offset", offset);

    pp.query("radius",     radius);
    pp.query("height",     height);
    pp.query("direction",  direction);
    pp.getarr("center",    centervec,  0, 3);
    Array<Real,3> center={centervec[0], centervec[1], centervec[2]};


    amrex::Print() << " " << std::endl;
    amrex::Print() << " Internal Flow: " << inside << std::endl;
    amrex::Print() << " Radius:    " << radius    << std::endl;
    amrex::Print() << " Height:    " << height    << std::endl;
    amrex::Print() << " Offset:    " << offset    << std::endl;
    amrex::Print() << " Direction: " << direction << std::endl;
    amrex::Print() << " Center:    " << center[0] << ", "
                                     << center[1] << ", "
                                     << center[2] << std::endl;


    /****************************************************************************
     *                                                                          *
     * Construct AMR level-set                                                  *
     *                                                                          *
     ***************************************************************************/

    if (use_amr_ls)
    {
        int lev_lowest = 0;

        const RealBox & rb = geom[lev_lowest].ProbDomain();
        Box domain = geom[lev_lowest].Domain();
        domain.coarsen(amr_ls_crse);

        const IntVect & dom_lo = domain.smallEnd();
        const IntVect & dom_hi = domain.bigEnd();

        // Picket-fence principle
        IntVect n_cells = dom_hi - dom_lo + IntVect{1, 1, 1};

        Vector<int> v_cells = {
            AMREX_D_DECL(n_cells[0], n_cells[1], n_cells[2])
        };

        amrex::Print() << "Declaring AMR levelset:" << std::endl
                       << "coarsest level: " << domain << " n_cells: " << n_cells;

        amrex::Print() << std::endl << rb << std::endl;

        EB2::CylinderIF if_cyl(radius, height, direction, center, inside);

        if (close_bottom)
        {
            Array<Real,3> point{0.0, 0.0, 0.0};
            Array<Real,3> normal{0.0, 0.0, 0.0};

            point[direction] = geom[lev_lowest].ProbLo(direction) + offset;
            normal[direction] = -1.0;

            EB2::PlaneIF if_plane(point, normal);
            auto if_union = EB2::makeUnion(if_cyl, if_plane);

            amr_level_set.reset(new LSCore<decltype(if_union)>(EB2::makeShop(if_union),
                                                               & rb, amr_ls_max_level, v_cells));
        }
        else
        {
            amr_level_set.reset(new LSCore<decltype(if_cyl)>(EB2::makeShop(if_cyl)));
        }

        amrex::Print() << "... done declaring AMR levelset" << std::endl;
    }
}
