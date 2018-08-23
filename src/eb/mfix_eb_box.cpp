#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Intersection.H>

#include <AMReX_EB_levelset.H>
#include <mfix_level.H>
#include <mfix_eb_F.H>


/****************************************************************************
 * Function to create a simple rectangular box with EB walls.               *
 *                                                                          *
 ****************************************************************************/
void
mfix_level::make_eb_box(int lev)
{
    ParmParse pp("box");

    int max_level_here = 0;

    /***************************************************************************
     *                                                                         *
     * Build standard EB Factories                                             *
     *                                                                         *
     ***************************************************************************/

    // set up ebfactory
    int m_eb_basic_grow_cells = nghost;
    int m_eb_volume_grow_cells = nghost;
    int m_eb_full_grow_cells = nghost;
    EBSupport m_eb_support_level = EBSupport::full;

    EB2::useEB2(true);

    int max_coarsening_level = 100;
    Real offset    = 1.0e-8;

    amrex::Print() << " " << std::endl;
    amrex::Print() << "Now making the ebfactory's ..." << std::endl;

    if (geom[lev].isAllPeriodic() )
    {
        make_eb_regular(lev);
    } 
    else 
    {
        Real xlo = geom[lev].ProbLo(0) + offset;
        Real xhi = geom[lev].ProbHi(0) - offset;

        if (geom[lev].isPeriodic(0))
        {
            xlo = 2.0*geom[lev].ProbLo(0) - geom[lev].ProbHi(0);
            xhi = 2.0*geom[lev].ProbHi(0) - geom[lev].ProbLo(0);
        }

        Real ylo = geom[lev].ProbLo(1) + offset;
        Real yhi = geom[lev].ProbHi(1) - offset;

        if (geom[lev].isPeriodic(1))
        {
            ylo = 2.0*geom[lev].ProbLo(1) - geom[lev].ProbHi(1);
            yhi = 2.0*geom[lev].ProbHi(1) - geom[lev].ProbLo(1);
        }

        Real zlo = geom[lev].ProbLo(2) + offset;
        Real zhi = geom[lev].ProbHi(2) - offset;

        if (geom[lev].isPeriodic(2))
        {
            zlo = 2.0*geom[lev].ProbLo(2) - geom[lev].ProbHi(2);
            zhi = 2.0*geom[lev].ProbHi(2) - geom[lev].ProbLo(2);
        }

        Array<Real,3>  point_lox{ xlo, 0.0, 0.0};
        Array<Real,3> normal_lox{-1.0, 0.0, 0.0};
        Array<Real,3>  point_hix{ xhi, 0.0, 0.0};
        Array<Real,3> normal_hix{ 1.0, 0.0, 0.0};

        Array<Real,3>  point_loy{0.0, ylo, 0.0};
        Array<Real,3> normal_loy{0.0,-1.0, 0.0};
        Array<Real,3>  point_hiy{0.0, yhi, 0.0};
        Array<Real,3> normal_hiy{0.0, 1.0, 0.0};

        Array<Real,3>  point_loz{0.0, 0.0, zlo};
        Array<Real,3> normal_loz{0.0, 0.0,-1.0};
        Array<Real,3>  point_hiz{0.0, 0.0, zhi};
        Array<Real,3> normal_hiz{0.0, 0.0, 1.0};

        EB2::PlaneIF plane_lox(point_lox,normal_lox);
        EB2::PlaneIF plane_hix(point_hix,normal_hix);

        EB2::PlaneIF plane_loy(point_loy,normal_loy);
        EB2::PlaneIF plane_hiy(point_hiy,normal_hiy);

        EB2::PlaneIF plane_loz(point_loz,normal_loz);
        EB2::PlaneIF plane_hiz(point_hiz,normal_hiz);

        auto gshop = EB2::makeShop(EB2::makeUnion(plane_lox,plane_hix,
                                                  plane_loy,plane_hiy,
                                                  plane_loz,plane_hiz ));
 
        EB2::Build(gshop, geom.back(), max_level_here,
                   max_level_here+max_coarsening_level);

        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geom[lev]);

        if (solve_fluid)
           ebfactory[lev].reset(new EBFArrayBoxFactory(eb_level,
                    geom[lev], grids[lev], dmap[lev],
                    {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                     m_eb_full_grow_cells}, m_eb_support_level));

        if (solve_dem)
           particle_ebfactory.reset(new EBFArrayBoxFactory(eb_level,
                    geom[lev], grids[lev], dmap[lev],
                    {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                     m_eb_full_grow_cells}, m_eb_support_level));

       amrex::Print() << "Done making the ebfactory's ..." << std::endl;
       amrex::Print() << " " << std::endl;
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
