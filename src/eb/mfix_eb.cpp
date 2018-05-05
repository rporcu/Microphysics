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
#include <eb_levelset.H>
#include <mfix_level.H>
#include <mfix_F.H>

void
mfix_level::make_eb_geometry(int lev)
{
    if (geom[lev].isAllPeriodic()) return;

    // Implicit functions for:
    //    * impfunc       -> all EBs in the domain
    //    * impfunc_poly2 -> EBs belonging to the (polynomial) walls
    //    * impfunc_walls -> EBs belonging to the (mfix.dat) flat walls
    std::unique_ptr<BaseIF> impfunc;
    std::unique_ptr<BaseIF> impfunc_poly2;
    std::unique_ptr<BaseIF> impfunc_walls;

    ParmParse pp("mfix");

    bool use_walls = true;
    bool use_poly2 = false;

    pp.query("use_walls", use_walls);
    pp.query("use_poly2", use_poly2);


   /****************************************************************************
    * Generate PolynomialIF representing the non-planar EB walls               *
    ****************************************************************************/

    if(use_poly2) {
      amrex::Print() << "Using poly2 geometry" << std::endl;
      impfunc_poly2 = get_poly(lev, SpaceDim, "poly2");


     /**************************************************************************
      * Generate BaseIF representing the planar EB walls                       *
      *          and intersect with with PolynomialIF                          *
      **************************************************************************/

      if(use_walls){
        bool has_walls;
        // Store wall implicit function separately (used by level-set)
        impfunc_walls = get_walls(lev, false, has_walls);

        Vector<BaseIF*> funcs(2);
        funcs[0] = impfunc_poly2.get();
        funcs[1] = impfunc_walls.get();
        IntersectionIF implicit(funcs);
        impfunc.reset(implicit.newImplicitFunction());

      } else {
        impfunc.reset(impfunc_poly2->newImplicitFunction());
      }


     /**************************************************************************
      * Generate BaseIF representing ONLY the planar EB walls                  *
      **************************************************************************/
    } else if(use_walls){
      bool has_walls;
      impfunc_walls = get_walls(lev, true, has_walls);
      impfunc.reset(impfunc_walls->newImplicitFunction());
    }



    int max_level = 0;
    int grid_size = 16;
    bool eb_verbosity = true;

   /****************************************************************************
    *                                                                          *
    * Fill Level-set using:                                                    *
    *      -> Planes (where the GeometryShop's implicit function is a signed   *
    *         distance): implicit function's value                             *
    *      -> Poly2 (where GeometryShop's implicit function is singed but not  *
    *         a distance): min distance to EB facets                           *
    * Note: this requires building and destroying the EBTower (twice), so any  *
    * EBTower data built before this will be lost...                           *
    *                                                                          *
    ****************************************************************************/

    Geometry geom_eb = LSUtility::make_eb_geometry(* level_set);
    if(use_walls){
        // Define components of the GeometryShop separately:
        GeometryShop gshop_walls(* impfunc_walls, eb_verbosity);

        // Define the EBIS first using only the walls...
        AMReX_EBIS::instance()->define(geom_eb.Domain(),
                                       RealVect::Zero,  // ......... origin of EBIndexSpace
                                       geom_eb.CellSize()[0],  // .. reference cell size of EBIndexSpace [1]
                                       gshop_walls,  // ............ GeometryShop object
                                       grid_size, max_level);
        // [1]: EBIndexSpace internally assumes an isotropic grid. Any
        // anisotropic implicit function (e.g AnisotrpicPlaneIF) uses dx as a
        // reference, and rescales dy and dz wrt dx. => dx goes here.

        EBTower::Build();
        // GeometryShop's Planes' implicit function is actually a signed distance function
        //      => it's just easier to fill the level-set this way
        level_set->intersection_ebis(* AMReX_EBIS::instance());
        EBTower::Destroy();
    }
    if(use_poly2){
        // Define components of the GeometryShop separately:
        GeometryShop gshop_poly2(* impfunc_poly2, eb_verbosity);

        // Define the EBIS using only the poly2 (after deleting the walls-only EBTower)...
        AMReX_EBIS::instance()->define(geom_eb.Domain(),
                                       RealVect::Zero,  // ......... origin of EBIndexSpace
                                       geom_eb.CellSize()[0],  // .. reference cell size of EBIndexSpace [1]
                                       gshop_poly2,  // ............ GeometryShop object
                                       grid_size, max_level);

        EBTower::Build();
        // GeometryShop's PolynomialIF is not a signed distance function...
        //      => it's easier to use PolynomialIF to build an
        //         EBFArrayBoxFactory which defines our EB surface now
        //          => define the level set as the (signed) distance to the
        //             closest point on the EB-facets
        int eb_pad = level_set->get_eb_pad();
        EBFArrayBoxFactory eb_factory_poly2(geom_eb, level_set->get_eb_ba(), dmap[lev],
                                            {eb_pad, eb_pad, eb_pad}, EBSupport::full);
        level_set->intersection_ebf(eb_factory_poly2, * AMReX_EBIS::instance());
        EBTower::Destroy();
    }

    // store copy of level set (for plotting).
    std::unique_ptr<MultiFab> ls_data = level_set->coarsen_data();
    ls[lev]->copy(* ls_data, 0, 0, 1, 0, 0);
    ls[lev]->FillBoundary(geom[lev].periodicity());

    /***************************************************************************
     *                                                                         *
     * Build standard EB Factories                                             *
     *                                                                         *
     ***************************************************************************/

    // dx : cell size used by EBIndexSpace
    Box domain(geom[lev].Domain());
    Real dx = geom[lev].CellSize()[0];

    GeometryShop gshop(*impfunc, eb_verbosity);
    AMReX_EBIS::instance()->define(domain, RealVect::Zero, dx, gshop, grid_size, max_level);

    // set up ebfactory
    int m_eb_basic_grow_cells = 2;
    int m_eb_volume_grow_cells = 2;
    int m_eb_full_grow_cells = 2;
    EBSupport m_eb_support_level = EBSupport::full;

    EBTower::Build();

    ebfactory          = std::unique_ptr<EBFArrayBoxFactory>(
                         new EBFArrayBoxFactory(geom[lev], grids[lev], dmap[lev],
                         {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                         m_eb_support_level));

    particle_ebfactory = std::unique_ptr<EBFArrayBoxFactory>(
                         new EBFArrayBoxFactory(geom[lev], grids[lev], dmap[lev],
                         {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                         m_eb_support_level));

    eb_normals         = pc->EBNormals(lev, particle_ebfactory.get(), dummy.get());
}



//std::unique_ptr<BaseIF>
//make_eb_wall(int dir, Real position, Real height, Real width, int lev, int water_tight) {
//    Vector<BaseIF *> planes;
//    planes.resize(0); // Make sure `planes` is empty (will use `push_back`)
//
//}



void
mfix_level::make_eb_hourglass(int lev)
{
    if (geom[lev].isAllPeriodic()) return;

    // Implicit functions for:
    //    * impfunc       -> all EBs in the domain
    //    * impfunc_poly2 -> EBs belonging to the (curved/bulb) walls
    //    * impfunc_walls -> EBs belonging to the (mfix.dat) flat walls
    std::unique_ptr<BaseIF> impfunc;
    std::unique_ptr<BaseIF> impfunc_unpolys;
    std::unique_ptr<BaseIF> impfunc_walls;

    amrex::Print() << "Using poly geometry" << std::endl;


   /****************************************************************************
    * Define components: poly1 and poly2 => Hourglass "bulbs"                  *
    *                    walls => Flat base intersecting with bottom bulb      *
    ****************************************************************************/

    // Construct the implicit function for the "top" bulb. Note the the curved
    // walls the defined using a PolynomialIF of order 5.
    std::unique_ptr<BaseIF> poly1t = get_poly(lev, 5, "poly1");
    // Construct the implicit function for the "bottom" bulb.
    std::unique_ptr<BaseIF> poly2t = get_poly(lev, 5, "poly2");

    // Construct any additional walls. The smooth hourglass has a flat base
    // which defined via the planar walls (in the mfix.dat).
    bool has_walls;
    impfunc_walls = get_walls(lev, false, has_walls);

    if (! has_walls)
        amrex::Print() << "WARNING: You have made an hourglass without end walls!"
            << std::endl
            << "If you're not careful, all the particles might fall out the bottom..."
            << std::endl;


   /****************************************************************************
    * Combine components: unbulbs => union of top and bottom PolynomialIF      *
    *                       +--- intersects  walls                             *
    ****************************************************************************/

    Vector<BaseIF *> bulbs{poly1t.get(),poly2t.get()};
    UnionIF unpolys(bulbs);

    impfunc_unpolys = std::unique_ptr<BaseIF>(unpolys.newImplicitFunction());

    Vector<BaseIF*> funcs{impfunc_unpolys.get(), impfunc_walls.get()};
    IntersectionIF implicit(funcs);

    impfunc.reset(implicit.newImplicitFunction());


    int max_level = 0;
    int grid_size = 16;
    bool eb_verbosity = true;



    /***************************************************************************
     *                                                                         *
     * Fill Level-set using:                                                   *
     *      -> Planes (where the GeometryShop's implicit function is a signed  *
     *         distance): implicit function's value                            *
     *      -> Poly (where GeometryShop's implicit function is singed but not  *
     *         a distance): min distance to EB facets                          *
     * Note: this requires building and destroying the EBTower (twice), so any *
     * EBTower data built before this will be lost...                          *
     *                                                                         *
     ***************************************************************************/

    // Define both components of the GeometryShop separately:
    GeometryShop gshop_upoly(* impfunc_unpolys, eb_verbosity);
    GeometryShop gshop_walls(* impfunc_walls, eb_verbosity);

    Geometry geom_eb = LSUtility::make_eb_geometry(* level_set);

    // Define the EBIS first using only the walls...
    AMReX_EBIS::instance()->define(geom_eb.Domain(),
                                   RealVect::Zero,  // ......... origin of EBIndexSpace
                                   geom_eb.CellSize()[0],  // .. reference cell size of EBIndexSpace [1]
                                   gshop_walls,  // ............ GeometryShop object
                                   grid_size, max_level);
    // [1]: EBIndexSpace internally assumes an isotropic grid. Any anisotropic
    // implicit function (e.g AnisotrpicPlaneIF) uses dx as a reference, and
    // rescales dy and dz wrt dx. => dx goes here.


    EBTower::Build();
    // GeometryShop's Planes' implicit function is actually a signed distance function
    //      => it's just easier to fill the level-set this way
    level_set->intersection_ebis(* AMReX_EBIS::instance());
    EBTower::Destroy();

    // Define the EBIS using only the poly (after deleting the walls-only EBTower)...
    AMReX_EBIS::instance()->define(geom_eb.Domain(),
                                   RealVect::Zero,  // ......... origin of EBIndexSpace
                                   geom_eb.CellSize()[0],  // .. reference cell size of EBIndexSpace [1, above]
                                   gshop_upoly,  // ............ GeometryShop object
                                   grid_size, max_level);

    EBTower::Build();
    // GeometryShop's PolynomialIF is not a signed distance function...
    //      => it's easier to use PolynomialIF to build an EBFArrayBoxFactory
    //         which defines our EB surface now
    //          => define the level set as the (signed) distance to the closest
    //             point on the EB-facets
    int eb_pad = level_set->get_eb_pad();
    EBFArrayBoxFactory eb_factory_poly(geom_eb, level_set->get_eb_ba(), dmap[lev],
                                       {eb_pad, eb_pad, eb_pad}, EBSupport::full);
    level_set->intersection_ebf(eb_factory_poly, * AMReX_EBIS::instance());
    EBTower::Destroy();

    // store copy of level set (for plotting).
    std::unique_ptr<MultiFab> ls_data = level_set->coarsen_data();
    ls[lev]->copy(* ls_data, 0, 0, 1, 0, 0 /*ls[lev]->nGrow(), ls[lev]->nGrow()*/);
    ls[lev]->FillBoundary(geom[lev].periodicity());



    /***************************************************************************
     *                                                                         *
     * Build standard EB Factories                                             *
     *                                                                         *
     ***************************************************************************/

    Box domain(geom[lev].Domain());
    Real dx = geom[lev].CellSize()[0];

    GeometryShop gshop(*impfunc, eb_verbosity);
    AMReX_EBIS::instance()->define(domain, RealVect::Zero, dx, gshop, grid_size, max_level);

    // set up ebfactory
    int m_eb_basic_grow_cells = 2;
    int m_eb_volume_grow_cells = 2;
    int m_eb_full_grow_cells = 2;
    EBSupport m_eb_support_level = EBSupport::full;

    EBTower::Build();

    ebfactory          = std::unique_ptr<EBFArrayBoxFactory>(
                         new EBFArrayBoxFactory(geom[lev], grids[lev], dmap[lev],
                         {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                         m_eb_support_level));

    particle_ebfactory = std::unique_ptr<EBFArrayBoxFactory>(
                         new EBFArrayBoxFactory(geom[lev], grids[lev], dmap[lev],
                         {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                         m_eb_support_level));

    eb_normals         = pc->EBNormals(lev, particle_ebfactory.get(), dummy.get());
}



void
mfix_level::make_eb_clr(int lev)
{
    std::unique_ptr<BaseIF> impfunc;

    ParmParse pp("mfix");

    amrex::Print() << "Creating CLR geometry\n";

    Vector<Real> transvec(3);
    RealVect translation;
    Real lradius, lheight;
    int cylinder_dir;

    // CLR can be constructed in either "water-tight" mode (each component is
    // assembled as a level-set intersection/union), or the much faster EB mode
    // (where the level-set is filled using the already assembled eb-factory).
    // Note that the EB mode can't have leaks, but it can have irregular edges.

    bool water_tight = false;
    pp.query("levelset__water-tight", water_tight);

    // If the mfix_level::make_cylinder uses a union
    //  => ensure that mfix_level::level_set is initialized to min
    if(water_tight) level_set->invert();

    //ct_ls_mf = 0;
    //std::unique_ptr<MultiFab> ls_mf = level_set->copy_data();
    //amrex::VisMF::Write(* ls_mf, "ls_empty");

    //------------------------------------------------------------- Riser
    pp.getarr("riser_translate", transvec,  0, 3);
    translation = RealVect(transvec);

    pp.query("riser_lower_radius", lradius);
    pp.query("riser_lower_height", lheight);

    cylinder_dir=1;
    std::unique_ptr<BaseIF> riser_lower = make_cylinder(cylinder_dir, lradius, lheight, translation,
                                                        lev, water_tight);


    pp.query("riser_upper_radius", lradius);
    pp.query("riser_upper_height", lheight);

    cylinder_dir=1;
    std::unique_ptr<BaseIF> riser_upper = make_cylinder(cylinder_dir, lradius, lheight, translation,
                                                        lev, water_tight);

    //------------------------------------------------------------- J-leg
    pp.getarr("jleg_htranslate", transvec,  0, 3);
    translation = RealVect(transvec);

    pp.query("jleg_horz_radius", lradius);
    pp.query("jleg_horz_height", lheight);

    cylinder_dir=0;
    std::unique_ptr<BaseIF> jleg_horz = make_cylinder(cylinder_dir, lradius, lheight, translation,
                                                      lev, water_tight);


    pp.getarr("jleg_vtranslate", transvec,  0, 3);
    translation = RealVect(transvec);

    pp.query("jleg_vert_radius", lradius);
    pp.query("jleg_vert_height", lheight);

    cylinder_dir=1;
    std::unique_ptr<BaseIF> jleg_vert = make_cylinder(cylinder_dir, lradius, lheight, translation,
                                                      lev, water_tight);

    //------------------------------------------------------------- reactor
    pp.getarr("reactor_translate", transvec,  0, 3);
    translation = RealVect(transvec);

    pp.query("reactor_radius", lradius);
    pp.query("reactor_height", lheight);

    cylinder_dir=1;
    std::unique_ptr<BaseIF> reactor = make_cylinder(cylinder_dir, lradius, lheight, translation,
                                                    lev, water_tight);

    //------------------------------------------------------------- loop-seal to reactor
    pp.getarr("ls2fbr_translate", transvec,  0, 3);
    translation = RealVect(transvec);

    pp.query("ls2fbr_radius", lradius);
    pp.query("ls2fbr_height", lheight);

    cylinder_dir=1;
    std::unique_ptr<BaseIF> ls2fbr = make_cylinder(cylinder_dir, lradius, lheight, translation,
                                                   lev, water_tight);


    //------------------------------------------------------------- loop-seal to reactor
    pp.getarr("loopseal_translate", transvec,  0, 3);
    translation = RealVect(transvec);

    pp.query("loopseal_radius", lradius);
    pp.query("loopseal_height", lheight);

    cylinder_dir=1;
    std::unique_ptr<BaseIF> loopseal = make_cylinder(cylinder_dir, lradius, lheight, translation,
                                                     lev, water_tight);



    //------------------------------------------------------------- loop-seal to reactor
    pp.getarr("cy2ls_translate", transvec,  0, 3);
    translation = RealVect(transvec);

    pp.query("cy2ls_radius", lradius);
    pp.query("cy2ls_height", lheight);

    cylinder_dir=1;
    std::unique_ptr<BaseIF> cy2ls = make_cylinder(cylinder_dir, lradius, lheight, translation,
                                                  lev, water_tight);


    //------------------------------------------------------------- cyclone
    pp.getarr("cyclone_translate", transvec,  0, 3);
    translation = RealVect(transvec);

    pp.query("cyclone_radius", lradius);
    pp.query("cyclone_height", lheight);

    cylinder_dir=1;
    std::unique_ptr<BaseIF> cyclone = make_cylinder(cylinder_dir, lradius, lheight, translation,
                                                    lev, water_tight);


    //------------------------------------------------------------- crossover
    pp.getarr("crossover_translate", transvec,  0, 3);
    translation = RealVect(transvec);

    pp.query("crossover_radius", lradius);
    pp.query("crossover_height", lheight);

    cylinder_dir=0;
    std::unique_ptr<BaseIF> crossover = make_cylinder(cylinder_dir, lradius, lheight, translation,
                                                      lev, water_tight);

    //--------------------------------------------------------------------

    Vector<BaseIF*> clr_parts(10);
    clr_parts[0] = riser_lower.get();
    clr_parts[1] = riser_upper.get();
    clr_parts[2] = jleg_horz.get();
    clr_parts[3] = jleg_vert.get();
    clr_parts[4] = reactor.get();
    clr_parts[5] = ls2fbr.get();
    clr_parts[6] = loopseal.get();
    clr_parts[7] = cy2ls.get();
    clr_parts[8] = cyclone.get();
    clr_parts[9] = crossover.get();

    UnionIF clr(clr_parts);
    impfunc.reset(clr.newImplicitFunction());



    int max_level = 0;
    int grid_size = 16;
    bool eb_verbosity = true;
    GeometryShop gshop(* impfunc, eb_verbosity);

   /****************************************************************************
    *                                                                          *
    * IF not using water-tight mode (i.e. level-set was not already filled     *
    * using the `make_cylinder` function):                                     *
    *      -> Construct an EBIndexSpace and using the refinement specified     *
    *         in `levelset__eb_refinement`.                                    *
    *      -> Fill level-set using an eb-factory defined on this EBIS.         *
    *      Note: This is much faster than "water-tight" mode, but does not     *
    *      resolve edges a well => *might have leaks*                          *
    *                                                                          *
    ****************************************************************************/

    if(! water_tight){

        // Construct EBIS geometry
        Geometry geom_eb = LSUtility::make_eb_geometry(* level_set);
        AMReX_EBIS::instance()->define(geom_eb.Domain(),
                                       RealVect::Zero,  // ......... origin of EBIndexSpace
                                       geom_eb.CellSize()[0],  // .. reference cell size of EBIndexSpace [1]
                                       gshop,  // .................. GeometryShop object
                                       grid_size, max_level);
        // [1]: EBIndexSpace internally assumes an isotropic grid. Any
        // anisotropic implicit function (e.g AnisotrpicPlaneIF) uses dx as a
        // reference, and rescales dy and dz wrt dx. => dx goes here.

        // Temporary EB tower
        EBTower::Build();

        // Fill level_set.
        int eb_grow = level_set->get_eb_pad();
        EBFArrayBoxFactory eb_factory(geom_eb, level_set->get_eb_ba(), dmap[lev],
                                      {eb_grow, eb_grow, eb_grow}, EBSupport::full);
        level_set->intersection_ebf(eb_factory, * AMReX_EBIS::instance());

        // Destroy temporary EB tower (rebuilt below)
        EBTower::Destroy();
    }



   /****************************************************************************
    *                                                                          *
    * Build standard EB Factories                                              *
    *                                                                          *
    ****************************************************************************/

    Box domain(geom[lev].Domain());
    Real dx = geom[lev].CellSize()[0];

    AMReX_EBIS::instance()->define(domain, RealVect::Zero, dx, gshop, grid_size, max_level);

    // set up ebfactory
    int m_eb_basic_grow_cells = 2;
    int m_eb_volume_grow_cells = 2;
    int m_eb_full_grow_cells = 2;
    EBSupport m_eb_support_level = EBSupport::full;

    EBTower::Build();

    ebfactory          = std::unique_ptr<EBFArrayBoxFactory>(
                         new EBFArrayBoxFactory(geom[lev], grids[lev], dmap[lev],
                         {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                         m_eb_support_level));

    particle_ebfactory = std::unique_ptr<EBFArrayBoxFactory>(
                         new EBFArrayBoxFactory(geom[lev], grids[lev], dmap[lev],
                         {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                         m_eb_support_level));

    eb_normals         = pc->EBNormals(lev, particle_ebfactory.get(), dummy.get());

    // Promote completed copy of level set into the mfix_level.
    std::unique_ptr<MultiFab> ls_data = level_set->coarsen_data();
    ls[lev]->copy(* ls_data, 0, 0, 1, 0, 0 /*ls[lev]->nGrow(), ls[lev]->nGrow()*/);
    ls[lev]->FillBoundary(geom[lev].periodicity());
}


void
mfix_level::make_eb_clr_riser(int lev)
{
    std::unique_ptr<BaseIF> impfunc;

    ParmParse pp("mfix");

    amrex::Print() << "Creating CLR riser geometry\n";

    Vector<Real> transvec(3);
    RealVect translation;
    int cylinder_dir;

    // CLR can be constructed in either "water-tight" mode (each component is
    // assembled as a level-set intersection/union), or the much faster EB mode
    // (where the level-set is filled using the already assembled eb-factory).
    // Note that the EB mode can't have leaks, but it can have irregular edges.

    bool water_tight = false;
    pp.query("levelset__water-tight", water_tight);

    // If the mfix_level::make_cylinder uses a union
    //  => ensure that mfix_level::level_set is initialized to min
    if(water_tight) level_set->invert();


    //------------------------------------------------------------- Riser
    pp.getarr("riser_translate", transvec,  0, 3);
    translation = RealVect(transvec);

    // cylinder direction (y-axis)
    cylinder_dir=1;

    Real lradius, uradius;
    Real lheight, uheight, cheight;
    pp.query("riser_lower_radius", lradius);
    pp.query("riser_lower_height", lheight);
    pp.query("riser_upper_radius", uradius);
    pp.query("riser_upper_height", uheight);
    pp.query("riser_c2c_height",   cheight);

    // Build lower cylinder
    std::unique_ptr<BaseIF> riser_lower =
      make_cylinder(cylinder_dir, lradius, lheight, translation, lev, water_tight);

    // Offset for lower cylinder height.
    translation[cylinder_dir] += lheight;

    // Build cylinder-to-cylinder connector.
    std::unique_ptr<BaseIF> riser_c2c =
      make_cone(cylinder_dir, uradius, lradius, cheight, translation, lev, water_tight);

    // Offset for connector height.
    translation[cylinder_dir] += cheight;

    // Build upper cylinder, offset for lower cylinder height
    std::unique_ptr<BaseIF> riser_upper =
      make_cylinder(cylinder_dir, uradius, uheight, translation, lev, water_tight);

    Vector<BaseIF*> clr_parts(3);
    clr_parts[0] = riser_lower.get();
    clr_parts[1] = riser_c2c.get();
    clr_parts[2] = riser_upper.get();

    UnionIF clr(clr_parts);
    impfunc.reset(clr.newImplicitFunction());


    int max_level = 0;
    int grid_size = 16;
    bool eb_verbosity = true;
    GeometryShop gshop(* impfunc, eb_verbosity);

   /****************************************************************************
    *                                                                          *
    * IF not using water-tight mode (i.e. level-set was not already filled     *
    * using the `make_cylinder` function):                                     *
    *      -> Construct an EBIndexSpace and using the refinement specified     *
    *         in `levelset__eb_refinement`.                                    *
    *      -> Fill level-set using an eb-factory defined on this EBIS.         *
    *      Note: This is much faster than "water-tight" mode, but does not     *
    *      resolve edges a well => *might have leaks*                          *
    *                                                                          *
    ****************************************************************************/

    if(! water_tight){

        // Construct EBIS geometry
        Geometry geom_eb = LSUtility::make_eb_geometry(* level_set);
        AMReX_EBIS::instance()->define(geom_eb.Domain(),
                                       RealVect::Zero,  // ......... origin of EBIndexSpace
                                       geom_eb.CellSize()[0],  // .. reference cell size of EBIndexSpace [1]
                                       gshop,  // .................. GeometryShop object
                                       grid_size, max_level);

        // [1]: EBIndexSpace internally assumes an isotropic grid. Any
        // anisotropic implicit function (e.g AnisotrpicPlaneIF) uses dx as a
        // reference, and rescales dy and dz wrt dx. => dx goes here.

        // Temporary EB tower
        EBTower::Build();

        // Fill level_set.
        int eb_grow = level_set->get_eb_pad();
        EBFArrayBoxFactory eb_factory(geom_eb, level_set->get_eb_ba(), dmap[lev],
                                      {eb_grow, eb_grow, eb_grow}, EBSupport::full);
        level_set->intersection_ebf(eb_factory, * AMReX_EBIS::instance());

        // Destroy temporary EB tower (rebuilt below)
        EBTower::Destroy();
    }



   /****************************************************************************
    *                                                                          *
    * Build standard EB Factories                                              *
    *                                                                          *
    ****************************************************************************/

    Box domain(geom[lev].Domain());
    Real dx = geom[lev].CellSize()[0];

    AMReX_EBIS::instance()->define(domain, RealVect::Zero, dx, gshop, grid_size, max_level);

    // set up ebfactory
    int m_eb_basic_grow_cells = 2;
    int m_eb_volume_grow_cells = 2;
    int m_eb_full_grow_cells = 2;
    EBSupport m_eb_support_level = EBSupport::full;

    EBTower::Build();

    ebfactory          = std::unique_ptr<EBFArrayBoxFactory>(
                         new EBFArrayBoxFactory(geom[lev], grids[lev], dmap[lev],
                         {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                         m_eb_support_level));

    particle_ebfactory = std::unique_ptr<EBFArrayBoxFactory>(
                         new EBFArrayBoxFactory(geom[lev], grids[lev], dmap[lev],
                         {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                         m_eb_support_level));

    eb_normals         = pc->EBNormals(lev, particle_ebfactory.get(), dummy.get());

    // Promote completed copy of level set into the mfix_level.
    std::unique_ptr<MultiFab> ls_data = level_set->coarsen_data();
    ls[lev]->copy(* ls_data, 0, 0, 1, 0, 0 /*ls[lev]->nGrow(), ls[lev]->nGrow()*/);
    ls[lev]->FillBoundary(geom[lev].periodicity());
}


std::unique_ptr<BaseIF> mfix_level::get_walls(int lev, bool anisotropic, bool & has_walls) {
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


std::unique_ptr<BaseIF> mfix_level::get_poly(int lev, int max_order, std::string field_prefix) {
    // Construct the ParamParse database field names based on the
    // `field_prefix` string:
    ParmParse pp("mfix");

    // Coefficients vector is stored in the inputs database with the field name:
    //      <field_prefix>_[x,y,z]_coeffs
    const std::array<const string, 3> var_names{"x", "y", "z"};
    std::array<string, 3> field_names;
    for(int i = 0; i < 3; i++) {
        std::stringstream field_name;
        field_name << field_prefix;
        field_name << "_" << var_names[i] << "_coeffs";
        field_names[i] = field_name.str();
    }

    // There are two more fields assoicated with the PolynomialIF:
    //      <field_prefix>_mirror    (true if fluid is inside the PolynomialIF)
    //      <field_prefix>_translate (vector representing center-axis position)
    std::stringstream mirror_field, translate_field;
    mirror_field << field_prefix << "_mirror";
    translate_field << field_prefix << "_translate";


    // Generate vector representing polynomial
    Vector<PolyTerm> poly;
    for(int idir = 0; idir < 3; idir++) {
        Vector<Real> coefvec(SpaceDim);

        if(idir == 0)      pp.getarr(field_names[idir].c_str(), coefvec, 0, max_order);
        else if(idir == 1) pp.getarr(field_names[idir].c_str(), coefvec, 0, max_order);
        else if(idir == 2) pp.getarr(field_names[idir].c_str(), coefvec, 0, max_order);

        for(int lc = 0; lc < max_order; lc++) {
            // x^(lc) term
            Real coef = coefvec[lc];
            IntVect powers = IntVect::Zero;
            powers[idir] = lc;

            PolyTerm mono = {.coef = coef, .powers = powers};
            poly.push_back(mono);
        }
    }


    /***************************************************************************
     * Construct PolynomialIF (and apply translation)                          *
     ***************************************************************************/

    bool flip = false;
    pp.query(mirror_field.str().c_str(), flip);

    PolynomialIF mirror(poly, flip);
    RealVect translation;

    Vector<Real> transvec(SpaceDim);
    pp.getarr(translate_field.str().c_str(), transvec, 0, SpaceDim);

    for(int idir = 0; idir < 3; idir++)
        translation[idir] = transvec[idir];

    TransformIF poly2(mirror);
    poly2.translate(translation);

    return std::unique_ptr<BaseIF>(poly2.newImplicitFunction());
}


std::unique_ptr<BaseIF> mfix_level::make_cylinder(int dir, Real radius, Real length, const RealVect & translation,
                                                  int lev, bool water_tight) {
    // Construct a cylinder implicit function with finite radius and axis
    // offset (translation). The cylinder can be oriented along any cartesian
    // axis (dir).
    std::unique_ptr<BaseIF> cylinder_IF;


    // Polynomial defining (curved) cylinder walls parallel to a given axis:
    //     IF = a^2 + b^2 - R^2
    // where a, b \in {a, x, z} - {axis} for example, if the cylinder lies on
    // the y-axis => IF = x^2 + z^2 - R^2
    Vector<PolyTerm> poly;
    for(int idir = 0; idir < 3; idir++) {
        // Constucts the coefficient vector describing a cylinder with
        // orientation given by axis `dir`:
        //    *  coefvec[0] = R^2 term
        //    *  coefvec[2] = {x,y,z}^2 term
        Vector<Real> coefvec(3);
        if( idir == dir ) coefvec = { - std::pow(radius, 2), 0. ,0.};
        else              coefvec = {0., 0., 1};

        for(int lc = 0; lc < 3; lc++) {
            // x^(lc) term
            IntVect powers = IntVect::Zero;
            powers[idir] = lc;
            PolyTerm mono = {.coef = coefvec[lc], .powers = powers};
            poly.push_back(mono);
        }
    }


    // Internal flow cylinder
    PolynomialIF cylinder0(poly, true);
    TransformIF cylinder1(cylinder0);
    cylinder1.translate(translation);


    // box to clip to correct length
    RealVect normal = RealVect::Zero , center = RealVect::Zero;
    Vector<std::unique_ptr<BaseIF>> planes;

    center[dir] = 0.0;
    normal[dir] = 1.0;
    planes.push_back(std::unique_ptr<BaseIF>(new PlaneIF(normal, center, true)));

    center[dir] = length;
    normal[dir] =-1.0;
    planes.push_back(std::unique_ptr<BaseIF>(new PlaneIF(normal, center, true)));

    // The IntersectionIF constructor requires a vector of pointers to the
    Vector<BaseIF *> plane_ptrs = {planes[0].get(), planes[1].get()};
    IntersectionIF bounding_box(plane_ptrs);
    TransformIF walls(bounding_box);
    walls.translate(translation);

    Vector<BaseIF * > funcs(2);
    funcs[0] = & cylinder0;
    funcs[1] = & bounding_box;

    IntersectionIF cylinder(funcs);

    TransformIF cylinder_trans(cylinder);
    cylinder_trans.translate(translation);

    cylinder_IF.reset(cylinder_trans.newImplicitFunction());


    // IF we are not using the water-tight mode, return now:

    if(! water_tight)
        return cylinder_IF;

    // ELSE construct the level-set by unioning each cylinder (intersected
    // component) of the CLR using the level-set
    //   => corners are much more cleanly resolved, but is much slower.



    int max_level = 0;
    int grid_size = 16;
    bool eb_verbosity = true;

   /****************************************************************************
    *                                                                          *
    * Fill Level-set using:                                                    *
    *      -> Walls (where the GeometryShop's implicit function is a signed    *
    *         distance): implicit function's value                             *
    *      -> Cylinder (where GeometryShop's implicit function is singed but   *
    *         not a distance): min distance to EB facets                       *
    *      Note: this requires building and destroying the EBTower (twice),    *
    *      so any EBTower data built before this will be lost...               *
    *                                                                          *
    ****************************************************************************/

    // Define both components of the GeometryShop separately:
    GeometryShop gshop_upoly(cylinder1, eb_verbosity);
    GeometryShop gshop_walls(walls, eb_verbosity);

    // Define temporary level-sets used for constructing the cylinder:
    LSFactory ls_cylinder(* level_set);

    Geometry geom_eb = LSUtility::make_eb_geometry(ls_cylinder);

    // Define the EBIS first using only the walls...
    AMReX_EBIS::instance()->define(geom_eb.Domain(),
                                   RealVect::Zero,  // ......... origin of EBIndexSpace
                                   geom_eb.CellSize()[0],  // .. reference cell size of EBIndexSpace [1]
                                   gshop_walls,  // ............ GeometryShop object
                                   grid_size, max_level);
    // [1]: EBIndexSpace internally assumes an isotropic grid. Any anisotropic
    // implicit function (e.g AnisotrpicPlaneIF) uses dx as a reference, and
    // rescales dy and dz wrt dx. => dx goes here.

    EBTower::Build();
    // GeometryShop's Planes' implicit function is actually a signed distance function
    //      => it's just easier to fill the level-set this way
    ls_cylinder.intersection_ebis(* AMReX_EBIS::instance());
    EBTower::Destroy();

    // Define the EBIS using only the poly (after deleting the walls-only EBTower)...
    AMReX_EBIS::instance()->define(geom_eb.Domain(),
                                   RealVect::Zero,  // ......... origin of EBIndexSpace
                                   geom_eb.CellSize()[0],  // .. reference cell size of EBIndexSpace [1, above]
                                   gshop_upoly,  // ............ GeometryShop object
                                   grid_size, max_level);

    EBTower::Build();
    // GeometryShop's PolynomialIF is not a signed distance function...
    //      => it's easier to use PolynomialIF to build an EBFArrayBoxFactory
    //         which defines our EB surface
    //          => define the level set as the (signed) distance to the closest
    //             point on the EB-facets
    int eb_grow = level_set->get_eb_pad();
    EBFArrayBoxFactory eb_factory_poly(geom_eb, level_set->get_eb_ba(), dmap[lev],
                                       {eb_grow, eb_grow, eb_grow}, EBSupport::full);

    // Only EB facets that are "in range" (within `n_pad` of the local
    // BoxArray) are considered for filling the EB level-set. flag_valid = {1,0}
    // indicates if a cell's BoxArray contained any EB facets (1), or if the
    // cell's value is invalid (0) because all EB facets where too far away in
    // order to be considered.
    std::unique_ptr<iMultiFab> flag_valid = ls_cylinder.intersection_ebf(eb_factory_poly,
                                                                         * AMReX_EBIS::instance());
    EBTower::Destroy();

    // Union temporary level set (ls_cylinder) with the main level set.
    level_set->update_union(* ls_cylinder.get_data(), * flag_valid);

    return cylinder_IF;
}


std::unique_ptr<BaseIF>
mfix_level::make_cone(int dir, Real radius1, Real radius2, Real height,
                      const RealVect & translation, int lev, bool water_tight)
{
    std::unique_ptr<BaseIF> cone_IF;


    // The cone is defined using 2 radii. Set rad[1,2] such that rad2 >= rad1
    bool rad2_larger = radius2 > radius1;
    Real rad1 = rad2_larger ? radius1 : radius2;
    Real rad2 = rad2_larger ? radius2 : radius1;
    Real H    = height*rad2 / (rad2 - rad1);
    Real k    = rad2 / H;


    // Cone polynomial parallel to a given axis:
    //     IF = a^2 + b^2 - (k c)^2
    // where a, b \in {a, x, z} - {axis} and c is the coordinate describing the
    // axis. For example, if the cone lies on the y-axis:
    //     IF = x^2 + z^2 - (k y)^2
    Vector<PolyTerm> poly;
    for(int idir = 0; idir < 3; idir++) {
        // Constucts the coefficient vector describing a con with orientation
        // given by axis `dir`:
        //    *  coefvec[2] = k^2 term
        //    *  coefvec[2] = {x,y,z}^2 term
        Vector<Real> coefvec(3);
        if( idir == dir ) coefvec = {0., 0., - std::pow(k, 2)};
        else              coefvec = {0., 0., 1.};

        for(int lc = 0; lc < 3; lc++) {
            // x^(lc) term
            IntVect powers = IntVect::Zero;
            powers[idir] = lc;
            PolyTerm mono = {.coef = coefvec[lc], .powers = powers};
            poly.push_back(mono);
        }
    }


    // Internal flow cone
    PolynomialIF cone0(poly, true);
    TransformIF cone1(cone0);
    cone1.translate(translation);


    // box to clip to correct length
    Real cmin = rad2_larger ? -H          : H - height;
    Real cmax = rad2_larger ? -H + height : H;

    RealVect normal = RealVect::Zero, center = RealVect::Zero;
    Vector<std::unique_ptr<BaseIF>> planes;

    center[dir] = cmin;
    normal[dir] = 1.0;
    planes.push_back(std::unique_ptr<BaseIF>(new PlaneIF(normal, center, true)));

    center[dir] = cmax;
    normal[dir] =-1.0;
    planes.push_back( std::unique_ptr<BaseIF>(new PlaneIF(normal, center, true)));

    // The IntersectionIF constructor requires a vector of pointers to the
    Vector<BaseIF *> plane_ptrs = {planes[0].get(), planes[1].get()};
    IntersectionIF bounding_box(plane_ptrs);
    TransformIF walls(bounding_box);
    walls.translate(translation);

    Vector<BaseIF*> funcs(2);
    funcs[0] = & cone0;
    funcs[1] = & bounding_box;

    IntersectionIF cone(funcs);

    // The shift needs to take into account offset from the origin
    RealVect shift = RealVect(translation);
    shift[dir] = shift[dir]-cmin;

    TransformIF cone_trans(cone);
    cone_trans.translate(shift);

    cone_IF.reset(cone_trans.newImplicitFunction());

    // IF we are not using the water-tight mode, return now:

    if(! water_tight)
        return cone_IF;

    // ELSE construct the level-set by unioning each cone (intersected
    // component) of the CLR using the level-set
    //   => corners are much more cleanly resolved, but is much slower.



    int max_level = 0;
    int grid_size = 16;
    bool eb_verbosity = true;


   /****************************************************************************
    *                                                                          *
    * Fill Level-set using:                                                    *
    *      -> Walls (where the GeometryShop's implicit function is a signed    *
    *         distance): implicit function's value                             *
    *      -> Cone (where GeometryShop's implicit function is singed but       *
    *         not a distance): min distance to EB facets                       *
    *      Note: this requires building and destroying the EBTower (twice),    *
    *      so any EBTower data built before this will be lost...               *
    *                                                                          *
    ****************************************************************************/

    // Define both components of the GeometryShop separately:
    GeometryShop gshop_upoly(cone1, eb_verbosity);
    GeometryShop gshop_walls(walls, eb_verbosity);

    // Define temporary level-sets used for constructing the cone:
    LSFactory ls_cone(* level_set);

    // Define the EBIS first using only the walls...
    Geometry geom_eb = LSUtility::make_eb_geometry(ls_cone);
    AMReX_EBIS::instance()->define(geom_eb.Domain(),
                                   RealVect::Zero,  // ......... origin of EBIndexSpace
                                   geom_eb.CellSize()[0],  // .. reference cell size of EBIndexSpace [1]
                                   gshop_walls,  // ............ GeometryShop object
                                   grid_size, max_level);
    // [1]: EBIndexSpace internally assumes an isotropic grid. Any anisotropic
    // implicit function (e.g AnisotrpicPlaneIF) uses dx as a reference, and
    // rescales dy and dz wrt dx. => dx goes here.

    EBTower::Build();
    // GeometryShop's Planes' implicit function is actually a signed distance function
    //      => it's just easier to fill the level-set this way
    ls_cone.intersection_ebis(* AMReX_EBIS::instance());
    EBTower::Destroy();

    // Define the EBIS using only the poly (after deleting the walls-only EBTower)...
    AMReX_EBIS::instance()->define(geom_eb.Domain(),
                                   RealVect::Zero,  // ......... origin of EBIndexSpace
                                   geom_eb.CellSize()[0],  // .. reference cell size of EBIndexSpace [1, above]
                                   gshop_upoly,  // ............ GeometryShop object
                                   grid_size, max_level);

    EBTower::Build();
    // GeometryShop's PolynomialIF is not a signed distance function...
    //      => it's easier to use PolynomialIF to build an EBFArrayBoxFactory
    //         which defines our EB surface
    //          => define the level set as the (signed) distance to the closest
    //             point on the EB-facets
    int eb_grow = level_set->get_eb_pad();
    EBFArrayBoxFactory eb_factory_poly(geom_eb, level_set->get_eb_ba(), dmap[lev],
                                       {eb_grow, eb_grow, eb_grow}, EBSupport::full);

    // Only EB facets that are "in range" (within `n_pad` of the local
    // BoxArray) are considered for filling the EB level-set. flag_valid =
    // {1,0} indicates if a cell's BoxArray contained any EB facets (1), or if
    // the cell's value is invalid (0) because all EB facets where too far away
    // in order to be considered.
    std::unique_ptr<iMultiFab> flag_valid = ls_cone.intersection_ebf(eb_factory_poly, * AMReX_EBIS::instance());
    EBTower::Destroy();

    // Union temporary level set (ls_cylinder) with the main level set.
    level_set->update_union(* ls_cone.get_data(), * flag_valid);

    return cone_IF;
}
