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

#include <AMReX_VisMF.H>
#include <sstream>

#include <eb_levelset.H>
#include <mfix_level.H>
#include <mfix_F.H>

void
mfix_level::make_eb_geometry(int lev)
{
    if (geom[lev].isAllPeriodic()) return;

    Box domain(geom[lev].Domain());
    Real dx = geom[lev].CellSize()[0];

    int exists;
    RealVect normal, center;
    PlaneIF* plane;
    Vector<BaseIF*> planes;
    planes.resize(0);

    std::unique_ptr<BaseIF> impfunc;
    std::unique_ptr<BaseIF> impfunc_poly2;
    std::unique_ptr<BaseIF> impfunc_walls;

    ParmParse pp("mfix");

    bool use_walls = true;
    bool use_poly2 = false;

    pp.query("use_walls", use_walls);
    pp.query("use_poly2", use_poly2);

    if(use_poly2){

      amrex::Print() << "Using poly2 geometry\n";

      Vector<PolyTerm> poly;

      PolyTerm mono;
      Real coef;
      IntVect powers;

      Vector<Real> coefvec(SpaceDim);
      Vector<int>  powersvec(SpaceDim);
      Vector<Real> transvec(SpaceDim);

      for(int idir = 0; idir < 3; idir++) {
        if( idir == 0) {
          pp.getarr("poly2_x_coeffs",  coefvec,   0, SpaceDim);
        } else if( idir == 1) {
          pp.getarr("poly2_y_coeffs",  coefvec,   0, SpaceDim);
        } else if( idir == 2) {
          pp.getarr("poly2_z_coeffs",  coefvec,   0, SpaceDim);
        }

        for(int lc = 0; lc < 3; lc++) {
          // x^(lc) term
          coef = coefvec[lc];
          powers = IntVect::Zero;
          powers[idir] = lc;

          mono.coef   = coef;
          mono.powers = powers;

          poly.push_back(mono);

        }
      }

      bool flip = false;
      pp.query("poly2_mirror", flip);

      PolynomialIF mirror(poly,flip);
      RealVect translation;

      pp.getarr("poly2_translate", transvec,  0, SpaceDim);

      for(int idir = 0; idir < 3; idir++)
        translation[idir] = transvec[idir];

      TransformIF poly2(mirror);
      poly2.translate(translation);

      impfunc_poly2 = std::unique_ptr<BaseIF>(poly2.newImplicitFunction());

      if(use_walls){ // Combine poly2 with walls
        for (int i = 1; i <= 500; i++) {
          mfix_get_walls(&i, &exists, &normal, &center);
          if(exists){
            //center[0] = 1.e-3;
            amrex::Print() << "Normal " << normal << std::endl;
            amrex::Print() << "Center " << center << std::endl;
            plane = new PlaneIF(normal,center,true);
            planes.push_back(plane);
          }
        }
        IntersectionIF all_planes(planes);

        impfunc_walls = std::unique_ptr<BaseIF>(all_planes.newImplicitFunction());

        Vector<BaseIF*> funcs(2);
        funcs[0] = &poly2;
        funcs[1] = &all_planes;
        IntersectionIF implicit(funcs);
        impfunc.reset(implicit.newImplicitFunction());

      } else {
        impfunc.reset(poly2.newImplicitFunction());
      }

    } else if(use_walls){ // Just walls

      RealVect dxVec;
      for(int idir = 0; idir < 3; idir++)
        dxVec[idir] = geom[lev].CellSize()[idir];

      for (int i = 1; i <= 500; i++) {
        mfix_get_walls(&i, &exists, &normal, &center);
        if(exists){
          amrex::Print() << "Normal " << normal << std::endl;
          amrex::Print() << "Center " << center << std::endl;
          plane = new AnisotropicDxPlaneIF(normal,center,true,dxVec);
          planes.push_back(plane);
        }
      }
      IntersectionIF all_planes(planes);

      impfunc_walls = std::unique_ptr<BaseIF>(all_planes.newImplicitFunction());

      impfunc.reset(all_planes.newImplicitFunction());
    }


    int max_level = 0;
    int grid_size = 16;
    bool eb_verbosity = true;

    /******************************************************************************************************************
     *                                                                                                                *
     * Fill Level-set using:                                                                                          *
     *      -> Planes (where the GeometryShop's implicit function is a signed distance): implicit function's value    *
     *      -> Poly2 (where GeometryShop's implicit function is singed but not a distance): min distance to EB facets *
     * Note: this requires building and destroying the EBTower (twice), so any EBTower data built before this will be *
     * lost...                                                                                                        *
     *                                                                                                                *
     ******************************************************************************************************************/

    if(use_walls){
        // Define both components of the GeometryShop separately:
        GeometryShop gshop_walls(* impfunc_walls, eb_verbosity);

        // Define the EBIS first using only the walls...
        // Note GeometryShop's behaviour wrt anisotropic cells: * use x-component of dx as reference length-scale
        //                                                      * rescale y, z- components wrt to dx[0] (dx(1))
        Geometry geom_ls = LSUtility::make_ls_geometry(* level_set);
        AMReX_EBIS::instance()->define(geom_ls.Domain(), RealVect::Zero, geom_ls.CellSize()[0], gshop_walls, grid_size, max_level);

        EBTower::Build();
        // GeometryShop's Planes' implicit function is actually a signed distance function
        //      => it's just easier to fill the level-set this way
        level_set->intersection_ebis(* AMReX_EBIS::instance());
        EBTower::Destroy();
    }
    if(use_poly2){
        // Define both components of the GeometryShop separately:
        GeometryShop gshop_poly2(* impfunc_poly2, eb_verbosity);

        // Define the EBIS using only the poly2 (after deleting the walls-only EBTower)...
        // Note GeometryShop's behaviour wrt anisotropic cells: * use x-component of dx as reference length-scale
        //                                                      * rescale y, z- components wrt to dx[0] (dx(1)
        Geometry geom_eb = LSUtility::make_eb_geometry(* level_set);
        AMReX_EBIS::instance()->define(geom_eb.Domain(), RealVect::Zero, geom_eb.CellSize()[0], gshop_poly2, grid_size, max_level);

        EBTower::Build();
        // GeometryShop's PolynomialIF is not a signed distance function...
        //      => it's easier to use PolynomialIF to build an EBFArrayBoxFactory which defines our EB surface now
        //          => define the level set as the (signed) distance to the closest point on the EB-facets
        int eb_pad = level_set->get_eb_pad();
        EBFArrayBoxFactory eb_factory_poly2(geom_eb, level_set->get_eb_ba(), dmap[lev], {eb_pad, eb_pad, eb_pad}, EBSupport::full);
        level_set->intersection_ebf(eb_factory_poly2, * AMReX_EBIS::instance());
        EBTower::Destroy();
    }

    // store copy of level set (for later use).
    ls[lev] = level_set->copy_data();

    /******************************************************************************************************************
     *                                                                                                                *
     * Build standard EB Factories                                                                                    *
     *                                                                                                                *
     ******************************************************************************************************************/


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

    eb_normals         = pc -> EBNormals(lev, particle_ebfactory.get(), dummy.get());
}



void
mfix_level::make_eb_hourglass(int lev)
{
    if (geom[lev].isAllPeriodic()) return;

    Box domain(geom[lev].Domain());
    Real dx = geom[lev].CellSize()[0];

    int exists;
    RealVect normal, center;
    PlaneIF* plane;
    Vector<BaseIF*> planes;
    planes.resize(0);

    std::unique_ptr<BaseIF> impfunc;
    std::unique_ptr<BaseIF> impfunc_unpolys;
    std::unique_ptr<BaseIF> impfunc_walls;

    ParmParse pp("mfix");

    amrex::Print() << "Using poly geometry\n";

    Vector<PolyTerm> poly1;

    PolyTerm mono;
    Real coef;
    IntVect powers;

    Vector<Real> coefvec(5);
    Vector<int>  powersvec(5);
    Vector<Real> transvec(SpaceDim);

    for(int idir = 0; idir < 3; idir++) {

      if( idir == 0) {
        pp.getarr("poly1_x_coeffs", coefvec, 0, 5);
      } else if( idir == 1) {
        pp.getarr("poly1_y_coeffs", coefvec, 0, 5);
      } else if( idir == 2) {
        pp.getarr("poly1_z_coeffs", coefvec, 0, 5);
      }

      for(int lc = 0; lc < 5; lc++) {

        // x^(lc) term
        coef = coefvec[lc];
        powers = IntVect::Zero;
        powers[idir] = lc;

        mono.coef   = coef;
        mono.powers = powers;

        poly1.push_back(mono);

      }
    }

    bool flip1 = false;
    pp.query("poly1_mirror", flip1);

    PolynomialIF mirror1(poly1,flip1);
    RealVect translation1;

    pp.getarr("poly1_translate", transvec,  0, SpaceDim);

    for(int idir = 0; idir < 3; idir++) {
      translation1[idir] = transvec[idir];
    }

    TransformIF poly1t(mirror1);
    poly1t.translate(translation1);

    Vector<PolyTerm> poly2;

    for(int idir = 0; idir < 3; idir++) {

      if( idir == 0) {
        pp.getarr("poly2_x_coeffs", coefvec, 0, 5);
      } else if( idir == 1) {
        pp.getarr("poly2_y_coeffs", coefvec, 0, 5);
      } else if( idir == 2) {
        pp.getarr("poly2_z_coeffs", coefvec, 0, 5);
      }

      for(int lc = 0; lc < 5; lc++) {

        // x^(lc) term
        coef = coefvec[lc];
        powers = IntVect::Zero;
        powers[idir] = lc;

        mono.coef   = coef;
        mono.powers = powers;

        poly2.push_back(mono);

      }
    }

    bool flip2 = false;
    pp.query("poly2_mirror", flip2);

    PolynomialIF mirror2(poly2,flip2);
    RealVect translation2;

    pp.getarr("poly2_translate", transvec,  0, SpaceDim);

    for(int idir = 0; idir < 3; idir++) {
      translation2[idir] = transvec[idir];
    }

    TransformIF poly2t(mirror2);
    poly2t.translate(translation2);

    for (int i = 1; i <= 500; i++) {
      mfix_get_walls(&i, &exists, &normal, &center);
      if(exists){
        amrex::Print() << "Normal " << normal << std::endl;
        amrex::Print() << "Center " << center << std::endl;
        plane = new PlaneIF(normal,center,true);
        planes.push_back(plane);
      }
    }
    IntersectionIF all_planes(planes);

    impfunc_walls = std::unique_ptr<BaseIF>(all_planes.newImplicitFunction());

    Vector<BaseIF*> bulbs(2);
    bulbs[0] = &poly1t;
    bulbs[1] = &poly2t;
    UnionIF unpolys(bulbs);

    impfunc_unpolys = std::unique_ptr<BaseIF>(unpolys.newImplicitFunction());


    Vector<BaseIF*> funcs(2);
    funcs[0] = &unpolys;
    funcs[1] = &all_planes;
    IntersectionIF implicit(funcs);
    impfunc.reset(implicit.newImplicitFunction());


    int max_level = 0;
    int grid_size = 16;
    bool eb_verbosity = true;



    /******************************************************************************************************************
     *                                                                                                                *
     * Fill Level-set using:                                                                                          *
     *      -> Planes (where the GeometryShop's implicit function is a signed distance): implicit function's value    *
     *      -> Poly (where GeometryShop's implicit function is singed but not a distance): min distance to EB facets  *
     * Note: this requires building and destroying the EBTower (twice), so any EBTower data built before this will be *
     * lost...                                                                                                        *
     *                                                                                                                *
     ******************************************************************************************************************/

    // Define both components of the GeometryShop separately:
    GeometryShop gshop_upoly(* impfunc_unpolys, eb_verbosity);
    GeometryShop gshop_walls(* impfunc_walls, eb_verbosity);


    // Define the EBIS first using only the walls...
    Geometry geom_ls = LSUtility::make_ls_geometry(* level_set);
    AMReX_EBIS::instance()->define(geom_ls.Domain(), RealVect::Zero, geom_ls.CellSize()[0], gshop_walls, grid_size, max_level);

    EBTower::Build();
    // GeometryShop's Planes' implicit function is actually a signed distance function
    //      => it's just easier to fill the level-set this way
    level_set->intersection_ebis(* AMReX_EBIS::instance());
    EBTower::Destroy();

    // Define the EBIS using only the poly (after deleting the walls-only EBTower)...
    Geometry geom_eb = LSUtility::make_eb_geometry(* level_set);
    AMReX_EBIS::instance()->define(geom_eb.Domain(), RealVect::Zero, geom_eb.CellSize()[0], gshop_upoly, grid_size, max_level);

    EBTower::Build();
    // GeometryShop's PolynomialIF is not a signed distance function...
    //      => it's easier to use PolynomialIF to build an EBFArrayBoxFactory which defines our EB surface now
    //          => define the level set as the (signed) distance to the closest point on the EB-facets
    int eb_pad = level_set->get_eb_pad();
    EBFArrayBoxFactory eb_factory_poly(geom_eb, level_set->get_eb_ba(), dmap[lev], {eb_pad, eb_pad, eb_pad}, EBSupport::full);
    level_set->intersection_ebf(eb_factory_poly, * AMReX_EBIS::instance());
    EBTower::Destroy();

    // store copy of level set (for later use).
    ls[lev] = level_set->copy_data();

    /******************************************************************************************************************
     *                                                                                                                *
     * Build standard EB Factories                                                                                    *
     *                                                                                                                *
     ******************************************************************************************************************/

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

    eb_normals         = pc -> EBNormals(lev, particle_ebfactory.get(), dummy.get());
}



void
mfix_level::make_eb_clr(int lev)
{
    Box domain(geom[lev].Domain());
    Real dx = geom[lev].CellSize()[0];

    RealVect normal, center;
    Vector<BaseIF*> planes;
    planes.resize(0);

    std::unique_ptr<BaseIF> impfunc;
    std::unique_ptr<BaseIF> impfunc_poly2;

    ParmParse pp("mfix");

    amrex::Print() << "Creating CLR geometry\n";

    Vector<PolyTerm> poly;

    PolyTerm mono;
    IntVect powers;

    Vector<Real> coefvec(3);
    Vector<int>  powersvec(3);
    Vector<Real> transvec(3);
    RealVect translation;

    Real lradius, lheight;

    int cylinder_dir;

    // mfix_level::make_cylinder uses a union
    //  => ensure that mfix_level::level_set is initialized to min
    level_set->invert();

    ct_ls_mf = 0;
    std::unique_ptr<MultiFab> ls_mf = level_set->copy_data();
    amrex::VisMF::Write(* ls_mf, "ls_empty");

    //------------------------------------------------------------- Riser
    pp.getarr("riser_translate", transvec,  0, 3);

    for(int idir = 0; idir < 3; idir++)
        translation[idir] = transvec[idir];

    pp.query("riser_lower_radius", lradius);
    pp.query("riser_lower_height", lheight);

    cylinder_dir=1;
    std::unique_ptr<BaseIF> riser_lower = make_cylinder(cylinder_dir, lradius, lheight, translation, lev);


    pp.query("riser_upper_radius", lradius);
    pp.query("riser_upper_height", lheight);

    cylinder_dir=1;
    std::unique_ptr<BaseIF> riser_upper = make_cylinder(cylinder_dir, lradius, lheight, translation, lev);

    //------------------------------------------------------------- J-leg
    pp.getarr("jleg_htranslate", transvec,  0, 3);

    for(int idir = 0; idir < 3; idir++)
        translation[idir] = transvec[idir];

    pp.query("jleg_horz_radius", lradius);
    pp.query("jleg_horz_height", lheight);

    cylinder_dir=0;
    std::unique_ptr<BaseIF> jleg_horz = make_cylinder(cylinder_dir, lradius, lheight, translation, lev);


    pp.getarr("jleg_vtranslate", transvec,  0, 3);

    for(int idir = 0; idir < 3; idir++)
        translation[idir] = transvec[idir];

    pp.query("jleg_vert_radius", lradius);
    pp.query("jleg_vert_height", lheight);

    cylinder_dir=1;
    std::unique_ptr<BaseIF> jleg_vert = make_cylinder(cylinder_dir, lradius, lheight, translation, lev);

    //------------------------------------------------------------- reactor
    pp.getarr("reactor_translate", transvec,  0, 3);

    for(int idir = 0; idir < 3; idir++)
        translation[idir] = transvec[idir];

    pp.query("reactor_radius", lradius);
    pp.query("reactor_height", lheight);

    cylinder_dir=1;
    std::unique_ptr<BaseIF> reactor = make_cylinder(cylinder_dir, lradius, lheight, translation, lev);

    //------------------------------------------------------------- loop-seal to reactor
    pp.getarr("ls2fbr_translate", transvec,  0, 3);

    for(int idir = 0; idir < 3; idir++)
        translation[idir] = transvec[idir];

    pp.query("ls2fbr_radius", lradius);
    pp.query("ls2fbr_height", lheight);

    cylinder_dir=1;
    std::unique_ptr<BaseIF> ls2fbr = make_cylinder(cylinder_dir, lradius, lheight, translation, lev);


    //------------------------------------------------------------- loop-seal to reactor
    pp.getarr("loopseal_translate", transvec,  0, 3);

    for(int idir = 0; idir < 3; idir++)
        translation[idir] = transvec[idir];

    pp.query("loopseal_radius", lradius);
    pp.query("loopseal_height", lheight);

    cylinder_dir=1;
    std::unique_ptr<BaseIF> loopseal = make_cylinder(cylinder_dir, lradius, lheight, translation, lev);



    //------------------------------------------------------------- loop-seal to reactor
    pp.getarr("cy2ls_translate", transvec,  0, 3);

    for(int idir = 0; idir < 3; idir++)
        translation[idir] = transvec[idir];

    pp.query("cy2ls_radius", lradius);
    pp.query("cy2ls_height", lheight);

    cylinder_dir=1;
    std::unique_ptr<BaseIF> cy2ls = make_cylinder(cylinder_dir, lradius, lheight, translation, lev);


    //------------------------------------------------------------- cyclone
    pp.getarr("cyclone_translate", transvec,  0, 3);

    for(int idir = 0; idir < 3; idir++)
        translation[idir] = transvec[idir];

    pp.query("cyclone_radius", lradius);
    pp.query("cyclone_height", lheight);

    cylinder_dir=1;
    std::unique_ptr<BaseIF> cyclone = make_cylinder(cylinder_dir, lradius, lheight, translation, lev);


    //------------------------------------------------------------- crossover
    pp.getarr("crossover_translate", transvec,  0, 3);

    for(int idir = 0; idir < 3; idir++)
        translation[idir] = transvec[idir];

    pp.query("crossover_radius", lradius);
    pp.query("crossover_height", lheight);

    cylinder_dir=0;
    std::unique_ptr<BaseIF> crossover = make_cylinder(cylinder_dir, lradius, lheight, translation, lev);

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

    // store copy of level set (for later use).
    ls[lev] = level_set->copy_data();


    impfunc_poly2 = std::unique_ptr<BaseIF>(clr.newImplicitFunction());
    impfunc.reset(clr.newImplicitFunction());


    int max_level = 0;
    int grid_size = 16;
    bool eb_verbosity = true;

    /**************************************************************************************
     *                                                                                    *
     * Fill Level-set using:                                                              *
     *      -> Poly2 (where GeometryShop's implicit function is singed but not a          *
     *         distance): min distance to EB facets.                                      *
     *                                                                                    *
     * Note: this requires building and destroying the EBTower (twice), so any EBTower    *
     *       data built before this will be lost....                                      *
     *                                                                                    *
     **************************************************************************************/

    // Define both components of the the GeometryShop seperately:
    //GeometryShop gshop_poly2(* impfunc_poly2, eb_verbosity);

    // Define the EBIS using only the poly2 (after deleting the walls-only EBTower)...
    //Geometry geom_eb = LSUtility::make_eb_geometry(* level_set);
    //AMReX_EBIS::instance()->define(geom_eb.Domain(), RealVect::Zero, geom_eb.CellSize()[0], gshop_poly2, grid_size, max_level);

    //EBTower::Build();
    // GeometryShop's PolynomialIF is not a signed distance function...
    //      => it's easier to use PolynomialIF to build an EBFArrayBoxFactory which defines our EB surface now
    //          => define the level set as the (signed) distance to the closest point on the EB-facets
    //int eb_pad = level_set->get_eb_pad();
    //EBFArrayBoxFactory eb_factory_poly2(geom_eb, level_set->get_eb_ba(), dmap[lev], {eb_pad, eb_pad, eb_pad}, EBSupport::full);
    // level_set->update_ebf(eb_factory_poly2, * AMReX_EBIS::instance());
    //EBTower::Destroy();



    /*****************************************************************************
     *                                                                           *
     * Build standard EB Factories                                               *
     *                                                                           *
     *****************************************************************************/


    GeometryShop gshop(*impfunc, eb_verbosity);
    AMReX_EBIS::instance()->define(domain, RealVect::Zero, dx, gshop, grid_size, max_level);

    // set up ebfactory
    int m_eb_basic_grow_cells = 16;
    int m_eb_volume_grow_cells = 16;
    int m_eb_full_grow_cells = 16;
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

    //level_set->intersection_ebf(* ebfactory, * AMReX_EBIS::instance());

    // store copy of level set (for later use).
    ls[lev] = level_set->copy_data();

    eb_normals         = pc -> EBNormals(lev, particle_ebfactory.get(), dummy.get());
}



std::unique_ptr<BaseIF>
mfix_level::make_cylinder(int dir, Real radius, Real length, const RealVect & translation, int lev)
{
    std::unique_ptr<BaseIF> cylinder_IF;
    Vector<PolyTerm> poly;

    PolyTerm mono;
    Real coef;
    IntVect powers;

    Vector<Real> coefvec(3);
    Vector<int>  powersvec(3);

    for(int idir = 0; idir < 3; idir++) {

        if( idir == dir) {
            coefvec[0] = -radius*radius;
            coefvec[1] = 0.0;
            coefvec[2] = 0.0;
        } else {
            coefvec[0] = 0.0;
            coefvec[1] = 0.0;
            coefvec[2] = 1.0;
        }

        for(int lc = 0; lc < 3; lc++) {
            // x^(lc) term
            coef = coefvec[lc];
            powers = IntVect::Zero;
            powers[idir] = lc;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);
        }
    }

    // Internal flow cylinder
    PolynomialIF cylinder0(poly, true);
    TransformIF cylinder1(cylinder0);
    cylinder1.translate(translation);

    // box to clip to correct length
    RealVect normal, center;
    PlaneIF* plane;
    Vector<BaseIF*> planes;
    planes.resize(0);

    for(int i=0; i<3; i++) {
        center[i] = 0.0;
        normal[i] = 0.0;
    }

    center[dir] = 0.0;
    normal[dir] = 1.0;

    // amrex::Print() << "Plane 1\n";
    // amrex::Print() << "Center " << center  << "\n";
    // amrex::Print() << "Normal " << normal  << "\n";

    plane = new PlaneIF(normal,center,true);
    planes.push_back(plane);

    center[dir] = length;
    normal[dir] =-1.0;

    // amrex::Print() << "Plane 2\n";
    // amrex::Print() << "Center " << center  << "\n";
    // amrex::Print() << "Normal " << normal  << "\n";

    plane = new PlaneIF(normal,center,true);
    planes.push_back(plane);

    IntersectionIF bounding_box(planes);
    TransformIF walls(bounding_box);
    walls.translate(translation);

    Vector<BaseIF*> funcs(2);
    funcs[0] = &cylinder0;
    funcs[1] = &bounding_box;

    IntersectionIF cylinder(funcs);

    TransformIF cylinder_trans(cylinder);
    cylinder_trans.translate(translation);

    cylinder_IF.reset(cylinder_trans.newImplicitFunction());


    int max_level = 0;
    int grid_size = 16;
    bool eb_verbosity = true;


    /******************************************************************************************************************
     *                                                                                                                *
     * Fill Level-set using:                                                                                          *
     *      -> Planes (where the GeometryShop's implicit function is a signed distance): implicit function's value    *
     *      -> Poly (where GeometryShop's implicit function is singed but not a distance): min distance to EB facets  *
     * Note: this requires building and destroying the EBTower (twice), so any EBTower data built before this will be *
     * lost...                                                                                                        *
     *                                                                                                                *
     ******************************************************************************************************************/

    // Define both components of the GeometryShop separately:
    GeometryShop gshop_upoly(cylinder1, eb_verbosity);
    GeometryShop gshop_walls(walls, eb_verbosity);

    // Define a temporary level-set used for constructing the cylinder:
    LSFactory ls_cylinder(* level_set);
    LSFactory ls_walls(ls_cylinder);

    // Define the EBIS first using only the walls...
    Geometry geom_ls = LSUtility::make_eb_geometry(ls_cylinder);
    AMReX_EBIS::instance()->define(geom_ls.Domain(), RealVect::Zero, geom_ls.CellSize()[0], gshop_walls, grid_size, max_level);

    EBTower::Build();
    // GeometryShop's Planes' implicit function is actually a signed distance function
    //      => it's just easier to fill the level-set this way
    ls_cylinder.intersection_ebis(* AMReX_EBIS::instance());
    EBTower::Destroy();

    ct_ls_mf ++;
    std::stringstream ss1;
    ss1 << "ls_" << ct_ls_mf << "_a";
    std::unique_ptr<MultiFab> ls_mf_a = ls_cylinder.copy_data();
    amrex::VisMF::Write(* ls_mf_a, ss1.str());

    // Define the EBIS using only the poly (after deleting the walls-only EBTower)...
    Geometry geom_eb = LSUtility::make_eb_geometry(ls_cylinder);
    AMReX_EBIS::instance()->define(geom_eb.Domain(), RealVect::Zero, geom_eb.CellSize()[0], gshop_upoly, grid_size, max_level);

    EBTower::Build();
    // GeometryShop's PolynomialIF is not a signed distance function...
    //      => it's easier to use PolynomialIF to build an EBFArrayBoxFactory which defines our EB surface now
    //          => define the level set as the (signed) distance to the closest point on the EB-facets
    int eb_pad = level_set->get_eb_pad();
    EBFArrayBoxFactory eb_factory_poly(geom_eb, level_set->get_eb_ba(), dmap[lev], {eb_pad, eb_pad, eb_pad}, EBSupport::full);
    ls_cylinder.intersection_ebf(eb_factory_poly, * AMReX_EBIS::instance());
    ls_walls.intersection_ebf(eb_factory_poly, * AMReX_EBIS::instance());
    EBTower::Destroy();

    std::stringstream ss2;
    ss2 << "ls_" << ct_ls_mf << "_b";
    std::unique_ptr<MultiFab> ls_mf_b = ls_walls.copy_data();
    amrex::VisMF::Write(* ls_mf_b, ss2.str());

    std::stringstream ss;
    ss << "ls_" << ct_ls_mf;
    std::unique_ptr<MultiFab> ls_mf = ls_cylinder.copy_data();
    amrex::VisMF::Write(* ls_mf, ss.str());

    level_set->update_union(* ls_cylinder.get_data());

    return cylinder_IF;
}
