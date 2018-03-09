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
            center[0] = 1.e-3;
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
        // Define both components of the the GeometryShop seperately:
        GeometryShop gshop_walls(* impfunc_walls, eb_verbosity);

        // Define the EBIS first using only the walls...
        Geometry geom_ls = LSUtility::make_ls_geometry(* level_set);
        AMReX_EBIS::instance()->define(geom_ls.Domain(), RealVect::Zero, geom_ls.CellSize()[0], gshop_walls, grid_size, max_level);

        EBTower::Build();
        // GeometryShop's Planes' implicit function is actually a signed distance function
        //      => it's just easier to fill the level-set this way
        level_set->update_ebis(* AMReX_EBIS::instance());
        EBTower::Destroy();
    }
    if(use_poly2){
        // Define both components of the the GeometryShop seperately:
        GeometryShop gshop_poly2(* impfunc_poly2, eb_verbosity);

        // Define the EBIS using only the poly2 (after deleting the walls-only EBTower)...
        Geometry geom_eb = LSUtility::make_eb_geometry(* level_set);
        AMReX_EBIS::instance()->define(geom_eb.Domain(), RealVect::Zero, geom_eb.CellSize()[0], gshop_poly2, grid_size, max_level);

        EBTower::Build();
        // GeometryShop's PolynomialIF is not a signed distance function... 
        //      => it's easier to use PolynomialIF to build an EBFArrayBoxFactory which defines our EB surface now
        //          => define the level set as the (signed) distance to the closest point on the EB-facets
        EBFArrayBoxFactory eb_factory_poly2(geom_eb, level_set->get_eb_ba(), dmap[lev], {2, 2, 2}, EBSupport::full);
        level_set->update_ebf(eb_factory_poly2, * AMReX_EBIS::instance());
        EBTower::Destroy();
    }



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

    // Define both components of the the GeometryShop seperately:
    GeometryShop gshop_upoly(* impfunc_unpolys, eb_verbosity);
    GeometryShop gshop_walls(* impfunc_walls, eb_verbosity);


    // Define the EBIS first using only the walls...
    Geometry geom_ls = LSUtility::make_ls_geometry(* level_set);
    AMReX_EBIS::instance()->define(geom_ls.Domain(), RealVect::Zero, geom_ls.CellSize()[0], gshop_walls, grid_size, max_level);

    EBTower::Build();
    // GeometryShop's Planes' implicit function is actually a signed distance function
    //      => it's just easier to fill the level-set this way
    level_set->update_ebis(* AMReX_EBIS::instance());
    EBTower::Destroy();

    // Define the EBIS using only the poly (after deleting the walls-only EBTower)...
    Geometry geom_eb = LSUtility::make_eb_geometry(* level_set);
    AMReX_EBIS::instance()->define(geom_eb.Domain(), RealVect::Zero, geom_eb.CellSize()[0], gshop_upoly, grid_size, max_level);

    EBTower::Build();
    // GeometryShop's PolynomialIF is not a signed distance function... 
    //      => it's easier to use PolynomialIF to build an EBFArrayBoxFactory which defines our EB surface now
    //          => define the level set as the (signed) distance to the closest point on the EB-facets
    EBFArrayBoxFactory eb_factory_poly(geom_eb, level_set->get_eb_ba(), dmap[lev], {2, 2, 2}, EBSupport::full);
    level_set->update_ebf(eb_factory_poly, * AMReX_EBIS::instance());
    EBTower::Destroy();



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
