//#include <AMReX_VisMF.H>  // amrex::VisMF::Write(MultiFab)
//#include <sstream>

#include <algorithm>
#include <AMReX_EB_levelset.H>
#include <mfix_level.H>
#include <mfix_eb_F.H>



void mfix_level::fill_levelset(int lev, bool use_walls, bool use_poly,
                               const WallsIF & impfunc_walls, const PolynomialIF & impfunc_poly,
                               int max_level_here, int grid_size, bool eb_verbosity)
{
    // Do nothing if loading level-set from restart file:
    if(levelset__restart) return;


    /****************************************************************************
    *                                                                           *
    * Fill Level-set using:                                                     *
    *      -> Planes (where the GeometryShop's implicit function is a signed    *
    *         distance): implicit function's value                              *
    *      -> Poly (where GeometryShop's implicit function is singed but not    *
    *         a distance): min distance to EB facets                            *
    * Note: this requires building and destroying the EBTower (twice), so any   *
    * EBTower data built before this will be lost...                            *
    *                                                                           *
    *****************************************************************************/

    Geometry geom_eb = LSUtility::make_eb_geometry(* level_set, geom[lev]);
    if(use_walls){
        // Define components of the GeometryShop separately:
        EB2::GeometryShop<WallsIF> gshop_walls(impfunc_walls);
        GShopLSFactory<WallsIF>    ls_gshop_walls(gshop_walls, * level_set);

        // Implicit function used by LSFactory => returned MF has the same DM as LSFactory
        std::unique_ptr<MultiFab>  mf_impfunc_walls = ls_gshop_walls.fill_impfunc();

        // Plane implicit function is already a signed distance function => it's
        // just easier to fill the level-set this way
        level_set->intersection_impfunc(* mf_impfunc_walls);
    }
    if(use_poly){
        // Define components of the GeometryShop separately:
        EB2::GeometryShop<PolynomialIF> gshop_poly(impfunc_poly);
        GShopLSFactory<PolynomialIF>    ls_gshop_poly(gshop_poly, * level_set);

        // Implicit function used by LSFactory => returned MF has the same DM as LSFactory
        std::unique_ptr<MultiFab>  mf_impfunc_poly = ls_gshop_poly.fill_impfunc();

        EB2::Build(gshop_poly, geom[lev], max_level_here, max_level_here);

        const EB2::IndexSpace & poly_ebis = EB2::IndexSpace::top();
        const EB2::Level &      poly_lev  = poly_ebis.getLevel(geom[lev]);

        // PolynomialIF is not a signed distance function... => it's easier to
        // use PolynomialIF to build an EBFArrayBoxFactory which defines our EB
        // surface now => define the level set as the (signed) distance to the
        // closest point on the EB-facets
        int eb_pad = level_set->get_eb_pad();
        EBFArrayBoxFactory eb_factory_poly(poly_lev, geom_eb,
                                           pc->ParticleBoxArray(lev), pc->ParticleDistributionMap(lev),
                                           /* level_set->get_eb_ba(), level_set->get_eb_dm(), */
                                           {eb_pad, eb_pad, eb_pad}, EBSupport::full);

        // According to Weiqun `fix_up_geometry` is deprecated.
        //fix_up_geometry(& eb_factory_poly, lev);

        level_set->intersection_ebf(eb_factory_poly, * mf_impfunc_poly);
    }
}

