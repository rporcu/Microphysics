
#include <algorithm>
#include <AMReX_EB_levelset.H>
#include <mfix_level.H>
#include <mfix_eb_F.H>

#if 0
void mfix_level::fill_levelset(int lev, bool use_walls, bool use_poly,
                               const BaseIF & impfunc_walls, const BaseIF & impfunc_poly,
                               int max_level_here, int grid_size, bool eb_verbosity) 
{
    // Do nothing if loading level-set from restart file:
    if(levelset__restart) return;


   /****************************************************************************
    *                                                                          *
    * Fill Level-set using:                                                    *
    *      -> Planes (where the GeometryShop's implicit function is a signed   *
    *         distance): implicit function's value                             *
    *      -> Poly (where GeometryShop's implicit function is singed but not   *
    *         a distance): min distance to EB facets                           *
    * Note: this requires building and destroying the EBTower (twice), so any  *
    * EBTower data built before this will be lost...                           *
    *                                                                          *
    ****************************************************************************/

    Geometry geom_eb = LSUtility::make_eb_geometry(* level_set, geom[lev]);
    if(use_walls){
        // Define components of the GeometryShop separately:
        GeometryShop gshop_walls(impfunc_walls, eb_verbosity);

        // Define the EBIS first using only the walls...
        AMReX_EBIS::instance()->define(geom_eb.Domain(),
                                       RealVect::Zero,  // ......... origin of EBIndexSpace
                                       geom_eb.CellSize()[0],  // .. reference cell size of EBIndexSpace [1]
                                       gshop_walls,  // ............ GeometryShop object
                                       grid_size, max_level_here);
        // [1]: EBIndexSpace internally assumes an isotropic grid. Any
        // anisotropic implicit function (e.g AnisotrpicPlaneIF) uses dx as a
        // reference, and rescales dy and dz wrt dx. => dx goes here.

        EBTower::Build();
        // GeometryShop's Planes' implicit function is actually a signed distance function
        //      => it's just easier to fill the level-set this way
        level_set->intersection_ebis(* AMReX_EBIS::instance()); EBTower::Destroy();
    }
    if(use_poly){
        // Define components of the GeometryShop separately:
        GeometryShop gshop_poly(impfunc_poly, eb_verbosity);

        // Define the EBIS using only the poly2 (after deleting the walls-only EBTower)...
        AMReX_EBIS::instance()->define(geom_eb.Domain(),
                                       RealVect::Zero,  // ......... origin of EBIndexSpace
                                       geom_eb.CellSize()[0],  // .. reference cell size of EBIndexSpace [1]
                                       gshop_poly,  // ............ GeometryShop object
                                       grid_size, max_level_here);

        EBTower::Build();
        // GeometryShop's PolynomialIF is not a signed distance function...
        //      => it's easier to use PolynomialIF to build an
        //         EBFArrayBoxFactory which defines our EB surface now
        //          => define the level set as the (signed) distance to the
        //             closest point on the EB-facets
        int eb_pad = level_set->get_eb_pad();
        EBFArrayBoxFactory eb_factory_poly2(geom_eb, pc->ParticleBoxArray(lev), pc->ParticleDistributionMap(lev),
                                            /* level_set->get_eb_ba(), level_set->get_eb_dm(), */
                                            {eb_pad, eb_pad, eb_pad}, EBSupport::full);
        fix_up_geometry(&eb_factory_poly2,lev);
        level_set->intersection_ebf(eb_factory_poly2, * AMReX_EBIS::instance());
        EBTower::Destroy();
    }
}
#endif
