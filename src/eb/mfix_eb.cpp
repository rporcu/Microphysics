#include <algorithm>
#include <AMReX_EB_levelset.H>
#include <mfix.H>
#include <mfix_eb_F.H>

void mfix::make_eb_geometry ()
{

    /****************************************************************************
     *                                                                          *
     * mfix.geometry=<string> specifies the EB geometry. <string> can be on of  *
     * box, cylinder, hopper, clr, clr_riser, general (or blank)                *
     *                                                                          *
     ***************************************************************************/


    ParmParse pp("mfix");

    std::string geom_type;
    pp.query("geometry", geom_type);

    /****************************************************************************
     *                                                                          *
     * Legacy inputs:                                                           *
     *   -- mfix.hourglass = true <=> mfix.geometry=box                         *
     *   -- mfix.clr       = true <=> mfix.geometry=clr                         *
     *   -- mfix.clr_riser = true <=> mfix.geometry=clr_riser                   *
     *   -- mfix.use_walls = true <=> mfix.geometry=general                     *
     *   -- mfix.use_poy2  = true <=> mfix.geometry=general                     *
     *                                                                          *
     ***************************************************************************/


    bool hourglass    = false;
    bool clr          = false;
    bool clr_riser    = false;
    bool eb_general   = false;

    pp.query("hourglass", hourglass);
    pp.query("clr", clr);
    pp.query("clr_riser", clr_riser);

    bool eb_poly2 = false;
    bool eb_walls = false;

    pp.query("use_poly2", eb_poly2);
    pp.query("use_walls", eb_walls);
    eb_general = eb_poly2 || eb_walls;

    // Avoid multiple (ambiguous) inputs
    if (hourglass || clr || clr_riser || eb_general) {
        if (! geom_type.empty()) {
            amrex::Abort("The input file cannot specify both:\n"
                         "mfix.<geom_type>=true and mfix.geometry=<geom_type>\n"
                         "at the same time."                                   );
        }
    }

    if (hourglass)  geom_type = "hourglass";
    if (clr)        geom_type = "clr";
    if (clr_riser)  geom_type = "clr_riser";
    if (eb_general) geom_type = "general";


    /****************************************************************************
     *                                                                          *
     *  CONSTRUCT EB                                                            *
     *                                                                          *
     ***************************************************************************/

    if (geom_type == "box") {
        amrex::Print() << "\n Building box geometry." << std::endl;
        make_eb_box();
    } else if (geom_type == "cylinder") {
        amrex::Print() << "\n Building cylinder geometry." << std::endl;
        make_eb_cylinder();
    } else if (geom_type == "hopper") {
        amrex::Print() << "\n Building hopper geometry." << std::endl;
        make_eb_hopper();
    } else if (geom_type == "cyclone") {
        amrex::Print() << "\n Building cyclone geometry." << std::endl;
        make_eb_cyclone();
    } else if(geom_type == "general") {
        amrex::Print() << "\n Building general geometry (poly2 with extra walls)." << std::endl;
        make_eb_general();
    } else {
        amrex::Print() << "\n No EB geometry declared in inputs => "
                       << " Will read walls from mfix.dat only."
                       << std::endl;
        make_eb_regular();
    }
}


void mfix::make_eb_factories () {

    /****************************************************************************
     *                                                                          *
     * Fill EB factories as an initial run. Since the particle container might  *
     * not have been created yet, use mfix grids instead.                       *
     *                                                                          *
     ***************************************************************************/

    for (int lev = 0; lev < nlev; lev++) {
        ebfactory[lev].reset(
            new EBFArrayBoxFactory(* eb_levels[lev], geom[lev], grids[lev], dmap[lev],
                                   {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                    m_eb_full_grow_cells}, m_eb_support_level)
            );

        // Grow EB factory by +2 in order to avoid edge cases. This is not
        // necessary for multi-level mfix.
        particle_ebfactory[lev].reset(
            new EBFArrayBoxFactory(* eb_levels[lev], geom[lev], grids[lev], dmap[lev],
                                   {levelset__eb_pad + 2, levelset__eb_pad + 2,
                                    levelset__eb_pad + 2}, m_eb_support_level)
            );
    }
}


void mfix::fill_eb_levelsets ()
{
    const DistributionMapping & part_dm = pc->ParticleDistributionMap(0);
    const BoxArray &            part_ba = pc->ParticleBoxArray(0);


    LSFactory lsf(0, levelset__refinement, levelset__eb_refinement,
                  levelset__pad, levelset__eb_pad, part_ba, geom[0], part_dm );

    MultiFab impfunc = MFUtil::regrid(lsf.get_ls_ba(), part_dm, * implicit_functions[1], true);
    lsf.Fill( * particle_ebfactory[0], impfunc);

    level_sets[1] = lsf.copy_data(part_dm);
    level_sets[0] = lsf.coarsen_data();

    // if (nlev == 1)
    // {
    //     const DistributionMapping & part_dm = pc->ParticleDistributionMap(0);
    //     const BoxArray &            part_ba = pc->ParticleBoxArray(0);

    //     //_______________________________________________________________________
    //     // Baseline Level-Set
    //     {
    //         // NOTE: reference BoxArray is not nodal
    //         BoxArray ba = amrex::convert(part_ba, IntVect::TheNodeVector());

    //         level_sets[0].reset(new MultiFab);
    //         level_sets[0]->define(ba, part_dm, 1, levelset__pad);
    //         iMultiFab valid(ba, part_dm, 1, levelset__pad);

    //         // NOTE: implicit function data might not be on the right grids
    //         MultiFab impfunc = MFUtil::regrid(ba, part_dm, * implicit_functions[0], true);

    //         LSFactory::fill_data(* level_sets[0], valid, * particle_ebfactory[0], impfunc,
    //                              32, 1, 1, geom[0], geom[0]);
    //     }

    //     //_______________________________________________________________________
    //     // Refined Level-Set
    //     // TODO: Don't actually refine this thing if levelset refinement is 1
    //     {
    //         // Set up refined geometry
    //         Box dom = geom[0].Domain();
    //         dom.refine(levelset__refinement);
    //         Geometry geom_lev(dom);

    //         // Set up refined BoxArray. NOTE: reference BoxArray is not nodal
    //         BoxArray ba = amrex::convert(part_ba, IntVect::TheNodeVector());
    //         ba.refine(levelset__refinement);

    //         level_sets[1].reset(new MultiFab);
    //         level_sets[1]->define(ba, part_dm, 1, levelset__pad);
    //         iMultiFab valid_ref(ba, part_dm, 1, levelset__pad);

    //         // NOTE: implicit function data might not be on the right grids
    //         MultiFab impfunc = MFUtil::regrid(ba, part_dm, * implicit_functions[1], true);

    //         LSFactory::fill_data(* level_sets[1], valid_ref, * particle_ebfactory[0], impfunc,
    //                              32, levelset__refinement, 1, geom_lev, geom[0]);
    //     }
    // }
    // else
    // {

    //     const DistributionMapping & part_dm = pc->ParticleDistributionMap(0);
    //     const BoxArray &            part_ba = pc->ParticleBoxArray(0);

    //     //_______________________________________________________________________
    //     // Multi-level level-set: build finer level using coarse level set

    //     EBFArrayBoxFactory eb_factory(* eb_levels[0], geom[0], part_ba, part_dm,
    //                                   {levelset__eb_pad + 2, levelset__eb_pad + 2,
    //                                    levelset__eb_pad + 2}, EBSupport::full);

    //     // NOTE: reference BoxArray is not nodal
    //     BoxArray ba = amrex::convert(part_ba, IntVect::TheNodeVector());
    //     level_sets[0].reset(new MultiFab);
    //     level_sets[0]->define(ba, part_dm, 1, levelset__pad);
    //     iMultiFab valid(ba, part_dm, 1, levelset__pad);

    //     // NOTE: implicit function data might not be on the right grids
    //     MultiFab impfunc = MFUtil::regrid(ba, part_dm, * implicit_functions[0], true);

    //     LSFactory::fill_data(* level_sets[0], valid, * particle_ebfactory[0], impfunc,
    //                          32, 1, 1, geom[0], geom[0]);

    //     for (int lev = 1; lev < nlev; lev++)
    //     {

    //         const DistributionMapping & part_dm = pc->ParticleDistributionMap(lev);
    //         const BoxArray &            part_ba = pc->ParticleBoxArray(lev);

    //         // NOTE: reference BoxArray is not nodal
    //         BoxArray ba = amrex::convert(part_ba, IntVect::TheNodeVector());
    //         level_sets[lev].reset(new MultiFab);
    //         iMultiFab valid(ba, part_dm, 1, levelset__pad);

    //         // Fills level-set[lev] with coarse data
    //         LSCoreBase::MakeNewLevelFromCoarse( * level_sets[lev], * level_sets[lev-1],
    //                                            part_ba, part_dm, geom[lev], geom[lev-1],
    //                                            bcs_ls, refRatio(lev-1));

    //         EBFArrayBoxFactory eb_factory(* eb_levels[lev], geom[lev], part_ba, part_dm,
    //                                       {levelset__eb_pad + 2, levelset__eb_pad + 2,
    //                                        levelset__eb_pad + 2}, EBSupport::full);

    //         // NOTE: implicit function data might not be on the right grids
    //         MultiFab impfunc = MFUtil::regrid(ba, part_dm, * implicit_functions[lev]);

    //         IntVect ebt_size{AMREX_D_DECL(32, 32, 32)}; // Fudge factors...
    //         LSCoreBase::FillLevelSet(* level_sets[lev], * level_sets[lev], eb_factory, impfunc,
    //                                  ebt_size, levelset__eb_pad, geom[lev]);
    //     }
    // }

    // Add walls (for instance MI) to levelset data
    intersect_ls_walls();
}


void mfix::intersect_ls_walls ()
{

    bool has_walls = false;
    std::unique_ptr<UnionListIF<EB2::PlaneIF>> walls = get_walls(has_walls);
    auto gshop = EB2::makeShop(* walls);

    if (has_walls == false)
        return;

    if (nlev == 1)
    {
        //_______________________________________________________________________
        // Baseline Level-Set
        {
            const int ng = level_sets[0]->nGrow();
            const BoxArray & ba = level_sets[0]->boxArray();
            const DistributionMapping & dm = level_sets[0]->DistributionMap();

            MultiFab wall_if(ba, dm, 1, ng);
            iMultiFab valid(ba, dm, 1, ng);
            valid.setVal(1);

            GShopLSFactory<UnionListIF<EB2::PlaneIF>> gshop_lsf(gshop, geom[0], ba, dm, ng);
            std::unique_ptr<MultiFab> impfunc = gshop_lsf.fill_impfunc();

            LSFactory::fill_data(wall_if, valid, *impfunc, levelset__eb_pad, geom[0]);
            LSFactory::intersect_data(*level_sets[0], valid, wall_if, valid, geom[0]);
        }

        //_______________________________________________________________________
        // Refined Level-Set
        // TODO: Don't actually refine this thing if levelset refinement is 1
        {
            const int ng = level_sets[1]->nGrow();
            const BoxArray & ba = level_sets[1]->boxArray();
            const DistributionMapping & dm = level_sets[1]->DistributionMap();

            MultiFab wall_if(ba, dm, 1, ng);
            iMultiFab valid(ba, dm, 1, ng);
            valid.setVal(1);

            GShopLSFactory<UnionListIF<EB2::PlaneIF>> gshop_lsf(gshop, geom[0], ba, dm, ng);
            std::unique_ptr<MultiFab> impfunc = gshop_lsf.fill_impfunc();

            LSFactory::fill_data(wall_if, valid, *impfunc, levelset__eb_pad, geom[0]);
            LSFactory::intersect_data(*level_sets[1], valid, wall_if, valid, geom[0]);
        }
    }
    else
    {
        //_______________________________________________________________________
        // Multi-level level-set: apply wall-intersection to each level

        for (int lev = 0; lev < nlev; lev++)
        {
            const int ng = level_sets[lev]->nGrow();
            const BoxArray & ba = level_sets[lev]->boxArray();
            const DistributionMapping & dm = level_sets[lev]->DistributionMap();

            MultiFab wall_if(ba, dm, 1, ng);
            iMultiFab valid(ba, dm, 1, ng);
            valid.setVal(1);

            GShopLSFactory<UnionListIF<EB2::PlaneIF>> gshop_lsf(gshop, geom[lev], ba, dm, ng);
            std::unique_ptr<MultiFab> impfunc = gshop_lsf.fill_impfunc();

            LSFactory::fill_data(wall_if, valid, *impfunc, levelset__eb_pad, geom[lev]);
            LSFactory::intersect_data(*level_sets[lev], valid, wall_if, valid, geom[lev]);
        }
    }


}
