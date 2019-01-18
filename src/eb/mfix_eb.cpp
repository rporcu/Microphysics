#include <algorithm>
#include <AMReX_EB_levelset.H>
#include <mfix.H>
#include <mfix_eb_F.H>

void mfix::make_eb_geometry ()
{

    /****************************************************************************
     * mfix.geometry=<string> specifies the EB geometry. <string> can be on of  *
     * box, cylinder, hopper, clr, clr_riser, general (or blank)                *
     ***************************************************************************/

    ParmParse pp("mfix");

    std::string geom_type;
    pp.query("geometry", geom_type);

    /****************************************************************************
     * Legacy inputs:                                                           *
     *   -- mfix.hourglass = true <=> mfix.geometry=box                         *
     *   -- mfix.clr       = true <=> mfix.geometry=clr                         *
     *   -- mfix.clr_riser = true <=> mfix.geometry=clr_riser                   *
     *   -- mfix.use_walls = true <=> mfix.geometry=general                     *
     *   -- mfix.use_poy2  = true <=> mfix.geometry=general                     *
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

    // // store copy of level set (for plotting and velocity reconstruction).
    // if (solve_dem && !levelset__restart) {
    //     std::unique_ptr<MultiFab> ls_data = level_set->coarsen_data();
    //     int ng = ls_data->nGrow();

    //     // TODO: currently this is a wee bit of a hack... as each level is
    //     // getting the same level-set function. Once the AMR levelset is fully
    //     // implemented, use that instead.
    //     for (int lev = 0; lev < nlev; lev++)
    //         ls[lev]->copy(* ls_data, 0, 0, 1, ng, ng);
    // }
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


void mfix::fill_eb_levelsets () {
    if (nlev == 1) {

        const DistributionMapping & part_dm = pc->ParticleDistributionMap(0);
        const BoxArray &            part_ba = pc->ParticleBoxArray(0);

        //_______________________________________________________________________
        // Baseline Level-Set
        {
            // NOTE: reference BoxArray is not nodal
            BoxArray ba = amrex::convert(part_ba, IntVect::TheNodeVector());

            level_sets[0].reset(new MultiFab);
            level_sets[0]->define(ba, part_dm, 1, levelset__pad);
            iMultiFab valid(ba, part_dm, 1, levelset__pad);

            // NOTE: implicit function data might not be on the right grids
            MultiFab impfunc = MFUtil::regrid(ba, part_dm, * implicit_functions[0]);

            LSFactory::fill_data(* level_sets[0], valid, * particle_ebfactory[0], impfunc,
                                 32, 1, 1, geom[0], geom[0]);
        }

        //_______________________________________________________________________
        // Refined Level-Set
        {
            // Set up refined geometry
            Box dom = geom[0].Domain();
            dom.refine(levelset__refinement);
            Geometry geom_lev(dom);

            // Set up refined BoxArray. NOTE: reference BoxArray is not nodal
            BoxArray ba = amrex::convert(part_ba, IntVect::TheNodeVector());
            ba.refine(levelset__refinement);

            level_sets[1].reset(new MultiFab);
            level_sets[1]->define(ba, part_dm, 1, levelset__pad);
            iMultiFab valid_ref(ba, part_dm, 1, levelset__pad);

            // NOTE: implicit function data might not be on the right grids
            MultiFab impfunc = MFUtil::regrid(ba, part_dm, * implicit_functions[1]);

            LSFactory::fill_data(* level_sets[1], valid_ref, * particle_ebfactory[0], impfunc,
                                 32, levelset__refinement, 1, geom_lev, geom[0]);

            VisMF::Write(*level_sets[1], "level_set");
        }
    } else {

        const DistributionMapping & part_dm = pc->ParticleDistributionMap(0);
        const BoxArray &            part_ba = pc->ParticleBoxArray(0);

        //_______________________________________________________________________
        // Multi-level level-set: build finer level using coarse level set

        EBFArrayBoxFactory eb_factory(* eb_levels[0], geom[0], part_ba, part_dm,
                                      {levelset__eb_pad + 1, levelset__eb_pad + 1, levelset__eb_pad + 1},
                                      EBSupport::full);

        // NOTE: reference BoxArray is not nodal
        BoxArray ba = amrex::convert(part_ba, IntVect::TheNodeVector());
        level_sets[0].reset(new MultiFab);
        level_sets[0]->define(ba, part_dm, 1, levelset__pad);
        iMultiFab valid(ba, part_dm, 1, levelset__pad);

        // NOTE: implicit function data might not be on the right grids
        MultiFab impfunc = MFUtil::duplicate<MultiFab,
                                             MFUtil::SymmetricGhost>(ba, part_dm, * implicit_functions[0]);

        LSFactory::fill_data(* level_sets[0], valid, * particle_ebfactory[0], impfunc,
                             32, 1, 1, geom[0], geom[0]);

        for (int lev = 1; lev < nlev; lev++) {

            const DistributionMapping & part_dm = pc->ParticleDistributionMap(lev);
            const BoxArray &            part_ba = pc->ParticleBoxArray(lev);

            // NOTE: reference BoxArray is not nodal
            BoxArray ba = amrex::convert(part_ba, IntVect::TheNodeVector());
            level_sets[lev].reset(new MultiFab);
            iMultiFab valid(ba, part_dm, 1, levelset__pad);

            // Fills level-set[lev] with coarse data
            LSCoreBase::MakeNewLevelFromCoarse( * level_sets[lev], * level_sets[lev-1],
                                               part_ba, part_dm, geom[lev], geom[lev-1],
                                               bcs_ls, refRatio(lev-1));

            EBFArrayBoxFactory eb_factory(* eb_levels[lev], geom[lev], part_ba, part_dm,
                                          {levelset__eb_pad + 1, levelset__eb_pad + 1, levelset__eb_pad + 1},
                                          EBSupport::full);

            // NOTE: implicit function data might not be on the right grids
            MultiFab impfunc = MFUtil::duplicate<MultiFab,
                                             MFUtil::SymmetricGhost>(ba, part_dm, * implicit_functions[lev]);

            IntVect ebt_size{AMREX_D_DECL(32, 32, 32)}; // Fudge factors...
            LSCoreBase::FillLevelSet(* level_sets[lev], * level_sets[lev], eb_factory, impfunc,
                                     ebt_size, levelset__eb_pad, geom[lev]);
        }
    }
}


void mfix::make_amr_geometry ()
{
    if (! use_amr_ls)
        return;


    /****************************************************************************
     * mfix.geometry=<string> specifies the EB geometry. <string> can be on of  *
     * box, cylinder, hopper, clr, clr_riser, general (or blank)                *
     ***************************************************************************/

    ParmParse pp("mfix");

    std::string geom_type;
    pp.query("geometry", geom_type);

    /****************************************************************************
     * Legacy inputs:                                                           *
     *   -- mfix.hourglass = true <=> mfix.geometry=box                         *
     *   -- mfix.clr       = true <=> mfix.geometry=clr                         *
     *   -- mfix.clr_riser = true <=> mfix.geometry=clr_riser                   *
     *   -- mfix.use_walls = true <=> mfix.geometry=general                     *
     *   -- mfix.use_poy2  = true <=> mfix.geometry=general                     *
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
        amrex::Print() << "\n Building box AMR level-set." << std::endl;
        make_amr_box();
    } else if (geom_type == "cylinder") {
        amrex::Print() << "\n Building cylinder AMR level-set." << std::endl;
        make_amr_cylinder();
    } else if (geom_type == "hopper") {
        amrex::Print() << "\n Building hopper AMR level-set." << std::endl;
        make_amr_hopper();
    } else if (geom_type == "cyclone") {
        amrex::Print() << "\n Building cyclone AMR level-set." << std::endl;
        make_amr_cyclone();
    } else if(geom_type == "general") {
        amrex::Print() << "\n Building general AMR level-set (poly2 with extra walls)." << std::endl;
        make_amr_general();
    } else {
        amrex::Print() << "\n No EB geometry declared in inputs => "
                       << " Will read walls from mfix.dat only."
                       << std::endl;
        make_amr_regular();
    }


    if (use_amr_ls)
    {
        amrex::Print() << "Constructing AMR levelset on "
                       << amr_ls_max_level + 1 << " levels" << std::endl;

        int ct_lev = 0;
        Vector<Real> ls_tags;
        for (Real tag = amr_ls_baseline_tag; tag > 0; tag -= amr_ls_tag_step){
            ls_tags.push_back(tag);
            ct_lev ++;

            if (ct_lev >= amr_ls_max_level + 1)
                break;
        }

        // amrex::Print() << "Level tags: ";
        // for (Real tag : ls_tags)
        //     amrex::Print() << tag << " ";
        // amrex::Print() << std::endl;

        //amr_level_set->InitData(ls_tags);
        amr_level_set->InitData(false); // Don't use levelset tagging, instead use volfrac
        amr_level_set->WritePlotFile();
        amrex::Print() << "... done constructing AMR levelset" << std::endl;
    }
}
