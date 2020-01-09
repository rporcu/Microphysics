#include <algorithm>
#include <AMReX_EB_levelset.H>
#include <mfix.H>
#include <mfix_F.H>
#include <mfix_eb_F.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Intersection.H>



void mfix::make_eb_geometry ()
{

    /****************************************************************************
     *                                                                          *
     * IMPORTANT: the MFIX BC types need to be set in order to correctly        *
     *            identify walls (which are needed to build the right EBs)      *
     *                                                                          *
     ***************************************************************************/

    MakeBCArrays();

    int cyc_x=0, cyc_y=0, cyc_z=0;

    if (geom[0].isPeriodic(0)) cyc_x = 1;
    if (geom[0].isPeriodic(1)) cyc_y = 1;
    if (geom[0].isPeriodic(2)) cyc_z = 1;

    mfix_set_cyclic(&cyc_x, &cyc_y, &cyc_z);

    for (int lev = 0; lev < nlev; lev++)
        mfix_set_bc_type(lev);

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
     *   -- mfix.use_walls = true <=> mfix.geometry=general                     *
     *   -- mfix.use_poy2  = true <=> mfix.geometry=general                     *
     *                                                                          *
     ***************************************************************************/


    bool hourglass    = false;
    bool eb_general   = false;

    pp.query("hourglass", hourglass);

    bool eb_poly2 = false;
    bool eb_walls = false;

    pp.query("use_poly2", eb_poly2);
    pp.query("use_walls", eb_walls);
    eb_general = eb_poly2 || eb_walls;

    // Avoid multiple (ambiguous) inputs
    if (hourglass || eb_general) {
        if (! geom_type.empty()) {
            amrex::Abort("The input file cannot specify both:\n"
                         "mfix.<geom_type>=true and mfix.geometry=<geom_type>\n"
                         "at the same time."                                   );
        }
    }

    if (hourglass)  geom_type = "hourglass";
    if (eb_general) geom_type = "general";


    /****************************************************************************
     *                                                                          *
     *  CONSTRUCT EB                                                            *
     *                                                                          *
     ***************************************************************************/

    contains_ebs = false; // default

    if (geom_type == "box") {
        amrex::Print() << "\n Building box geometry." << std::endl;
        make_eb_box();
        contains_ebs = true;
    } else if (geom_type == "cylinder") {
        amrex::Print() << "\n Building cylinder geometry." << std::endl;
        make_eb_cylinder();
        contains_ebs = true;
    } else if (geom_type == "hopper") {
        amrex::Print() << "\n Building hopper geometry." << std::endl;
        make_eb_hopper();
        contains_ebs = true;
    } else if (geom_type == "cyclone") {
        amrex::Print() << "\n Building cyclone geometry." << std::endl;
        make_eb_cyclone();
        contains_ebs = true;
    } else if(geom_type == "air-reactor") {
      amrex::Print() << "\n Building air-reactor geometry." << std::endl;
        make_eb_air_reactor();
        contains_ebs = true;
    } else if(geom_type == "prototype clr") {
      amrex::Print() << "\n Building full-loop clr." << std::endl;
      make_eb_proto_clr();
      contains_ebs = true;
    } else if(geom_type == "general") {
      amrex::Print() << "\n Building general geometry (poly2 with extra walls)." << std::endl;
      // TODO: deal with inflow volfrac
      make_eb_general();
      contains_ebs = true;
    } else {
        amrex::Print() << "\n No EB geometry declared in inputs => "
                       << " Will read walls from mfix.dat only."
                       << std::endl;
        make_eb_regular();

        // NOTE: `mfix::contains_ebs` this is set internally depending if EBs
        // where detected in the mfix.dat file.
    }
}


void mfix::make_eb_factories () {

    /****************************************************************************
     *                                                                          *
     * Fill EB factories as an initial run. Since the particle container might  *
     * not have been created yet, use mfix grids instead.                       *
     *                                                                          *
     ***************************************************************************/

    for (int lev = 0; lev < nlev; lev++)
    {
        ebfactory[lev].reset(
            new EBFArrayBoxFactory(* eb_levels[lev], geom[lev], grids[lev], dmap[lev],
                                   {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                    m_eb_full_grow_cells}, m_eb_support_level)
            );

        // Grow EB factory by +2 in order to avoid edge cases. This is not
        // necessary for multi-level mfix.
        particle_ebfactory[lev].reset(
            new EBFArrayBoxFactory(* particle_eb_levels[lev], geom[lev], grids[lev], dmap[lev],
                                   {levelset__eb_pad + 2, levelset__eb_pad + 2,
                                    levelset__eb_pad + 2}, m_eb_support_level)
            );
    }
}


void mfix::fill_eb_levelsets ()
{
    if (ooo_debug) amrex::Print() << "fill_eb_levelsets" << std::endl;
    /****************************************************************************
     *                                                                          *
     * Fill levels either as a single level with refinement, or as a proper     *
     * multi-level hierarchy                                                    *
     *                                                                          *
     ***************************************************************************/

    if (nlev == 1)
    {

        //_______________________________________________________________________
        // COMPATIBILITY: in order to be compatible with benchmarks, fill using
        // the level-set factory (which is then thrown away).

        const DistributionMapping & part_dm = pc->ParticleDistributionMap(0);
        const BoxArray &            part_ba = pc->ParticleBoxArray(0);

        LSFactory lsf(0, levelset__refinement, levelset__eb_refinement,
                      levelset__pad, levelset__eb_pad, part_ba, geom[0], part_dm );


        //___________________________________________________________________________
        // NOTE: Boxes are different (since we're not refining, we need to treat
        // corners this way). IMPORTANT: the box case is assembled from planes
        // => fill the level-set with these (otherwise corners will be
        // inaccurately resolved.)

        ParmParse pp("mfix");

        std::string geom_type;
        pp.query("geometry", geom_type);


        if (geom_type == "box")
        {
            ParmParse pp("box");

            if ( geom[0].isAllPeriodic() )
            {
                make_eb_regular();
            }
            else
            {
                Vector<Real> boxLo(3), boxHi(3);
                Real offset    = 1.0e-15;

                for (int i = 0; i < 3; i++)
                {
                    boxLo[i] = geom[0].ProbLo(i);
                    boxHi[i] = geom[0].ProbHi(i);
                }

                pp.queryarr("Lo", boxLo,  0, 3);
                pp.queryarr("Hi", boxHi,  0, 3);

                pp.query("offset", offset);

                Real xlo = boxLo[0] + offset;
                Real xhi = boxHi[0] - offset;

                // This ensures that the walls won't even touch the ghost cells.
                // By putting them one domain width away
                if (geom[0].isPeriodic(0))
                {
                    xlo = 2.0*geom[0].ProbLo(0) - geom[0].ProbHi(0);
                    xhi = 2.0*geom[0].ProbHi(0) - geom[0].ProbLo(0);
                }

                Real ylo = boxLo[1] + offset;
                Real yhi = boxHi[1] - offset;

                // This ensures that the walls won't even touch the ghost cells.
                // By putting them one domain width away
                if (geom[0].isPeriodic(1))
                {
                    ylo = 2.0*geom[0].ProbLo(1) - geom[0].ProbHi(1);
                    yhi = 2.0*geom[0].ProbHi(1) - geom[0].ProbLo(1);
                }

                Real zlo = boxLo[2] + offset;
                Real zhi = boxHi[2] - offset;

                // This ensures that the walls won't even touch the ghost cells.
                // By putting them one domain width away
                if (geom[0].isPeriodic(2))
                {
                    zlo = 2.0*geom[0].ProbLo(2) - geom[0].ProbHi(2);
                    zhi = 2.0*geom[0].ProbHi(2) - geom[0].ProbLo(2);
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

                auto if_box = EB2::makeUnion(plane_lox, plane_hix, plane_loy,
                                             plane_hiy, plane_loz, plane_hiz );
                auto gshop = EB2::makeShop(if_box);

                GShopLSFactory<decltype(if_box)> gshop_lsfactory(gshop, lsf);
                std::unique_ptr<MultiFab> mf_impfunc_box = gshop_lsfactory.fill_impfunc();
                lsf.Intersect(* mf_impfunc_box);
            }

            level_sets[1] = lsf.copy_data(part_dm);
            level_sets[0] = lsf.coarsen_data();

            return;
        }


        //___________________________________________________________________________
        // If there is an eb in the system, fill everything else using the
        // LSFactory EBF filling routine

        if (contains_ebs)
        {
            MultiFab impfunc(lsf.get_ls_ba(), lsf.get_dm(), 1, lsf.get_ls_pad());
            eb_levels[1]->fillLevelSet(impfunc, lsf.get_ls_geom());
            impfunc.FillBoundary(lsf.get_ls_geom().periodicity());

            lsf.Fill( * ebfactory[0], impfunc);
        }

        level_sets[1] = lsf.copy_data(part_dm);
        level_sets[0] = lsf.coarsen_data();
    }
    else
    {

        const DistributionMapping & part_dm = pc->ParticleDistributionMap(0);
        const BoxArray &            part_ba = pc->ParticleBoxArray(0);

        //_______________________________________________________________________
        // Multi-level level-set: build finer level using coarse level set

        EBFArrayBoxFactory eb_factory(* eb_levels[0], geom[0], part_ba, part_dm,
                                      {levelset__eb_pad + 2, levelset__eb_pad + 2,
                                       levelset__eb_pad + 2}, EBSupport::full);

        // NOTE: reference BoxArray is not nodal
        BoxArray ba = amrex::convert(part_ba, IntVect::TheNodeVector());
        level_sets[0].reset(new MultiFab);
        level_sets[0]->define(ba, part_dm, 1, levelset__pad);
        iMultiFab valid(ba, part_dm, 1, levelset__pad);

        MultiFab impfunc(ba, part_dm, 1, levelset__pad);
        eb_levels[0]->fillLevelSet(impfunc, geom[0]);
        impfunc.FillBoundary(geom[0].periodicity());


        LSFactory::fill_data(* level_sets[0], valid, * particle_ebfactory[0], impfunc,
                             32, 1, 1, geom[0], geom[0]);

        for (int lev = 1; lev < nlev; lev++)
        {

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
                                          {levelset__eb_pad + 2, levelset__eb_pad + 2,
                                           levelset__eb_pad + 2}, EBSupport::full);

            MultiFab impfunc(ba, part_dm, 1, levelset__pad);
            eb_levels[lev]->fillLevelSet(impfunc, geom[lev]);
            impfunc.FillBoundary(geom[lev].periodicity());

            IntVect ebt_size{AMREX_D_DECL(32, 32, 32)}; // Fudge factors...
            LSCoreBase::FillLevelSet(* level_sets[lev], * level_sets[lev], eb_factory, impfunc,
                                     ebt_size, levelset__eb_pad, geom[lev]);
        }
    }

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

            const auto & flags = particle_ebfactory[0]->getMultiEBCellFlagFab();
            int ngrow_eb = flags.nGrow();
            LSFactory::fill_data(wall_if, valid, *impfunc, ngrow_eb, geom[0]);
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

            // Set up refined geometry
            Box dom = geom[0].Domain();
            dom.refine(levelset__refinement);
            Geometry geom_lev(dom);

            GShopLSFactory<UnionListIF<EB2::PlaneIF>> gshop_lsf(gshop, geom_lev, ba, dm, ng);
            std::unique_ptr<MultiFab> impfunc = gshop_lsf.fill_impfunc();

            const auto & flags = particle_ebfactory[0]->getMultiEBCellFlagFab();
            int ngrow_eb = flags.nGrow();
            LSFactory::fill_data(wall_if, valid, *impfunc, ngrow_eb, geom_lev);
            LSFactory::intersect_data(*level_sets[1], valid, wall_if, valid, geom_lev);
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
