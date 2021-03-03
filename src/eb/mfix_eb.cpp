#include <algorithm>
#include <mfix.H>
#include <AMReX_EB2.H>
#include <AMReX_EB_utils.H>
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

    // nghost_tmp = 4 is enough to make the bc arrays for now; we will remake
    // them later when we know what nghost_state() really is
    const int nghost_tmp = 4;
    MakeBCArrays(nghost_tmp);

    for (int lev = 0; lev < nlev; lev++)
        mfix_set_bc_type(lev,nghost_tmp);

    /****************************************************************************
     *                                                                          *
     * mfix.geometry=<string> specifies the EB geometry. <string> can be on of  *
     * box, cylinder, hopper, clr, clr_riser, general (or blank)                *
     *                                                                          *
     ***************************************************************************/

    ParmParse pp("mfix");

    std::string geom_type;
    std::string csg_file;
    pp.query("geometry", geom_type);
    pp.query("geometry_filename", csg_file);
    amrex::Print() << "mfix.geometry_filename: " << csg_file;

#ifndef MFIX_GEOMETRY_CSG
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( csg_file.empty(), "CSG Geometry defined in input deck but solver not built with CSG support!");
#endif

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

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(geom_type.empty() || csg_file.empty(),
                                     "The input file cannot specify both:\n"
                                     "mfix.<geom_type> and mfix.geometry_filename\n"
                                     "at the same time.");

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

#ifdef MFIX_GEOMETRY_CSG
    } else if(!csg_file.empty()) {
      amrex::Print() << "\n Building geometry from .csg file:  " << csg_file << std::endl;
      make_eb_csg(csg_file);
      contains_ebs = true;
#endif

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
            new EBFArrayBoxFactory(*eb_levels[lev], geom[lev], grids[lev], dmap[lev],
                                   {nghost_eb_basic(), nghost_eb_volume(),
                                    nghost_eb_full()}, m_eb_support_level));

        // Grow EB factory by +2 in order to avoid edge cases. This is not
        // necessary for multi-level mfix.
        particle_ebfactory[lev].reset(
            new EBFArrayBoxFactory(*particle_eb_levels[lev], geom[lev], grids[lev], dmap[lev],
                                   {levelset_eb_pad + 2, levelset_eb_pad + 2,
                                    levelset_eb_pad + 2}, m_eb_support_level));
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

        const DistributionMapping& ls_dm = part_dm;
        const Geometry& ls_geom = amrex::refine(geom[0], levelset_refinement);

        BoxArray ls_ba = amrex::convert(part_ba, IntVect::TheNodeVector());
        level_sets[0].reset(new MultiFab(ls_ba, ls_dm, 1, levelset_pad/levelset_refinement));

        if (levelset_refinement != 1) ls_ba.refine(levelset_refinement);
        level_sets[1].reset(new MultiFab(ls_ba, ls_dm, 1, levelset_pad));

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
            if ( geom[0].isAllPeriodic() )
            {
                make_eb_regular();
                level_sets[1]->setVal(std::numeric_limits<Real>::max());
                level_sets[0]->setVal(std::numeric_limits<Real>::max());
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

                ParmParse pp_box("box");

                pp_box.queryarr("Lo", boxLo,  0, 3);
                pp_box.queryarr("Hi", boxHi,  0, 3);

                pp_box.query("offset", offset);

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

                amrex::FillImpFunc(*level_sets[1], gshop, ls_geom);
                level_sets[1]->negate(level_sets[1]->nGrow()); // signed distance f = - imp. f.
                if (levelset_refinement == 1) {
                    MultiFab::Copy(*level_sets[0], *level_sets[1], 0, 0, 1, level_sets[0]->nGrow());
                } else {
                    amrex::average_down_nodal(*level_sets[1], *level_sets[0],
                                              IntVect(levelset_refinement),
                                              level_sets[0]->nGrow(), true);
                }
            }

            return;
        }

        if (contains_ebs) {
            AMREX_ALWAYS_ASSERT(levelset_eb_refinement == 1); // Not sure about its purpose anyway
            amrex::FillSignedDistance(*level_sets[1], *eb_levels[1], *ebfactory[0],
                                      levelset_refinement);
            if (levelset_refinement == 1) {
                MultiFab::Copy(*level_sets[0], *level_sets[1], 0, 0, 1, level_sets[0]->nGrow());
            } else {
                amrex::average_down_nodal(*level_sets[1], *level_sets[0],
                                          IntVect(levelset_refinement),
                                          level_sets[0]->nGrow(), true);
            }
        } else {
            level_sets[1]->setVal(std::numeric_limits<Real>::max());
            level_sets[0]->setVal(std::numeric_limits<Real>::max());
        }
    }
    else
    {
        amrex::Abort("xxxxx fill_eb_levelsets todo");
#if 0
        const DistributionMapping & part_dm = pc->ParticleDistributionMap(0);
        const BoxArray &            part_ba = pc->ParticleBoxArray(0);

        //_______________________________________________________________________
        // Multi-level level-set: build finer level using coarse level set

        EBFArrayBoxFactory eb_factory(* eb_levels[0], geom[0], part_ba, part_dm,
                                      {levelset_eb_pad + 2, levelset_eb_pad + 2,
                                       levelset_eb_pad + 2}, EBSupport::full);

        // NOTE: reference BoxArray is not nodal
        BoxArray ba = amrex::convert(part_ba, IntVect::TheNodeVector());
        level_sets[0].reset(new MultiFab(ba, part_dm, 1, levelset_pad));
        iMultiFab valid(ba, part_dm, 1, levelset_pad);

        MultiFab impfunc(ba, part_dm, 1, levelset_pad);
        eb_levels[0]->fillLevelSet(impfunc, geom[0]);
        impfunc.FillBoundary(geom[0].periodicity());


        LSFactory::fill_data(*level_sets[0], valid, *particle_ebfactory[0], impfunc,
                             32, 1, 1, geom[0], geom[0]);

        for (int lev = 1; lev < nlev; lev++)
        {
            const DistributionMapping & part_dm_lev = pc->ParticleDistributionMap(lev);
            const BoxArray &            part_ba_lev = pc->ParticleBoxArray(lev);

            // NOTE: reference BoxArray is not nodal
            BoxArray ba_lev = amrex::convert(part_ba_lev, IntVect::TheNodeVector());
            if (level_sets[lev] != nullptr) delete level_sets[lev];
            level_sets[lev] = new MultiFab();
            // iMultiFab valid_lev(ba_lev, part_dm_lev, 1, levelset_pad);

            // Fills level-set[lev] with coarse data
            LSCoreBase::MakeNewLevelFromCoarse( * level_sets[lev], * level_sets[lev-1],
                                               part_ba_lev, part_dm_lev, geom[lev], geom[lev-1],
                                               bcs_ls, refRatio(lev-1));

            EBFArrayBoxFactory eb_factory_lev(* eb_levels[lev], geom[lev], part_ba_lev, part_dm_lev,
                                              {levelset_eb_pad + 2, levelset_eb_pad + 2,
                                               levelset_eb_pad + 2}, EBSupport::full);

            MultiFab impfunc_lev(ba_lev, part_dm_lev, 1, levelset_pad);
            eb_levels[lev]->fillLevelSet(impfunc_lev, geom[lev]);
            impfunc_lev.FillBoundary(geom[lev].periodicity());

            IntVect ebt_size{AMREX_D_DECL(32, 32, 32)}; // Fudge factors...
            LSCoreBase::FillLevelSet(* level_sets[lev], * level_sets[lev], eb_factory_lev, impfunc_lev,
                                     ebt_size, levelset_eb_pad, geom[lev]);
        }
#endif
    }

    // Add walls (for instance MI) to levelset data
    intersect_ls_walls();
}


void mfix::intersect_ls_walls ()
{

    bool has_walls = false;
    std::shared_ptr<UnionListIF<EB2::PlaneIF>> walls = get_walls(has_walls);
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

            amrex::FillImpFunc(wall_if, gshop, geom[0]);

            // The difference between this and the previous code is we no
            // longer clamp the level sets.
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*level_sets[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.growntilebox();
                Array4<Real> const& sdf = level_sets[0]->array(mfi);
                Array4<Real const> const& wif = wall_if.const_array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    sdf(i,j,k) = amrex::min(sdf(i,j,k),-wif(i,j,k));
                });
            }
        }

        //_______________________________________________________________________
        // Refined Level-Set
        if (levelset_refinement == 1) {
            MultiFab::Copy(*level_sets[1], *level_sets[0], 0, 0, 1, level_sets[0]->nGrow());
        } else {
            const int ng = level_sets[1]->nGrow();
            const BoxArray & ba = level_sets[1]->boxArray();
            const DistributionMapping & dm = level_sets[1]->DistributionMap();
            MultiFab wall_if(ba, dm, 1, ng);

            amrex::FillImpFunc(wall_if, gshop, amrex::refine(geom[0],levelset_refinement));

            // The difference between this and the previous code is we no
            // longer clamp the level sets.
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*level_sets[1],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.growntilebox();
                Array4<Real> const& sdf = level_sets[1]->array(mfi);
                Array4<Real const> const& wif = wall_if.const_array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    sdf(i,j,k) = amrex::min(sdf(i,j,k),-wif(i,j,k));
                });
            }
        }
    }
    else
    {
        amrex::Abort("xxxxx intersect_ls_walls todo");
#if 0
        //_______________________________________________________________________
        // Multi-level level-set: apply wall-intersection to each level

        for (int lev = 0; lev < nlev; lev++)
        {
            const int ng = level_sets[lev]->nGrow();
            const BoxArray & ba = level_sets[lev]->boxArray();
            const DistributionMapping & dm = level_sets[lev]->DistributionMap();

            MultiFab wall_if(ba, dm, 1, ng);
            iMultiFab valid_lev(ba, dm, 1, ng);
            valid_lev.setVal(1);

            GShopLSFactory<UnionListIF<EB2::PlaneIF>> gshop_lsf(gshop, geom[lev], ba, dm, ng);
            std::unique_ptr<MultiFab> impfunc = gshop_lsf.fill_impfunc();

            LSFactory::fill_data(wall_if, valid_lev, *impfunc, levelset_eb_pad, geom[lev]);
            LSFactory::intersect_data(*level_sets[lev], valid_lev, wall_if, valid_lev, geom[lev]);
        }
#endif
    }
}
