#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Polynomial.H>
#include <AMReX_EB2_IF_Translation.H>

#include <mfix_eb_if.H>

//#include <AMReX_VisMF.H>  // amrex::VisMF::Write(MultiFab)
//#include <sstream>

#include <algorithm>
#include <type_traits>
#include <AMReX_EB_levelset.H>
#include <mfix_level.H>
#include <mfix_eb_F.H>


void
mfix_level::make_eb_general(int lev) {

    ParmParse pp("mfix");

    int max_level_here = 0;

    bool use_walls   = true;
    bool use_poly2   = false;
    bool use_divider = false;

    pp.query("use_walls",   use_walls);
    pp.query("use_poly2",   use_poly2);
    pp.query("use_divider", use_divider);


    /****************************************************************************
     * IF using divider => Extract parameters                                   *
     ****************************************************************************/

    int div_dir = 0;
    Real div_pos, div_height, div_width;

    if(use_divider) {
        pp.query("divider_position", div_pos);
        pp.query("divider_height",   div_height);
        pp.query("divider_width",    div_width);
        pp.query("div_dir",          div_dir);
    }

    // Non-planar side wall
    std::unique_ptr<EB2::TranslationIF<EB2::PolynomialIF>> impfunc_poly2;

    // Planar side walls (particles can see additional walls)
    std::unique_ptr<UnionListIF<EB2::PlaneIF>> impfunc_walls_part;
    std::unique_ptr<UnionListIF<EB2::PlaneIF>> impfunc_walls_fluid;

    // Planar dividing wall
    std::unique_ptr<
        EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF
                            >> impfunc_divider;

    /****************************************************************************
     * Generate PolynomialIF representing the non-planar EB walls               *
     ****************************************************************************/
    if (use_poly2) {
        amrex::Print() << "Using poly2 geometry" << std::endl;

        impfunc_poly2 = get_poly(lev, SpaceDim, "poly2");
    }

    /****************************************************************************
     * Generate UnionListIF representing the planar EB walls                    *
     *          this IF represents the union of a list of planes (walls)        *
     ****************************************************************************/
    // Flags checking if mfix.dat even has walls
    // IMPORTANT NOTE: has_real_walls => has_walls <=> ! has_walls => ! has_real_walls
    bool has_walls = false, has_real_walls = false;
    if (use_walls) {
        amrex::Print() << "Using wall geometry from mfix.dat" << std::endl;

        impfunc_walls_part  = get_walls(lev, has_walls);
        impfunc_walls_fluid = get_real_walls(lev, has_real_walls);
    }

    /****************************************************************************
     * Generate IntersectionIF representing a planar (partial) partition        *
     *          this partition has finite thickness and penetrates the "bottom" *
     *          wall => it is constructed using 3 PlaneIFs                      *
     ****************************************************************************/
    if (use_divider) {
        amrex::Print() << "Using divider-wall geometry" << std::endl;

        impfunc_divider = make_wall(div_dir, div_pos, div_height, div_width);
    }


    /***************************************************************************
     *                                                                         *
     * Build EB Factories                                                      *
     *                                                                         *
     ***************************************************************************/

    // Stores implicit function for the combined particle IF
    std::unique_ptr<MultiFab> mf_impfunc;
    // Stores implicit function for the particle walls IF only
    std::unique_ptr<MultiFab> mf_impfunc_walls;
    // Stores implicit function representing the polynomial "walls"
    std::unique_ptr<MultiFab> mf_impfunc_poly2;
    // For DEM: save the polynomial level separately (to allow "water-tight"
    // intersection with walls).
    const EB2::Level * poly_lev;

    EBSupport m_eb_support_level = EBSupport::full;

    int max_coarsening_level = 100;

    if (use_poly2) {

        // For DEM: generate polynomial IF separately (to allow "water-tight"
        // intersection with walls).
        if(solve_dem){
            auto gshop = EB2::makeShop(* impfunc_poly2);
            EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
            const EB2::IndexSpace & poly_ebis = EB2::IndexSpace::top();
            poly_lev = & poly_ebis.getLevel(geom[lev]);

            GShopLSFactory<std::decay<decltype(* impfunc_poly2)>::type
                           > ls_gshop_poly(gshop, * level_set);
            mf_impfunc_poly2 = ls_gshop_poly.fill_impfunc();
        }

        if (has_walls && use_divider) { // ............................... poly2 + walls + divider

            if (solve_dem) {
                amrex::Print() << "Making the particle ebfactory ..." << std::endl;

                auto eb_if = EB2::makeUnion(* impfunc_poly2, * impfunc_walls_part, * impfunc_divider);
                auto gshop = EB2::makeShop(eb_if);

                EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
                const EB2::IndexSpace & eb_is = EB2::IndexSpace::top();
                eb_level_particles = & eb_is.getLevel(geom.back());

                GShopLSFactory<decltype(eb_if)> gshop_lsfactory(gshop, * level_set);
                mf_impfunc = gshop_lsfactory.fill_impfunc();

                auto eb_if_walls = EB2::makeUnion(* impfunc_walls_part, * impfunc_divider);
                auto gshop_walls = EB2::makeShop(eb_if_walls);

                GShopLSFactory<decltype(eb_if_walls)> walls_lsfactory(gshop_walls, * level_set);
                mf_impfunc_walls = walls_lsfactory.fill_impfunc();

                amrex::Print() << "Done making the particle ebfactory." << std::endl;
            }

            if (solve_fluid) {
                amrex::Print() << "Making the fluid ebfactory ..." << std::endl;

                if (has_real_walls) { // since ! has_walls => ! has_real_walls
                    auto gshop = EB2::makeShop(EB2::makeUnion(* impfunc_poly2,
                                                              * impfunc_walls_fluid,
                                                              * impfunc_divider)
                                               );
                    EB2::Build(gshop, geom.back(), max_level_here, max_level_here +
                               max_coarsening_level);
                } else {
                    auto gshop = EB2::makeShop(EB2::makeUnion(* impfunc_poly2, * impfunc_divider));
                    EB2::Build(gshop, geom.back(), max_level_here, max_level_here +
                               max_coarsening_level);
                }
                const EB2::IndexSpace & eb_is = EB2::IndexSpace::top();
                eb_level_fluid = & eb_is.getLevel(geom.back());

                amrex::Print() << "Done making the fluid ebfactory." << std::endl;
            }

        } else if (has_walls) { // ....................................... poly2 + walls + ! divider

            if (solve_dem) {
                amrex::Print() << "Making the particle ebfactory ..." << std::endl;

                auto eb_if = EB2::makeUnion(* impfunc_poly2, * impfunc_walls_part);
                auto gshop = EB2::makeShop(eb_if);

                EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
                const EB2::IndexSpace & eb_is = EB2::IndexSpace::top();
                eb_level_particles = & eb_is.getLevel(geom.back());

                GShopLSFactory<decltype(eb_if)> gshop_lsfactory(gshop, * level_set);
                mf_impfunc = gshop_lsfactory.fill_impfunc();

                auto gshop_walls = EB2::makeShop(* impfunc_walls_part);

                GShopLSFactory<std::decay<decltype(* impfunc_walls_part)>::type
                               > walls_lsfactory(gshop_walls, * level_set);
                mf_impfunc_walls = walls_lsfactory.fill_impfunc();

                amrex::Print() << "Done making the particle ebfactory." << std::endl;
            }

            if (solve_fluid) {
                amrex::Print() << "Making the fluid ebfactory ..." << std::endl;

                if (has_real_walls) { // since ! has_walls => ! has_real_walls
                    auto gshop = EB2::makeShop(EB2::makeUnion(* impfunc_poly2,
                                                              * impfunc_walls_fluid)
                                               );
                    EB2::Build(gshop, geom.back(), max_level_here, max_level_here +
                               max_coarsening_level);
                } else {
                    auto gshop = EB2::makeShop(* impfunc_poly2);
                    EB2::Build(gshop, geom.back(), max_level_here, max_level_here +
                               max_coarsening_level);
                }
                const EB2::IndexSpace & eb_is = EB2::IndexSpace::top();
                eb_level_fluid = & eb_is.getLevel(geom.back());

                amrex::Print() << "Done making the fluid ebfactory." << std::endl;
            }

        } else if (use_divider) { // ..................................... poly2 + ! walls + divider

            // NOTE: the divider applies to both, particles _and_ fluid

            amrex::Print() << "Making the particle and fluid ebfactory ..." << std::endl;

            auto eb_if = EB2::makeUnion(* impfunc_poly2, * impfunc_divider);
            auto gshop = EB2::makeShop(eb_if);

            EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
            const EB2::IndexSpace & eb_is = EB2::IndexSpace::top();
            eb_level_particles = & eb_is.getLevel(geom.back());
            eb_level_fluid     = eb_level_particles;

            if (solve_dem) {
                GShopLSFactory<decltype(eb_if)> gshop_lsfactory(gshop, * level_set);
                mf_impfunc = gshop_lsfactory.fill_impfunc();

                auto gshop_walls = EB2::makeShop(* impfunc_divider);

                GShopLSFactory<std::decay<decltype(* impfunc_divider)>::type
                               > walls_lsfactory(gshop_walls, * level_set);
                mf_impfunc_walls = walls_lsfactory.fill_impfunc();
            }

            amrex::Print() << "Done making the particle and fluid ebfactory." << std::endl;

        } else { // ...................................................... poly2 + ! walls + ! divider

            amrex::Print() << "Making the particle and fluid ebfactory ..." << std::endl;

            auto gshop = EB2::makeShop(* impfunc_poly2);

            EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
            const EB2::IndexSpace & eb_is = EB2::IndexSpace::top();
            eb_level_particles = & eb_is.getLevel(geom.back());
            eb_level_fluid     = eb_level_particles;

            if (solve_dem) {
                GShopLSFactory<std::decay<decltype(* impfunc_poly2)>::type
                               > gshop_lsfactory(gshop, * level_set);
                mf_impfunc = gshop_lsfactory.fill_impfunc();
            }

            amrex::Print() << "Done making the particle and fluid ebfactory." << std::endl;
        }

    } else {
        if (has_walls && use_divider) { // ............................... ! poly2 + walls + divider

            if (solve_dem) {
                amrex::Print() << "Making the particle ebfactory ..." << std::endl;

                auto eb_if = EB2::makeUnion(* impfunc_walls_part, * impfunc_divider);
                auto gshop = EB2::makeShop(eb_if);

                EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
                const EB2::IndexSpace & eb_is = EB2::IndexSpace::top();
                eb_level_particles = & eb_is.getLevel(geom.back());

                GShopLSFactory<decltype(eb_if)> gshop_lsfactory(gshop, * level_set);
                mf_impfunc       = gshop_lsfactory.fill_impfunc();
                mf_impfunc_walls = gshop_lsfactory.fill_impfunc(); // without poly2 IF walls == IF

                amrex::Print() << "Done making the particle ebfactory." << std::endl;
            }

            if (solve_fluid) {
                amrex::Print() << "Making the fluid ebfactory ..." << std::endl;

                if (has_real_walls) { // since ! has_walls => ! has_real_walls
                    auto gshop = EB2::makeShop(EB2::makeUnion(* impfunc_walls_fluid,
                                                              * impfunc_divider)
                                               );
                    EB2::Build(gshop, geom.back(), max_level_here, max_level_here +
                               max_coarsening_level);
                } else {
                    auto gshop = EB2::makeShop(* impfunc_divider);
                    EB2::Build(gshop, geom.back(), max_level_here, max_level_here +
                               max_coarsening_level);
                }
                const EB2::IndexSpace & eb_is = EB2::IndexSpace::top();
                eb_level_fluid = & eb_is.getLevel(geom.back());

                amrex::Print() << "Done making the fluid ebfactory." << std::endl;
            }

        } else if (has_walls) { // ....................................... ! poly2 + walls + ! divider

            if (solve_dem) {
                auto gshop = EB2::makeShop(* impfunc_walls_part);

                EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
                const EB2::IndexSpace & eb_is = EB2::IndexSpace::top();
                eb_level_particles = & eb_is.getLevel(geom.back());

                GShopLSFactory<std::decay<decltype(* impfunc_walls_part)>::type
                               > gshop_lsfactory(gshop, * level_set);
                mf_impfunc = gshop_lsfactory.fill_impfunc();
                mf_impfunc_walls = gshop_lsfactory.fill_impfunc(); // without poly2 IF walls == IF
            }

            if (solve_fluid) {
                amrex::Print() << "Making the fluid ebfactory ..." << std::endl;

                if (has_real_walls) { // since ! has_walls => ! has_real_walls
                    auto gshop = EB2::makeShop(* impfunc_walls_fluid);
                    EB2::Build(gshop, geom.back(), max_level_here, max_level_here +
                               max_coarsening_level);
                } else {
                    EB2::AllRegularIF my_regular;
                    auto gshop = EB2::makeShop(my_regular);
                    EB2::Build(gshop, geom.back(), max_level_here, max_level_here +
                               max_coarsening_level);
                }

                const EB2::IndexSpace & eb_is = EB2::IndexSpace::top();
                eb_level_fluid = & eb_is.getLevel(geom.back());

                amrex::Print() << "Done making the fluid ebfactory." << std::endl;
            }

        } else if (use_divider) { // ..................................... ! poly2 + ! walls + divider

            amrex::Print() << "Making the particle and fluid ebfactory ..." << std::endl;

            auto gshop = EB2::makeShop(* impfunc_divider);

            EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
            const EB2::IndexSpace & eb_is = EB2::IndexSpace::top();
            eb_level_particles = & eb_is.getLevel(geom.back());
            eb_level_fluid     = eb_level_particles;

            if (solve_dem){
                GShopLSFactory<std::decay<decltype(* impfunc_divider)>::type
                               > gshop_lsfactory(gshop, * level_set);
                mf_impfunc = gshop_lsfactory.fill_impfunc();
                mf_impfunc_walls = gshop_lsfactory.fill_impfunc(); // without poly2 IF walls == IF
            }

            amrex::Print() << "Done making the particle and fluid ebfactory." << std::endl;

        } else { // .................................................... ! poly2 + ! walls + ! divider
            // Do nothing... (this will never happen)
        }

    }

    if (solve_fluid)
        ebfactory[lev].reset(new EBFArrayBoxFactory(* eb_level_fluid,
                                                    geom[lev], grids[lev], dmap[lev],
                                                    {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                                     m_eb_full_grow_cells}, m_eb_support_level));

    if (solve_dem) {
        particle_ebfactory[lev].reset(new EBFArrayBoxFactory(* eb_level_particles,
                                                             geom[lev], grids[lev], dmap[lev],
                                                             {m_eb_basic_grow_cells,
                                                              m_eb_volume_grow_cells,
                                                              m_eb_full_grow_cells},
                                                             m_eb_support_level));

        eb_normals = pc->EBNormals(lev, particle_ebfactory[lev].get(), dummy.get());

        /************************************************************************
         *                                                                      *
         * Fill level-set:                                                      *
         *                                                                      *
         ************************************************************************/
        if (!levelset__restart) {
            amrex::Print() << "Creating the levelset ..." << std::endl;

            if (use_walls){
                level_set->intersection_impfunc(* mf_impfunc_walls);
            }

            if (use_poly2) {
                int eb_pad = level_set->get_eb_pad();
                Geometry geom_eb = LSUtility::make_eb_geometry(* level_set, geom[lev]);
                EBFArrayBoxFactory eb_factory_poly(* poly_lev, geom_eb,
                                                   level_set->get_eb_ba(),
                                                   level_set->get_dm(),
                                                   {eb_pad, eb_pad, eb_pad},
                                                   EBSupport::full);
                level_set->intersection_ebf(eb_factory_poly, * mf_impfunc_poly2);
            }

            // store copy of level set (for plotting).
            std::unique_ptr<MultiFab> ls_data = level_set->coarsen_data();
            ls[lev]->copy(* ls_data, 0, 0, 1, 0, 0);
            ls[lev]->FillBoundary(geom[lev].periodicity());

            amrex::Print() << "Done making the levelset ..." << std::endl;
        } else {
            amrex::Print() << "Loaded level-set is fine => skipping levelset calculation."
                           << std::endl;
        }
    }
}


std::unique_ptr<EB2::TranslationIF<EB2::PolynomialIF>>
mfix_level::get_poly(int lev, int max_order, std::string field_prefix) {

    /****************************************************************************
     * Read polynomial data from inputs database                                *
     *      => Generate Vector<EB2::PolyTerm> used                              *
     ****************************************************************************/

    // Construct the ParamParse database field names based on the
    // `field_prefix` string:
    ParmParse pp("mfix");

    // Coefficients vector is stored in the inputs database with the field name:
    //      <field_prefix>_[x,y,z]_coeffs
    const std::array<const string, 3> var_names{"x", "y", "z"};
          std::array<      string, 3> field_names;
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
    mirror_field    << field_prefix << "_mirror";
    translate_field << field_prefix << "_translate";

    // Generate vector representing polynomial
    Vector<EB2::PolyTerm> poly;
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

            EB2::PolyTerm mono = {.coef = coef, .powers = powers};
            poly.push_back(mono);
        }
    }

    /***************************************************************************
     * Construct PolynomialIF (and apply translation)                          *
     ***************************************************************************/

    bool flip = false;
    pp.query(mirror_field.str().c_str(), flip);

    EB2::PolynomialIF mirror(poly, flip);

    Vector<Real> offset_vec(SpaceDim); RealArray offset;
    pp.getarr(translate_field.str().c_str(), offset_vec, 0, SpaceDim);
    std::copy_n(offset_vec.begin(), SpaceDim, offset.begin());

    std::unique_ptr<EB2::TranslationIF<EB2::PolynomialIF>> ret =
        std::unique_ptr<EB2::TranslationIF<EB2::PolynomialIF>>(
            new EB2::TranslationIF<EB2::PolynomialIF>(mirror, offset)
        );

    return ret;
}


std::unique_ptr<EB2::IntersectionIF<EB2::PlaneIF,EB2::PlaneIF,EB2::PlaneIF>>
mfix_level::make_wall( int dir, // direction (long edge) of wall
                       Real position, Real height, Real width )
{
    RealArray normal, center;

    // Upward-facing plane:
    normal = {0.,   0.,   1.};
    center = {0.,   0., height};
    EB2::PlaneIF plane_up(center, normal, false);

    // First side-facing plane
    normal = {0., 0., 0.};
    center = {0., 0., 0.};
    normal[dir] = 1.;
    center[dir] = position + width;
    EB2::PlaneIF plane_1(center, normal, false);

    // Second side-facing plane
    normal = {0., 0., 0.};
    center = {0., 0., 0.};
    normal[dir] = -1.;
    center[dir] = position - width;
    EB2::PlaneIF plane_2(center, normal, false);

    std::unique_ptr<EB2::IntersectionIF<EB2::PlaneIF,EB2::PlaneIF,EB2::PlaneIF>> ret =
        std::unique_ptr<EB2::IntersectionIF<EB2::PlaneIF,EB2::PlaneIF,EB2::PlaneIF>>
        (new EB2::IntersectionIF<EB2::PlaneIF,EB2::PlaneIF,EB2::PlaneIF>(plane_up, plane_1, plane_2));

    return ret;
}
