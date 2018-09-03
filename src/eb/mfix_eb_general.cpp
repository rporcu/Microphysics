#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Polynomial.H>
#include <AMReX_EB2_IF_Translation.H>

#include <mfix_eb_if.H>

//#include <AMReX_VisMF.H>  // amrex::VisMF::Write(MultiFab)
//#include <sstream>

#include <algorithm>
#include <AMReX_EB_levelset.H>
#include <mfix_level.H>
#include <mfix_eb_F.H>


void
mfix_level::make_eb_general(int lev) {

    ParmParse pp("mfix");

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
    std::unique_ptr<UnionListIF<EB2::PlaneIF>>        impfunc_walls_part;
    std::unique_ptr<UnionListIF<EB2::PlaneIF>>        impfunc_walls_fluid;

    // Planar dividing wall
    std::unique_ptr<EB2::IntersectionIF<EB2::PlaneIF,EB2::PlaneIF,EB2::PlaneIF>> divider;

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

        divider = make_wall(div_dir, div_pos, div_height, div_width);
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

