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

    // store copy of level set (for plotting and velocity reconstruction).
    if (solve_dem && !levelset__restart) {
        std::unique_ptr<MultiFab> ls_data = level_set->coarsen_data();
        int ng = ls_data->nGrow();

        // TODO: currently this is a wee bit of a hack... as each level is
        // getting the same level-set function. Once the AMR levelset is fully
        // implemented, use that instead.
        for (int lev = 0; lev < nlev; lev++)
            ls[lev]->copy(* ls_data, 0, 0, 1, ng, ng);
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
                       << std::endl;
        make_amr_regular();
    }


    if (use_amr_ls)
    {
        amrex::Print() << "Constructing AMR levelset on "
                       << amr_ls_max_level + 1 << " levels" << std::endl;

        Vector<Real> ls_tags;
        for (Real tag = amr_ls_baseline_tag; tag > 0; tag -= amr_ls_tag_step)
            ls_tags.push_back(tag);

        amrex::Print() << "Level tags: ";
        for (Real tag : ls_tags)
            amrex::Print() << tag << " ";
        amrex::Print() << std::endl;

        amr_level_set->InitData(ls_tags);
        amr_level_set->WritePlotFile();
        amrex::Print() << "... done constructing AMR levelset" << std::endl;
    }
}

