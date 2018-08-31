#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

#include <mfix_eb_if.H>

//#include <AMReX_VisMF.H>  // amrex::VisMF::Write(MultiFab)
//#include <sstream>

#include <algorithm>
// //#include <AMReX_EB_levelset.H>
//#include <AMReX_EB_levelset.H>
#include <mfix_level.H>
#include <mfix_eb_F.H>


std::unique_ptr<UnionListIF<EB2::PlaneIF>>
mfix_level::get_walls(int lev, bool & has_walls) {
    // Extracts all walls from the mfix.dat

    has_walls = false;  // will be set to true if there are any walls

    // Walls can be defined per phase => Itterarte over all phases and check
    // each for walls in the mfix.dat
    Vector<EB2::PlaneIF> planes;
    for (int i = 1; i <= 500; i++) {
        int exists;
        RealVect normal, center;
        mfix_get_walls(& i, & exists, & normal, & center);
        if(exists) {
            has_walls = true;
            amrex::Print() << "Normal " << normal << std::endl;
            amrex::Print() << "Center " << center << std::endl;

            RealArray p = {AMREX_D_DECL(center[0], center[1], center[2])};
            RealArray n = {AMREX_D_DECL(normal[0], normal[1], normal[2])};

            planes.emplace_back(p, n, false);
        }
    }

    std::unique_ptr<UnionListIF<EB2::PlaneIF>> ret =
        std::unique_ptr<UnionListIF<EB2::PlaneIF>>(new UnionListIF<EB2::PlaneIF>(planes));
    return ret;
}


std::unique_ptr<UnionListIF<EB2::PlaneIF>>
mfix_level::get_real_walls(int lev, bool & has_real_walls) {
    // Extracts all walls from the mfix.dat

    has_real_walls = false;  // will be set to true if there are any walls

    // Walls can be defined per phase => Itterarte over all phases and check
    // each for walls in the mfix.dat
    Vector<EB2::PlaneIF> planes;
    for (int i = 1; i <= 500; i++) {
        int exists;
        RealVect normal, center;
        mfix_get_real_walls(& i, & exists, & normal, & center);
        if(exists) {
            has_real_walls = true;
            amrex::Print() << "Normal " << normal << std::endl;
            amrex::Print() << "Center " << center << std::endl;

            RealArray p = {AMREX_D_DECL(center[0], center[1], center[2])};
            RealArray n = {AMREX_D_DECL(normal[0], normal[1], normal[2])};

            planes.emplace_back(p, n, false);
        }
    }

    std::unique_ptr<UnionListIF<EB2::PlaneIF>> ret =
        std::unique_ptr<UnionListIF<EB2::PlaneIF>>(new UnionListIF<EB2::PlaneIF>(planes));
    return ret;
}
