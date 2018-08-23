#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

//#include <AMReX_GeometryShop.H>
//#include <AMReX_SphereIF.H>
//#include <AMReX_PlaneIF.H>
//#include <AMReX_AllRegularService.H>
//#include <AMReX_FlatPlateGeom.H>
//#include <AMReX_EBISLayout.H>
//#include <AMReX_EBGraph.H>
//#include <AMReX_EBDebugOut.H>
//#include <AMReX_EBCellFAB.H>
//#include <AMReX_EBCellFactory.H>
//#include <AMReX_EBIndexSpace.H>
//#include <AMReX_UnionIF.H>
//#include <AMReX_TransformIF.H>
//#include <AMReX_ComplementIF.H>
//#include <AMReX_IntersectionIF.H>
//#include <AMReX_LatheIF.H>
//#include <AMReX_PolynomialIF.H>
//#include <AMReX_AnisotropicDxPlaneIF.H>
//#include <AMReX_AnisotropicIF.H>

//#include <AMReX_VisMF.H>  // amrex::VisMF::Write(MultiFab)
//#include <sstream>

#include <algorithm>
#include <AMReX_EB_levelset.H>
#include <mfix_level.H>
#include <mfix_eb_F.H>

std::unique_ptr<BaseIF> mfix_level::get_walls(int lev, bool anisotropic, bool & has_walls)
{
    // Extracts all walls from the mfix.dat

    has_walls = false;  // will be set to true if there are any walls

    // Walls can be defined per phase => Itterarte over all phases and check
    // each for walls in the mfix.dat
    Vector<std::unique_ptr<BaseIF>> planes;
    for (int i = 1; i <= 500; i++) {
        int exists;
        RealVect normal, center;
        mfix_get_walls(& i, & exists, & normal, & center);
        if(exists) {
            has_walls = true;
            amrex::Print() << "Normal " << normal << std::endl;
            amrex::Print() << "Center " << center << std::endl;
            std::unique_ptr<PlaneIF> plane;
            if(anisotropic) {
                RealVect dxVec;
                for(int idir = 0; idir < 3; idir++)
                    dxVec[idir] = geom[lev].CellSize()[idir];
                plane = std::unique_ptr<PlaneIF>(new AnisotropicDxPlaneIF(normal, center, true, dxVec));
            } else {
                plane = std::unique_ptr<PlaneIF>(new PlaneIF(normal, center, true));
            }
            planes.push_back(std::move(plane));
        }
    }

    // The IntersectionIF constructor requires a vector of pointers to the
    // individual PlaneIFs => Construct a pointer-vector. Note that these
    // pointers become invalid once the `planes` vector goes out of scope.
    Vector<BaseIF *> plane_ptrs;
    for(std::unique_ptr<BaseIF> & pl : planes)
        plane_ptrs.push_back(pl.get());
    IntersectionIF all_planes(plane_ptrs);

    return std::unique_ptr<BaseIF>(all_planes.newImplicitFunction());
}

std::unique_ptr<BaseIF> mfix_level::get_real_walls(int lev, bool anisotropic, bool & has_real_walls)
{
    // Extracts all walls from the mfix.dat

    has_real_walls = false;  // will be set to true if there are any walls

    // Walls can be defined per phase => Itterarte over all phases and check
    // each for walls in the mfix.dat
    Vector<std::unique_ptr<BaseIF>> planes;
    for (int i = 1; i <= 500; i++) {
        int exists;
        RealVect normal, center;
        mfix_get_real_walls(& i, & exists, & normal, & center);
        if(exists) {
            has_real_walls = true;
            amrex::Print() << "Normal " << normal << std::endl;
            amrex::Print() << "Center " << center << std::endl;
            std::unique_ptr<PlaneIF> plane;
            if(anisotropic) {
                RealVect dxVec;
                for(int idir = 0; idir < 3; idir++)
                    dxVec[idir] = geom[lev].CellSize()[idir];
                plane = std::unique_ptr<PlaneIF>(new AnisotropicDxPlaneIF(normal, center, true, dxVec));
            } else {
                plane = std::unique_ptr<PlaneIF>(new PlaneIF(normal, center, true));
            }
            planes.push_back(std::move(plane));
        }
    }

    // The IntersectionIF constructor requires a vector of pointers to the
    // individual PlaneIFs => Construct a pointer-vector. Note that these
    // pointers become invalid once the `planes` vector goes out of scope.
    Vector<BaseIF *> plane_ptrs;
    for(std::unique_ptr<BaseIF> & pl : planes)
        plane_ptrs.push_back(pl.get());
    IntersectionIF all_planes(plane_ptrs);

    return std::unique_ptr<BaseIF>(all_planes.newImplicitFunction());
}
