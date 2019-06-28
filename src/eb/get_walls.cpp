#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

#include <mfix_eb_if.H>

#include <algorithm>
#include <mfix.H>
#include <mfix_F.H>
#include <mfix_eb_F.H>


std::unique_ptr<UnionListIF<EB2::PlaneIF>>
mfix::get_walls(bool & has_walls) {
    // Extracts all walls from the mfix.dat

    has_walls = false;  // will be set to true if there are any walls

    // Walls can be defined per phase => Itterarte over all phases and check
    // each for walls in the mfix.dat
    Vector<EB2::PlaneIF> planes;
    for (int i = 1; i <= 6; i++) {
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
mfix::get_real_walls(bool & has_real_walls) {
    // Extracts all walls from the mfix.dat

    has_real_walls = false;  // will be set to true if there are any walls

    // Walls can be defined per phase => Itterarte over all phases and check
    // each for walls in the mfix.dat
    Vector<EB2::PlaneIF> planes;
    for (int i = 1; i <= 6; i++) {
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

void
mfix::get_input_bcs(){

  // Extracts all walls from the inputs file
  int cyclic;

  cyclic = geom[0].isPeriodic(0) ? 1 : 0;
  set_input_bcs("xlo", 1, cyclic, geom[0].ProbLo(0));
  set_input_bcs("xhi", 2, cyclic, geom[0].ProbHi(0));

  cyclic = geom[0].isPeriodic(1) ? 1 : 0;
  set_input_bcs("ylo", 3, cyclic, geom[0].ProbLo(1));
  set_input_bcs("yhi", 4, cyclic, geom[0].ProbHi(1));

  cyclic = geom[0].isPeriodic(2) ? 1 : 0;
  set_input_bcs("zlo", 5, cyclic, geom[0].ProbLo(2));
  set_input_bcs("zhi", 6, cyclic, geom[0].ProbHi(2));

}

void
mfix::set_input_bcs(const std::string bcID, const int index,
                    const int cyclic, const Real domloc) {

  const int und_  =   0;
  const int ig_   =   9;
  const int pinf_ =  10;
  const int pout_ =  11;
  const int minf_ =  20;
  const int nsw_  = 100;

  // Default a BC to ignore.
  int itype = ig_;

  Real pressure = -1.0;
  Real velocity =  0.0;
  Real location = domloc;
  Vector<Real> velocity_in;

  std::string bc_type = "null";

  ParmParse pp(bcID);

  pp.query("type", bc_type);

  if (bc_type == "null"){

    // Unspecified domain extents are assumed to be walls flush with
    itype = (cyclic == 1) ? und_ : ig_;

    if(cyclic == 0)
      amrex::Print() << bcID <<" set to ignore. "  << std::endl;

  } else if(bc_type == "pressure_inflow"  || bc_type == "pi" ||
            bc_type == "PRESSURE_INFLOW"  || bc_type == "PI" ) {

    amrex::Print() << bcID <<" set to pressure inflow. "  << std::endl;
    itype = pinf_;

    pp.get("pressure", pressure);

  } else if(bc_type == "pressure_outflow" || bc_type == "po" ||
            bc_type == "PRESSURE_OUTFLOW" || bc_type == "PO" ) {

    amrex::Print() << bcID <<" set to pressure outflow. "  << std::endl;
    itype = pout_;

    pp.get("pressure", pressure);


  } else if (bc_type == "mass_inflow"     || bc_type == "mi" ||
             bc_type == "MASS_INFLOW"     || bc_type == "MI" ) {

    // Flag that this is a mass inflow.
    amrex::Print() << bcID <<" set to mass inflow. "  << std::endl;
    itype = minf_;

    pp.query("pressure", pressure);
    pp.queryarr("velocity", velocity_in );


  } else if (bc_type == "no_slip_wall"    || bc_type == "nsw" ||
             bc_type == "NO_SLIP_WALL"    || bc_type == "NSW" ) {

    // Flag that this is a no-slip wall.
    amrex::Print() << bcID <<" set to no-slip wall. "  << std::endl;
    itype = nsw_;

    pp.query("location", location);

  }

  if ( cyclic == 1 && itype != und_){
    amrex::Abort("Cannot mix periodic BCs and Wall/Flow BCs.\n");
  }


  // Treat the first velocity as-usual.
  int velocity_size = velocity_in.size();
  if(velocity_size > 0) velocity = velocity_in[0];


  const Real* plo = geom[0].ProbLo();
  const Real* phi = geom[0].ProbHi();

  mfix_set_bc_mod(&index, &itype, plo, phi,
                  &location, &pressure, &velocity);

  // If more than one velocity is provided for a mass inflow, then there
  // sub-regions have been defined to overlap. These are appended to the
  // BC list.

  if(velocity_in.size() > 1){

    Vector<Real> xlo, xhi, ylo, yhi, zlo, zhi;

    if( index == 1 || index == 2){

      pp.queryarr("ylo", ylo );
      pp.queryarr("yhi", yhi );
      pp.queryarr("zlo", zlo );
      pp.queryarr("zhi", zhi );

      // Copy in values for x.
      Real val = (index == 1) ? plo[0] : phi[0];
      for (int bc = 0; bc < velocity_size - 1; bc++) {
        xlo.push_back(val);
        xhi.push_back(val);
      }

    } else if( index == 3 || index == 4) {

      pp.queryarr("xlo", xlo );
      pp.queryarr("xhi", xhi );
      pp.queryarr("zlo", zlo );
      pp.queryarr("zhi", zhi );

      // Copy in values for y.
      Real val = (index == 3) ? plo[1] : phi[1];
      for (int bc = 0; bc < velocity_size - 1; bc++) {
        ylo.push_back(val);
        yhi.push_back(val);
      }

    } else if( index == 5 || index == 6) {

      pp.queryarr("xlo", xlo );
      pp.queryarr("xhi", xhi );
      pp.queryarr("ylo", ylo );
      pp.queryarr("yhi", yhi );

      // Copy in values for z.
      Real val = (index == 5) ? plo[2] : phi[2];
      for (int bc = 0; bc < velocity_size - 1; bc++) {
        zlo.push_back(val);
        zhi.push_back(val);
      }

    }

    // Make sure all the arrays are the correct size.
    AMREX_ALWAYS_ASSERT(xlo.size() == velocity_size-1);
    AMREX_ALWAYS_ASSERT(xhi.size() == velocity_size-1);
    AMREX_ALWAYS_ASSERT(ylo.size() == velocity_size-1);
    AMREX_ALWAYS_ASSERT(yhi.size() == velocity_size-1);
    AMREX_ALWAYS_ASSERT(zlo.size() == velocity_size-1);
    AMREX_ALWAYS_ASSERT(zhi.size() == velocity_size-1);

    for (int bc = 0; bc < velocity_size -1; bc++) {

      velocity = velocity_in[bc+1];

      mfix_set_bc_mod_add_mi(&index,
                      &xlo[bc], &ylo[bc], &zlo[bc],
                      &xhi[bc], &yhi[bc], &zhi[bc],
                      &pressure, &velocity);
    }
  }



}
