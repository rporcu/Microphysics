#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_ParmParse.H>
#include <MFIX_BC_Parms.H>

namespace BC
{

  // Direction of pressure drop (0:x, 1:y, 2:z)
  int delp_dir = -1;

  // Specified constant gas density
  amrex::Real delp[3];


  void Initialize ()
  {

    // Initialize the periodic pressure drop to zero in all directions
    for (int dir=0; dir < 3; dir ++ )
      delp[dir] = 0.0;

    amrex::ParmParse pp("bc");

    amrex::Real delp_in(0.0);
    pp.query("delp", delp_in);

    if( delp_in != 0.0) {

      pp.query("delp_dir", delp_dir);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE( 0 <= delp_dir && delp_dir <= 2,
         "Direction of pressure drop (bc.delp_dir) must be specified.");

      delp[delp_dir] = delp_in;

    }

  }

}
