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
  amrex::Real delp = 0.0;


  void Initialize ()
  {

    amrex::ParmParse pp("bc");

    pp.query("delp", delp);

    if( delp != 0.0) {

      pp.query("delp_dir", delp_dir);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE( 0 <= delp_dir && delp_dir <= 2,
         "Direction of pressure drop (bc.delp_dir) must be specified.");

    }

  }

}
