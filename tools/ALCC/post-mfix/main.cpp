//
// Average a variable with optional density weighting (Favre average)

#include <regex>
#include <algorithm>

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <post_mfix.H>
#include <pm_fluid.H>

using namespace amrex;

std::string inputs_name = "";

void PrintHelp ();


int main(int argc, char* argv[])
{

  amrex::Initialize(argc, argv, false);

  {
    // timer for profiling
    BL_PROFILE_VAR("main()", pmain);

    const int narg = amrex::command_argument_count();

    if( narg == 0) {
      PrintHelp();
      return  0;
    }

    post_mfix post(amrex::get_command_argument(narg));

    if ( post.error() ) {
      PrintHelp();
      amrex::Abort("Error detected!");
    }

    post.print_vars();
    post.report();

    // destroy timer for profiling
    BL_PROFILE_VAR_STOP(pmain);

  }
  amrex::Finalize();
}



void PrintHelp ()
{
  Print() << "\n"
          << " Compute fluid avarges for variable in a plotfile\n"
          << " usage: \n"
          << "    post-mfix [OPTION]... plotfile\n"
          << "\n"
          << "   -f, --fluid    : comma separated list of fluid variables to average\n"
          << "   -p, --particle : comma separated list of particle variables to average\n"
          << "   --print-vars   : list all fluid and particle variables in plotfile\n"
          << "\n"
          << "\n"
          << "\n"
          << " Example: \n"
          << " ./post-mfix3d.gnu.ex -f u_g,v_g,w_g -p velx,vely,velz plt43668\n"
          << std::endl;
}
