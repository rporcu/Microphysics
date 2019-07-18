#include <fstream>
#include <iomanip>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_iMultiFab.H>

#include <AMReX_buildInfo.H>

#include <mfix.H>
#include <mfix_F.H>

int   max_step    = -1;
int   regrid_int  = -1;
Real stop_time    = -1.0;

bool write_eb_surface = false;
bool write_ls         = false;

std::string restart_file {""};

int repl_x = 1;
int repl_y = 1;
int repl_z = 1;

int  check_int = -1;
int  last_chk  = -1;
std::string check_file {"chk"};

int   plot_int = -1;
int   last_plt = -1;
std::string plot_file {"plt"};
std::string static_plt_file {"plt_ls"};

bool plotfile_on_restart = false;

int par_ascii_int = -1;
int last_par_ascii  = -1;
std::string par_ascii_file {"par"};

int avg_int = -1;
int last_avg  = -1;
std::string avg_file {"avg_region"};

std::string mfix_dat {"mfix.dat"};

void set_ptr_to_mfix(mfix& my_mfix);

void ReadParameters ()
{
  {
     ParmParse pp("amr");

     pp.query("stop_time", stop_time);
     pp.query("max_step", max_step);

     pp.query("check_file", check_file);
     pp.query("check_int", check_int);

     pp.query("plot_file", plot_file);
     pp.query("plot_int", plot_int);

     pp.query("plotfile_on_restart", plotfile_on_restart);

     pp.query("avg_int", avg_int );
     pp.query("avg_file", avg_file);

     pp.query("par_ascii_file", par_ascii_file);
     pp.query("par_ascii_int", par_ascii_int);

     pp.query("restart", restart_file);

     pp.query("repl_x", repl_x);
     pp.query("repl_y", repl_y);
     pp.query("repl_z", repl_z);
     pp.query("regrid_int",regrid_int);

     if ( regrid_int == 0 )
       amrex::Abort("regrid_int must be > 0 or < 0");
  }

  {
     ParmParse pp("mfix");

     pp.query("input_deck", mfix_dat);
     pp.query("write_eb_surface", write_eb_surface);
     pp.query("write_ls", write_ls);
  }
}

int main (int argc, char* argv[])
{
    // Issue an error if AMR input file is not given
    if ( argc < 2 )
       amrex::Abort("AMReX input file missing");

    // AMReX will now read the inputs file and the command line arguments, but the
    //        command line arguments are in mfix-format so it will just ignore them.
    amrex::Initialize(argc,argv);
    { // This start bracket and the end bracket before Finalize are essential so
      // that the mfix object is deleted before Finalize

    BL_PROFILE_VAR("main()", pmain)
    BL_PROFILE_REGION_START("mfix::main()");

    // Write out the MFIX git hash (the AMReX git hash is already written)
    const char* githash_mfix = buildInfoGetGitHash(1);
    amrex::Print() << "MFiX git hash: " << githash_mfix<< "\n";

    amrex::Gpu::setLaunchRegion(false);

    // Setting format to NATIVE rather than default of NATIVE_32
    FArrayBox::setFormat(FABio::FAB_NATIVE);

    // Copy arguments into MFIX -- note that the first argument is now the name of the
    //      inputs file to be read by AMReX, so we only pass the arguments after that
    for(int i=2; i < argc; i++) {
       int nlen = strlen(argv[i]);

       // If-statement avoids passing the name of the mfix input file if it is
       // specified on the command line or any AMReX command.
       if ( (strstr(argv[i], "input_file") == NULL) && (strstr(argv[i], "amr") == NULL)
                                                    && (strstr(argv[i], "mfix") == NULL) )
         mfix_add_argument(argv[i], &nlen);
    }

    Real strt_time = ParallelDescriptor::second();

    ReadParameters();

    int solve_fluid;
    int solve_dem;
    int call_udf;
    Real time=0.0L;
    int nstep = 0;  // Current time step

    Real dt = -1.;

    const char *cmfix_dat = mfix_dat.c_str();
    int name_len=mfix_dat.length();

    // Loads parameters (data) from fortran backend. Most notably this
    // subroutine loads the parameters from the `mfix.dat` file:
    //     mfix_get_data -> get_data -> read_namelist
    //                                        |
    //      (loads `mfix.dat`) ---------------+
    mfix_get_data( &solve_fluid, &solve_dem, &call_udf, &name_len, cmfix_dat);

    // Default constructor. Note inheritance: mfix : AmrCore : AmrMesh
    //                                                             |
    //  => Geometry is constructed here: (constructs Geometry) ----+
    mfix my_mfix;
    solve_fluid=1;
    my_mfix.get_input_bcs();

    if ( ParallelDescriptor::IOProcessor() )
      check_inputs();

    // Set global static pointer to mfix object. Used by fill-patch utility
    set_ptr_to_mfix(my_mfix);

    // Initialize internals from ParamParse database
    my_mfix.InitParams(solve_fluid, solve_dem, call_udf);

    // Initialize memory for data-array internals
    my_mfix.ResizeArrays();

    // Initialize EB geometry. This needs to be done before grid creation (in
    // mfix::Init), as the grids are created using each EB-level's volfrac.
    my_mfix.make_eb_geometry();

    // Initialize derived internals
    my_mfix.Init(time);

    // Create EB factories on new grids
    my_mfix.make_eb_factories();

    my_mfix.InitLevelData(time);
    my_mfix.calcVol();

    } // This end bracket and the start bracket after Initialize are essential so
      // that the mfix object is deleted before Finalize

    amrex::Finalize();
    return 0;
}
