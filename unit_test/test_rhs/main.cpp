#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_BCRec.H>


using namespace amrex;

#include <test_react.H>
#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <rhs_zones.H>
#include <AMReX_buildInfo.H>
#include <variables.H>
#include <unit_test.H>
#include <react_util.H>

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        main_main();
    }
    amrex::Finalize();
    return 0;
}

void main_main ()
{

    // AMREX_SPACEDIM: number of dimensions
    int n_cell, max_grid_size;

    std::string prefix = "plt";

#ifdef AMREX_USE_GPU
    IntVect tile_size(1024, 1024, 1024);
#else
    IntVect tile_size(1024, 8, 8);
#endif

    // inputs parameters
    {
        // ParmParse is way of reading inputs from the inputs file
        ParmParse pp;

        // We need to get n_cell from the inputs file - this is the
        // number of cells on each side of a square (or cubic) domain.
        pp.get("n_cell", n_cell);

        // The domain is broken into boxes of size max_grid_size
        max_grid_size = 32;
        pp.query("max_grid_size", max_grid_size);

        pp.query("prefix", prefix);

    }

    Vector<int> is_periodic(AMREX_SPACEDIM,0);
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
      is_periodic[idim] = 1;
    }

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
        IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);

        // Break up boxarray "ba" into chunks no larger than
        // "max_grid_size" along a direction
        ba.maxSize(max_grid_size);

        // This defines the physical box, [0, 1] in each direction.
        RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                         {AMREX_D_DECL(1.0, 1.0, 1.0)});

        // This defines a Geometry object
        geom.define(domain, &real_box,
                    CoordSys::cartesian, is_periodic.data());
    }

    // Nghost = number of ghost cells for each array
    int Nghost = 0;

    // do the runtime parameter initializations and microphysics inits
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "reading extern runtime parameters ..." << std::endl;
    }

    ParmParse ppa("amr");

    init_unit_test();

    if (integrator_rp::jacobian != 1) {
        amrex::Error("test_rhs only works for integrator.jacobian = 1");
    }

    // C++ EOS initialization (must be done after init_extern_parameters)
    eos_init(small_temp, small_dens);

    // C++ Network, RHS, screening, rates initialization
    network_init();

    plot_t vars;
    vars = init_variables();
    amrex::Vector<std::string> names;
    get_varnames(vars, names);

    // time = starting time in the simulation
    Real time = 0.0;

    // How Boxes are distributed among MPI processes
    DistributionMapping dm(ba);

    // we allocate our main multifabs
    MultiFab state(ba, dm, vars.n_plot_comps, Nghost);

    // Initialize the state

    Real dlogrho;
    Real dlogT;

    if (n_cell > 1) {
        dlogrho = (std::log10(dens_max) - std::log10(dens_min))/(n_cell - 1);
        dlogT   = (std::log10(temp_max) - std::log10(temp_min))/(n_cell - 1);
    } else {
        dlogrho = 0.0_rt;
        dlogT   = 0.0_rt;
    }

    init_t comp_data = setup_composition(n_cell);

    for ( MFIter mfi(state); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        auto state_arr = state.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {

            state_arr(i, j, k, vars.itemp) =
                std::pow(10.0_rt, (std::log10(temp_min) + static_cast<Real>(j)*dlogT));
            state_arr(i, j, k, vars.irho) =
                std::pow(10.0_rt, (std::log10(dens_min) + static_cast<Real>(i)*dlogrho));

            Real xn[NumSpec];
            get_xn(k, comp_data, xn, uniform_xn);

            for (int n = 0; n < NumSpec; n++) {
                state_arr(i, j, k, vars.ispec_old+n) =
                    amrex::max(xn[n], 1.e-10_rt);
            }
        });

    }

    // What time is it now?  We'll use this to compute total react time.
    Real strt_time = ParallelDescriptor::second();

    // Call the rhs
#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(state, tile_size); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.tilebox();

        auto s = state.array(mfi);

        AMREX_PARALLEL_FOR_3D(bx, i, j, k,
        {
            do_rhs(i, j, k, s, vars);
        });

    }

    // Call the timer again and compute the maximum difference between
    // the start time and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time, IOProc);

    // get the name of the integrator from the build info functions
    // written at compile time.  We will append the name of the
    // integrator to the output file name
    int nmodules = buildInfoGetNumModules();

    int int_idx = -1;
    std::string INT_KEY = "INTEGRATOR";
    for (int i=1; i<=nmodules; i++) {
      if (INT_KEY == buildInfoGetModuleName(i)) {
        int_idx = i;
      }
    }

    std::string name = "test_rhs.";
    std::string integrator = buildInfoGetModuleVal(int_idx);

    // Write a plotfile

    WriteSingleLevelPlotfile(prefix + name + integrator, state, names, geom, time, 0);

    write_job_info(prefix + name + integrator);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;

}
