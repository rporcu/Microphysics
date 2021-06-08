#ifndef MFIX_SOLIDS_PARMS_H_
#define MFIX_SOLIDS_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

#include <mfix_species_parms.H>

struct SolidParms {

  amrex::Real* h_cp_s0;
  amrex::Real* d_cp_s0;
  amrex::Real* h_cp_sn0;
  amrex::Real* d_cp_sn0;

  //AMREX_GPU_HOST_DEVICE
  //SolidParms () noexcept : cp_s0(nullptr), cp_sn0(nullptr) {}

  AMREX_GPU_HOST_DEVICE
  SolidParms (amrex::Real* h_cp_s,
              amrex::Real* d_cp_s,
              amrex::Real* h_cp_sn,
              amrex::Real* d_cp_sn)
    : h_cp_s0(h_cp_s)
    , d_cp_s0(d_cp_s)
    , h_cp_sn0(h_cp_sn)
    , d_cp_sn0(d_cp_sn)
  {}

#if defined(AMREX_USE_GPU)
  template <amrex::RunOn run_on>
#else
  template <amrex::RunOn run_on=amrex::RunOn::Host>
#endif
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real calc_cp_s (const int phase,
                         const amrex::Real /*T_s*/) const
  {
    if (run_on == amrex::RunOn::Gpu)
      return d_cp_s0[phase];
    else if (run_on == amrex::RunOn::Cpu)
      return h_cp_s0[phase];
    else
      amrex::Abort("Unrecongized RunOn option");

    return 0;
  }

#if defined(AMREX_USE_GPU)
  template <amrex::RunOn run_on>
#else
  template <amrex::RunOn run_on=amrex::RunOn::Host>
#endif
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real calc_partial_cp_s (const int /*phase*/,
                                 const amrex::Real /*T_s*/) const
  {
    return 0;
  }

#if defined(AMREX_USE_GPU)
  template <amrex::RunOn run_on>
#else
  template <amrex::RunOn run_on=amrex::RunOn::Host>
#endif
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real calc_h_s (const int phase,
                        const amrex::Real T_s) const
  {
    return calc_cp_s<run_on>(phase,T_s) * T_s;
  }

#if defined(AMREX_USE_GPU)
  template <amrex::RunOn run_on>
#else
  template <amrex::RunOn run_on=amrex::RunOn::Host>
#endif
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real calc_partial_h_s (const int phase,
                                const amrex::Real T_s) const
  {
    return calc_partial_cp_s<run_on>(phase,T_s) * T_s + calc_cp_s<run_on>(phase,T_s);
  }

#if defined(AMREX_USE_GPU)
  template <amrex::RunOn run_on>
#else
  template <amrex::RunOn run_on=amrex::RunOn::Host>
#endif
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real calc_cp_sn (const amrex::Real /*T_s*/,
                          const int n) const
  {
    if (run_on == amrex::RunOn::Gpu)
      return d_cp_sn0[n];
    else if (run_on == amrex::RunOn::Cpu)
      return h_cp_sn0[n];
    else
      amrex::Abort("Unrecongized RunOn option");

    return 0;
  }

#if defined(AMREX_USE_GPU)
  template <amrex::RunOn run_on>
#else
  template <amrex::RunOn run_on=amrex::RunOn::Host>
#endif
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real calc_partial_cp_sn (const amrex::Real /*T_s*/,
                                  const int /*n*/) const
  {
    return 0;
  }

#if defined(AMREX_USE_GPU)
  template <amrex::RunOn run_on>
#else
  template <amrex::RunOn run_on=amrex::RunOn::Host>
#endif
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real calc_h_sn (const amrex::Real T_s,
                         const int n) const
  {
    return calc_cp_sn<run_on>(T_s,n) * T_s;
  }

#if defined(AMREX_USE_GPU)
  template <amrex::RunOn run_on>
#else
  template <amrex::RunOn run_on=amrex::RunOn::Host>
#endif
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real calc_partial_h_sn (const amrex::Real T_s,
                                 const int n) const
  {
    return calc_partial_cp_sn<run_on>(T_s,n) * T_s + calc_cp_sn<run_on>(T_s,n);
  }
};


class SolidsPhase
{
  public:
    // Constructor
    SolidsPhase();

    // Destructor
    ~SolidsPhase();

    static constexpr int NMAX = 10;

    int NTYPES;

    // List of recognized fluid specific heat models
    struct SPECIFICHEATMODEL { enum {Invalid=-1, Constant, NASA9Polynomials}; };

    int SpecificHeatModel;

    // List of recognized fluid thermal conductivity models
    struct THERMALCONDUCTIVITYMODEL{ enum {Invalid=-1, Constant}; };

    int ThermalConductivityModel;

    // Names of solids used in IC/BC setups
    amrex::Vector<std::string> names;

    // Specified constant molecular weight of solids
    amrex::Real MW_s0;

    // Specified constant specific heat of solids
    amrex::Gpu::HostVector<amrex::Real> cp_s0;
    amrex::Gpu::DeviceVector<amrex::Real> d_cp_s0;

    // Specified constant enthalpy of formation of solids
    amrex::Real H_f0;

    // Specified constant gas reference temperature
    amrex::Real T_ref;

    // Flag to solve species solids equations
    int solve_species;

    // Particle species
    amrex::Vector<std::string> species;

    // Species unique identifying code
    amrex::Gpu::HostVector<int> species_id;
    amrex::Gpu::DeviceVector<int> d_species_id;

    // Total number of dem species
    int nspecies;

    // Specified constant solid phase species molecular weight
    amrex::Gpu::HostVector<amrex::Real> MW_sn0;
    amrex::Gpu::DeviceVector<amrex::Real> d_MW_sn0;

    // Flag to set solids as a mixture of solids species
    int is_a_mixture;

    // Specified constant species enthalpy of formation
    amrex::Gpu::HostVector<amrex::Real> H_fn0;
    amrex::Gpu::DeviceVector<amrex::Real> d_H_fn0;

    // Specified constant specific heat
    amrex::Gpu::HostVector<amrex::Real> cp_sn0;
    amrex::Gpu::DeviceVector<amrex::Real> d_cp_sn0;

    // Solids parameters for GPU access
    SolidParms* parameters;

    amrex::Real enthalpy_source;

    struct INPUT_DIST_t {
      std::string distribution;
      amrex::Real mean; // Mean value of distriubtion (or constant value)
      amrex::Real std;  // Standard deviation
      amrex::Real min;  // Minimum value of distribution
      amrex::Real max;  // Max value of distribution

      INPUT_DIST_t() :
        mean(-1.0),
        std(-1.0),
        min(-1.0),
        max(-1.0) { }

    };

    struct SOLIDS_t {

      std::string name;
      amrex::Real volfrac;
      amrex::Real statwt;
      amrex::Real temperature;
      amrex::Vector<amrex::Real> velocity;

      amrex::Vector<SPECIES::SPECIES_t> species;

      INPUT_DIST_t diameter;
      INPUT_DIST_t density;

      SOLIDS_t() :
        volfrac(0.0),
        statwt(0.0),
        temperature(0.0),
        velocity(0),
        species(0),
        diameter(),
        density()
      {}
    };

    void Initialize ();
};

#endif
