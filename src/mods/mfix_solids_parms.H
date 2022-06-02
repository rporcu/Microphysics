#ifndef MFIX_SOLIDS_PARMS_H_
#define MFIX_SOLIDS_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_RandomEngine.H>

#include <mfix_species_parms.H>
#include <mfix_reactions_parms.H>


class SolidsParms {

  public:
    static constexpr amrex::Real R =  8.31446261815324;

    const amrex::Real  m_T_ref;
    const int          m_nspecies;
    const int*         m_h_species_id;
    const int*         m_d_species_id;

  private:
    const amrex::Real* m_h_MW_sn;
    const amrex::Real* m_d_MW_sn;
//    const int m_ncoefficients;
    const amrex::Real* m_h_cp_sn;
    const amrex::Real* m_d_cp_sn;
    const amrex::Real* m_h_H_fn;
    const amrex::Real* m_d_H_fn;
    const int          m_nreactions;
    const amrex::Real* m_h_stoich_coeffs;
    const amrex::Real* m_d_stoich_coeffs;
    const amrex::Real* m_h_kp_sn;
    const amrex::Real* m_d_kp_sn;
    const amrex::Real  m_flpc;
    const amrex::Real  m_rough;
    const int          m_do_pfp_cond;


    // List of recognized particle specific heat models
    struct SPECIFICHEATMODEL { enum {Invalid=-1, Constant, NASA7Polynomials, NASA9Polynomials};  };
    int m_specific_heat_model;

    // List of recognized particle conductivity models
    struct CONDUCTIVITYMODEL { enum {Invalid=-1, Constant};  };
    int m_conductivity_model;

  public:
    AMREX_GPU_HOST_DEVICE
    SolidsParms ()
      : m_T_ref(0)
      , m_nspecies(0)
      , m_h_species_id(nullptr)
      , m_d_species_id(nullptr)
      , m_h_MW_sn(nullptr)
      , m_d_MW_sn(nullptr)
//      , m_ncoefficients(0)
      , m_h_cp_sn(nullptr)
      , m_d_cp_sn(nullptr)
      , m_h_H_fn(nullptr)
      , m_d_H_fn(nullptr)
      , m_nreactions(0)
      , m_h_stoich_coeffs(nullptr)
      , m_d_stoich_coeffs(nullptr)
      , m_h_kp_sn(nullptr)
      , m_d_kp_sn(nullptr)
      , m_flpc(0.4)
      , m_rough(2E-8)
      , m_do_pfp_cond(0)
      , m_specific_heat_model(SPECIFICHEATMODEL::Invalid)
      , m_conductivity_model(CONDUCTIVITYMODEL::Invalid)
    {}

    AMREX_GPU_HOST_DEVICE
    SolidsParms (const amrex::Real T_ref,
                 const int nspecies,
                 const int* h_species_id,
                 const int* d_species_id,
                 const amrex::Real* h_MW_sn,
                 const amrex::Real* d_MW_sn,
                 const int /*ncoefficients*/,
                 const amrex::Real* h_cp_sn,
                 const amrex::Real* d_cp_sn,
                 const amrex::Real* h_H_fn,
                 const amrex::Real* d_H_fn,
                 const int nreactions,
                 const amrex::Real* h_stoich_coeffs,
                 const amrex::Real* d_stoich_coeffs,
                 const amrex::Real* h_kp_sn,
                 const amrex::Real* d_kp_sn,
                 const amrex::Real flpc,
                 const amrex::Real rough,
                 const int do_pfp_cond,
                 const int specific_heat_model,
                 const int conductivity_model)
      : m_T_ref(T_ref)
      , m_nspecies(nspecies)
      , m_h_species_id(h_species_id)
      , m_d_species_id(d_species_id)
      , m_h_MW_sn(h_MW_sn)
      , m_d_MW_sn(d_MW_sn)
//      , m_ncoefficients(ncoefficients)
      , m_h_cp_sn(h_cp_sn)
      , m_d_cp_sn(d_cp_sn)
      , m_h_H_fn(h_H_fn)
      , m_d_H_fn(d_H_fn)
      , m_nreactions(nreactions)
      , m_h_stoich_coeffs(h_stoich_coeffs)
      , m_d_stoich_coeffs(d_stoich_coeffs)
      , m_h_kp_sn(h_kp_sn)
      , m_d_kp_sn(d_kp_sn)
      , m_flpc(flpc)
      , m_rough(rough)
      , m_do_pfp_cond(do_pfp_cond)
      , m_specific_heat_model(specific_heat_model)
      , m_conductivity_model(conductivity_model)
    {}

    // Driver molecular weight
    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_MW_s () const
    {
      return get_MW_sn<run_on>(0);
    }
    
    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_MW_sn (const int n) const
    {
      if (run_on == amrex::RunOn::Gpu)
        return m_d_MW_sn[n];
      else if (run_on == amrex::RunOn::Cpu)
        return m_h_MW_sn[n];
      else
        amrex::Abort("Unrecognized RunOn option");

      return 0;
    }

    // Get reaction coefficients
    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_stoich_coeff (const int n, const int q) const
    {
      if (run_on == amrex::RunOn::Gpu)
        return m_d_stoich_coeffs[n*m_nreactions+q];
      else if (run_on == amrex::RunOn::Cpu)
        return m_h_stoich_coeffs[n*m_nreactions+q];
      else
        amrex::Abort("Unrecognized RunOn option");

      return 0;
    }
    
    // Driver heat capacity
    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_cp_s (const int phase,
                           const amrex::Real T_s) const
    {
      return calc_cp_sn<run_on>(T_s, phase);
    }

    // Driver enthalpy
    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_h_s (const int phase,
                          const amrex::Real T_s) const
    {
      return calc_h_sn<run_on>(T_s, phase);
    }

    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_partial_h_s (const int phase,
                                  const amrex::Real T_s) const
    {
      return calc_partial_h_sn<run_on>(T_s, phase);
    }

    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_partial_h_sn (const amrex::Real T_s,
                                   const int n) const
    {
      return calc_cp_sn<run_on>(T_s, n);
    }

    // Get FLPC
    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_flpc () const
    {
      return m_flpc;
    }

    // Get min conduction distance
    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_min_cond () const
    {
      return m_rough;
    }

    // Do we enter PFP conduction?
    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int get_do_pfp_cond () const
    {
      return m_do_pfp_cond;
    }

    //************************************************************

    // Calc heat capacity
    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_cp_sn (const amrex::Real T_s,
                            const int n) const
    {
      amrex::Real result(0);

      if (m_specific_heat_model == SPECIFICHEATMODEL::Constant) {

        if (run_on == amrex::RunOn::Gpu) {

          result = m_d_cp_sn[n];

        } else if (run_on == amrex::RunOn::Cpu) {

          result = m_h_cp_sn[n];
        }

      } else if (m_specific_heat_model == SPECIFICHEATMODEL::NASA7Polynomials) {

        amrex::Real coeff(1);

        for (int m(0); m < 5; ++m) {

          int idx = T_s < 1000 ? n*12 + m : n*12 + m+6;

          if (run_on == amrex::RunOn::Gpu) {

            result += m_d_cp_sn[idx]*coeff;

          } else if (run_on == amrex::RunOn::Cpu) {

            result += m_h_cp_sn[idx]*coeff;
          }

          coeff *= T_s;
        }

      }

      return result;
    }

    // Calc enthalpy
    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_h_sn (const amrex::Real T_s,
                           const int n) const
    {
      amrex::Real result(0);

      if (m_specific_heat_model == SPECIFICHEATMODEL::Constant) {

        if (run_on == amrex::RunOn::Gpu) {

          result = m_d_H_fn[n] + m_d_cp_sn[n]*(T_s - m_T_ref);

        } else if (run_on == amrex::RunOn::Cpu) {

          result = m_h_H_fn[n] + m_h_cp_sn[n]*(T_s - m_T_ref);
        }

      } else if (m_specific_heat_model == SPECIFICHEATMODEL::NASA7Polynomials) {

        amrex::Real coeff(T_s);

        for (int m(0); m < 5; ++m) {

          int idx = T_s < 1000 ? n*12 + m : n*12 + m+6;

          if (run_on == amrex::RunOn::Gpu) {

            result += (1./(m+1.))*m_d_cp_sn[idx]*coeff;

          } else if (run_on == amrex::RunOn::Cpu) {

            result += (1./(m+1.))*m_h_cp_sn[idx]*coeff;

          }

          coeff *= T_s;
        }

        int idx = T_s < 1000 ? n*12 + 5 : n*12 + 11;

        if (run_on == amrex::RunOn::Gpu) {

          result += m_d_cp_sn[idx];

        } else if (run_on == amrex::RunOn::Cpu) {

          result += m_h_cp_sn[idx];

        }

      }

      return result;
    }

    // Calc conductivity
    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_kp_sn (const amrex::Real,
                            const int n) const
    {
      amrex::Real result(0);

      if (m_conductivity_model == CONDUCTIVITYMODEL::Constant) {

        if (run_on == amrex::RunOn::Gpu) {
          
          result = m_d_kp_sn[n];

        } else if (run_on == amrex::RunOn::Cpu) {

          result = m_h_kp_sn[n];
        }

      }

      return result;
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
    struct SPECIFICHEATMODEL { enum {Invalid=-1, Constant, NASA7Polynomials, NASA9Polynomials}; };
    int SpecificHeatModel;

    // List of recognized fluid thermal conductivity models
    struct THERMALCONDUCTIVITYMODEL{ enum {Invalid=-1, Constant}; };
    int ThermalConductivityModel;

    // Names of solids used in IC/BC setups
    amrex::Vector<std::string> names;

    // Flag to turn on/off particles mass update
    int solve_mass;

    // Flag to turn on/off particles velocity update
    int solve_momentum;

    // Specified constant gas reference temperature
    amrex::Real T_ref;

    // User-defined enthalpy volumetric source for particles
    amrex::Real enthalpy_source;

    // Flag to turn on/off particles enthalpy update
    int solve_enthalpy;

    // Flag to solve species solids equations
    int solve_species;

    // Particle species
    amrex::Vector<std::string> species_names;

    // Species unique identifying code
    amrex::Gpu::HostVector<int> species_IDs;
    amrex::Gpu::DeviceVector<int> d_species_IDs;

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

    // Tables for stoichiometric coefficients for each species in each reaction
    amrex::Gpu::HostVector<amrex::Real> stoich_coeffs;
    amrex::Gpu::DeviceVector<amrex::Real> d_stoich_coeffs;

    // Specified conductivity
    amrex::Gpu::HostVector<amrex::Real> kp_sn0;
    amrex::Gpu::DeviceVector<amrex::Real> d_kp_sn0;

    // Specified fluid-lens constant
    amrex::Real flpc;

    // Specified roughness constant
    amrex::Real rough;

    // Flag to do pfp conduction
    int do_pfp_cond;

    // Solids parameters for GPU access
    SolidsParms* parameters;

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
        max(-1.0)
      {}

      bool is_constant () { return (distribution.compare("constant") == 0);} ;
      bool is_uniform  () { return (distribution.compare("uniform")  == 0);} ;
      bool is_normal   () { return (distribution.compare("normal")   == 0);} ;

    };

    int is_initialized;

    void Initialize (const Species& species,
                     const Reactions& reactions);
};

struct SOLIDS_t {

  std::string name;
  amrex::Real volfrac;
  amrex::Real statwt;
  amrex::Real temperature;
  amrex::Vector<amrex::Real> velocity;
  amrex::Real velmag;
  amrex::Real volflow;
  amrex::Real vol_remainder;
  //int inflow_count;

  amrex::Vector<SPECIES_t> species;

  SolidsPhase::INPUT_DIST_t diameter;
  SolidsPhase::INPUT_DIST_t density;

  SOLIDS_t() :
    volfrac(0.0),
    statwt(0.0),
    temperature(0.0),
    velocity(0),
    velmag(0),
    volflow(0.),
    vol_remainder(0.),
    //inflow_count(0),
    species(0),
    diameter(),
    density()
  {}
};

#endif
