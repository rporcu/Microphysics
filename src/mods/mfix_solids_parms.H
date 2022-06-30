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


// Forward declaration
class SolidsPhase;

struct SOLIDSTHERMALCONDUCTIVITYMODEL { enum {Invalid=-1, Constant}; };


class SolidsParms {

  private:
    friend class SolidsPhase;

    static constexpr amrex::Real R =  8.31446261815324;

    amrex::Real  m_T_ref;

    int* m_h_species_id;
    int* m_d_species_id;

    amrex::Real* m_h_MW_sn;
    amrex::Real* m_d_MW_sn;

    amrex::Real* m_h_cp_sn;
    amrex::Real* m_d_cp_sn;

    amrex::Real* m_h_H_fn;
    amrex::Real* m_d_H_fn;

    int  m_nreactions;

    amrex::Real* m_h_stoich_coeffs;
    amrex::Real* m_d_stoich_coeffs;

    amrex::Real* m_h_kp_sn;
    amrex::Real* m_d_kp_sn;

    amrex::Real m_flpc;

    amrex::Real m_rough;

    int m_do_pfp_cond;

    int m_specific_heat_model;
    int m_thermal_conductivity_model;

  public:
    AMREX_GPU_HOST_DEVICE
    SolidsParms (const amrex::Real T_ref)
      : m_T_ref(T_ref)
      , m_h_species_id(nullptr)
      , m_d_species_id(nullptr)
      , m_h_MW_sn(nullptr)
      , m_d_MW_sn(nullptr)
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
      , m_thermal_conductivity_model(SOLIDSTHERMALCONDUCTIVITYMODEL::Invalid)
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
    amrex::Real calc_cp_s (const amrex::Real T_s) const
    {
      return calc_cp_sn<run_on>(T_s, 0);
    }

    // Driver enthalpy
    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_h_s (const amrex::Real T_s) const
    {
      return calc_h_sn<run_on>(T_s, 0);
    }

    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_partial_h_s (const amrex::Real T_s) const
    {
      return calc_partial_h_sn<run_on>(T_s, 0);
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

      if (m_thermal_conductivity_model == SOLIDSTHERMALCONDUCTIVITYMODEL::Constant) {

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

    int ntypes;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static int phase_to_index (const int phase)
    { return phase-1; }

    // List of recognized fluid specific heat models
    int specific_heat_model;

    // List of recognized fluid thermal conductivity models
    int thermal_conductivity_model;

    // Names of solids used in IC/BC setups
    amrex::Vector<std::string> names;

    // Phases IDs for each solids type
    amrex::Vector<int> phases;

    // Return phase ID given a solids type
    int name_to_phase (const std::string& name) const;

    // Flag to turn on/off particles mass update
    int update_mass;

    // Flag to solve species
    int solve_species;

    // Flag to turn on/off particles velocity update
    int update_momentum;

    // Specified constant gas reference temperature
    amrex::Real T_ref;

    // User-defined enthalpy volumetric source for particles
    amrex::Real enthalpy_source;

    // Flag to turn on/off particles enthalpy update
    int update_enthalpy;

    // Flag to solve enthalpy
    int solve_enthalpy;

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

    // Specified constant specific heat
    amrex::Gpu::HostVector<amrex::Real> cp_sn0;
    amrex::Gpu::DeviceVector<amrex::Real> d_cp_sn0;

    // Specified constant species enthalpy of formation
    amrex::Gpu::HostVector<amrex::Real> H_fn0;
    amrex::Gpu::DeviceVector<amrex::Real> d_H_fn0;

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

    int is_initialized;

    void Initialize (const Species& species,
                     const Reactions& reactions);
};


struct INPUT_DIST_t {

  INPUT_DIST_t()
    : m_is_constant(0)
    , m_is_uniform(0)
    , m_is_normal(0)
    , m_mean(-1.0)
    , m_stddev(-1.0)
    , m_min(-1.0)
    , m_max(-1.0)
  {}

  int set (const std::string a_field, const std::string a_pprop);

  int is_constant () const { return m_is_constant;};
  int is_uniform  () const { return m_is_uniform; };
  int is_normal   () const { return m_is_normal;  };

  amrex::Real get_mean ()   const { return m_mean;   };
  amrex::Real get_stddev () const { return m_stddev; };
  amrex::Real get_max ()    const { return m_max;    };
  amrex::Real get_min ()    const { return m_min;    };

  private:
    int m_is_constant;
    int m_is_uniform;
    int m_is_normal;

    amrex::Real m_mean;   // Mean value of distriubtion (or constant value)
    amrex::Real m_stddev; // Standard deviation
    amrex::Real m_min;    // Minimum value of distribution
    amrex::Real m_max;    // Max value of distribution
};

struct SOLIDS_t {

  std::string name;
  int phase;
  amrex::Real volfrac;
  amrex::Real statwt;
  amrex::Real temperature;
  amrex::Vector<amrex::Real> velocity;
  amrex::Real velmag;
  amrex::Real volflow;
  amrex::Real vol_remainder;

  amrex::Vector<SPECIES_t> species;

  INPUT_DIST_t diameter;
  INPUT_DIST_t density;

  SOLIDS_t() :
    phase(0),
    volfrac(0.0),
    statwt(0.0),
    temperature(0.0),
    velocity(0),
    velmag(0),
    volflow(0.),
    vol_remainder(0.),
    species(0),
    diameter(),
    density()
  {}
};

#endif
