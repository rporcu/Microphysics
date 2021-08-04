#ifndef MFIX_SOLIDS_PARMS_H_
#define MFIX_SOLIDS_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

#include <mfix_species_parms.H>

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
    const int m_ncoefficients;
    const amrex::Real* m_h_cp_sn;
    const amrex::Real* m_d_cp_sn;
    const amrex::Real* m_h_H_fn;
    const amrex::Real* m_d_H_fn;

  public:
    AMREX_GPU_HOST_DEVICE
    SolidsParms ()
      : m_T_ref(0)
      , m_nspecies(0)
      , m_h_species_id(nullptr)
      , m_d_species_id(nullptr)
      , m_h_MW_sn(nullptr)
      , m_d_MW_sn(nullptr)
      , m_ncoefficients(0)
      , m_h_cp_sn(nullptr)
      , m_d_cp_sn(nullptr)
      , m_h_H_fn(nullptr)
      , m_d_H_fn(nullptr)
    {}

    AMREX_GPU_HOST_DEVICE
    SolidsParms (const amrex::Real  T_ref,
                 const int nspecies,
                 const int* h_species_id,
                 const int* d_species_id,
                 const amrex::Real* h_MW_sn,
                 const amrex::Real* d_MW_sn,
                 const int ncoefficients,
                 const amrex::Real* h_cp_sn,
                 const amrex::Real* d_cp_sn,
                 const amrex::Real* h_H_fn,
                 const amrex::Real* d_H_fn)
      : m_T_ref(T_ref)
      , m_nspecies(nspecies)
      , m_h_species_id(h_species_id)
      , m_d_species_id(d_species_id)
      , m_h_MW_sn(h_MW_sn)
      , m_d_MW_sn(d_MW_sn)
      , m_ncoefficients(ncoefficients)
      , m_h_cp_sn(h_cp_sn)
      , m_d_cp_sn(d_cp_sn)
      , m_h_H_fn(h_H_fn)
      , m_d_H_fn(d_H_fn)
    {}

#if defined(AMREX_USE_GPU)
    template <amrex::RunOn run_on>
#else
    template <amrex::RunOn run_on=amrex::RunOn::Host>
#endif
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_cp_s (const int phase,
                           const amrex::Real T_s) const
    {
      amrex::Real result(0);
      amrex::Real coeff(1);

      if (run_on == amrex::RunOn::Gpu) {
        for (int m(0); m < m_ncoefficients; ++m) {
          const int idx = phase*m_ncoefficients + m;
          result += m_d_cp_sn[idx]*coeff;
          coeff *= T_s;
        }

        //result *= (R/m_d_MW_sn[phase]);

      } else if (run_on == amrex::RunOn::Cpu) {
        for (int m(0); m < m_ncoefficients; ++m) {
          const int idx = phase*m_ncoefficients + m;
          result += m_h_cp_sn[idx]*coeff;
          coeff *= T_s;
        }

        //result *= (R/m_h_MW_sn[phase]);

      } else {
        amrex::Abort("Unrecognized RunOn option");
      }

      return result;
    }

#if defined(AMREX_USE_GPU)
    template <amrex::RunOn run_on>
#else
    template <amrex::RunOn run_on=amrex::RunOn::Host>
#endif
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_integral_cp_s (const int phase,
                                    const amrex::Real T_s) const
    {
      amrex::Real result(0);
      amrex::Real coeff(T_s);

      if (run_on == amrex::RunOn::Gpu) {
        for (int m(0); m < m_ncoefficients; ++m) {
          const int idx = phase*m_ncoefficients + m;
          result += (1./(m+1))*m_d_cp_sn[idx]*coeff;
          coeff *= T_s;
        }

        //result *= (R/m_d_MW_sn[phase]);

      } else if (run_on == amrex::RunOn::Cpu) {
        for (int m(0); m < m_ncoefficients; ++m) {
          const int idx = phase*m_ncoefficients + m;
          result += (1./(m+1))*m_h_cp_sn[idx]*coeff;
          coeff *= T_s;
        }

        //result *= (R/m_h_MW_sn[phase]);

      } else {
        amrex::Abort("Unrecognized RunOn option");
      }

      return result;
    }

#if defined(AMREX_USE_GPU)
    template <amrex::RunOn run_on>
#else
    template <amrex::RunOn run_on=amrex::RunOn::Host>
#endif
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_partial_cp_s (const int phase,
                                   const amrex::Real T_s) const
    {
      amrex::Real result(0);
      amrex::Real coeff(1);

      if (run_on == amrex::RunOn::Gpu) {
        for (int m(1); m < m_ncoefficients; ++m) {
          const int idx = phase*m_ncoefficients + m;
          result += m*m_d_cp_sn[idx]*coeff;
          coeff *= T_s;
        }

        //result *= (R/m_d_MW_sn[phase]);

      } else if (run_on == amrex::RunOn::Cpu) {
        for (int m(1); m < m_ncoefficients; ++m) {
          const int idx = phase*m_ncoefficients + m;
          result += m*m_h_cp_sn[idx]*coeff;
          coeff *= T_s;
        }

        //result *= (R/m_h_MW_sn[phase]);

      } else {
        amrex::Abort("Unrecognized RunOn option");
      }

      return result;
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
      return calc_integral_cp_s<run_on>(phase,T_s) - calc_integral_cp_s<run_on>(phase,m_T_ref);
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
      return calc_cp_s<run_on>(phase,T_s);
    }

#if defined(AMREX_USE_GPU)
    template <amrex::RunOn run_on>
#else
    template <amrex::RunOn run_on=amrex::RunOn::Host>
#endif
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_cp_sn (const amrex::Real T_s,
                            const int n) const
    {
      amrex::Real result(0);
      amrex::Real coeff(1);

      if (run_on == amrex::RunOn::Gpu) {
        for (int m(0); m < m_ncoefficients; ++m) {
          const int idx = n*m_ncoefficients + m;
          result += m_d_cp_sn[idx]*coeff;
          coeff *= T_s;
        }

        //result *= (R/m_d_MW_sn[n]);

      } else if (run_on == amrex::RunOn::Cpu) {
        for (int m(0); m < m_ncoefficients; ++m) {
          const int idx = n*m_ncoefficients + m;
          result += m_h_cp_sn[idx]*coeff;
          coeff *= T_s;
        }

        //result *= (R/m_h_MW_sn[n]);

      } else {
        amrex::Abort("Unrecognized RunOn option");
      }

      return result;
    }

#if defined(AMREX_USE_GPU)
    template <amrex::RunOn run_on>
#else
    template <amrex::RunOn run_on=amrex::RunOn::Host>
#endif
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_integral_cp_sn (const amrex::Real T_s,
                                     const int n) const
    {
      amrex::Real result(0);
      amrex::Real coeff(T_s);

      if (run_on == amrex::RunOn::Gpu) {
        for (int m(0); m < m_ncoefficients; ++m) {
          const int idx = n*m_ncoefficients + m;
          result += (1./(m+1))*m_d_cp_sn[idx]*coeff;
          coeff *= T_s;
        }

        //result *= (R/m_d_MW_sn[n]);

      } else if (run_on == amrex::RunOn::Cpu) {
        for (int m(0); m < m_ncoefficients; ++m) {
          const int idx = n*m_ncoefficients + m;
          result += (1./(m+1))*m_h_cp_sn[idx]*coeff;
          coeff *= T_s;
        }

        //result *= (R/m_h_MW_sn[n]);

      } else {
        amrex::Abort("Unrecognized RunOn option");
      }

      return result;
    }

#if defined(AMREX_USE_GPU)
    template <amrex::RunOn run_on>
#else
    template <amrex::RunOn run_on=amrex::RunOn::Host>
#endif
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_partial_cp_sn (const amrex::Real T_s,
                                    const int n) const
    {
      amrex::Real result(0);
      amrex::Real coeff(1);

      if (run_on == amrex::RunOn::Gpu) {
        for (int m(1); m < m_ncoefficients; ++m) {
          const int idx = n*m_ncoefficients + m;
          result += m*m_d_cp_sn[idx]*coeff;
          coeff *= T_s;
        }

        //result *= (R/m_d_MW_sn[n]);

      } else if (run_on == amrex::RunOn::Cpu) {
        for (int m(1); m < m_ncoefficients; ++m) {
          const int idx = n*m_ncoefficients + m;
          result += m*m_h_cp_sn[idx]*coeff;
          coeff *= T_s;
        }

        //result *= (R/m_h_MW_sn[n]);

      } else {
        amrex::Abort("Unrecognized RunOn option");
      }

      return result;
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
      return calc_integral_cp_sn<run_on>(T_s,n) - calc_integral_cp_sn<run_on>(m_T_ref,n);
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
      return calc_cp_sn<run_on>(T_s,n);
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

    // List of recognized species heat of formation models
    struct ENTHALPYOFFORMATIONMODEL { enum {Invalid=-1, Constant}; };

    int EnthalpyOfFormationModel;

    // Names of solids used in IC/BC setups
    amrex::Vector<std::string> names;

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
    SolidsParms* parameters;

    amrex::Real enthalpy_source;

    int update_mass;

    int update_momentum;

    int update_enthalpy;

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
