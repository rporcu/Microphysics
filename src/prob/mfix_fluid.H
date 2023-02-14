#ifndef _MFIX_FLUID_H_
#define _MFIX_FLUID_H_

#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuControl.H>

#include <mfix_species.H>
#include <mfix_reactions.H>

// Forward declaration
class MFIXFluidParms;


class MFIXFluidPhase
{
  public:
    // 
    struct ViscosityModel { enum {Invalid=-1, Constant, Sutherland}; };

    // 
    struct ThermalConductivityModel { enum {Invalid=-1, Constant}; };

    // Projection constraint type
    struct ConstraintType { enum {Invalid=-1, IncompressibleFluid, IdealGasOpenSystem, IdealGasClosedSystem}; };

    // Ideal gas constant
    static constexpr amrex::Real R = 8.31446261815324;

    // Constructor
    MFIXFluidPhase();

    // Destructor
    ~MFIXFluidPhase();

    void Initialize (const MFIXSpecies& species, const MFIXReactions& reactions);

    const std::string& names (const int i) const { return m_names[i]; }

    bool solve () const { return (m_solve == 1);  }
    bool solve_density () const { return (m_solve_density == 1);  }
    bool solve_tracer () const { return (m_solve_tracer == 1);  }
    bool solve_enthalpy () const { return (m_solve_enthalpy == 1);  }
    bool solve_species () const { return (m_solve_species == 1);  }

    const amrex::Vector<std::string>& species_names () const { return m_species_names; }
    const std::string& species_names (const int n_g) const { return m_species_names[n_g]; }
    int nspecies () const { return m_nspecies; }

    int isInitialized () const { return m_is_initialized; }
    int isMixture () const { return m_is_a_mixture; }

    const MFIXFluidParms& parameters () const { return *m_parameters; }

    int constraint_type () const { return m_constraint_type; }

    amrex::Real mu_g () const { return m_mu_g0; }
    amrex::Real tracer () const { return m_trac_0; }
    amrex::Real cold_flow_temperature () const { return m_T_g0; }
    int p_therm_defined () const { return m_p_therm_defined; }
    amrex::Real thermodynamic_pressure () const { return m_thermodynamic_pressure; }

    void set_p_therm_defined (const int i) { m_p_therm_defined = i; }
    void set_thermodynamic_pressure (const amrex::Real p) { m_thermodynamic_pressure = p; }

  private:
    int m_ntypes;

    // List of recognized fluid viscosity models
    int m_viscosity_model;

    // List of recognized fluid specific heat models
    int m_specific_heat_model;

    // List of recognized fluid thermal conductivity models
    int m_thermal_conductivity_model;

    // Constraint type
    int m_constraint_type;

    // Names of fluids
    amrex::Vector<std::string> m_names;

    // Flag to solve fluid equations
    int m_solve;

    // Flag to solve fluid density equation
    int m_solve_density;

    // Flag to solve fluid tracer equation
    int m_solve_tracer;

    // Specified constant tracer value
    amrex::Real m_trac_0;

    // Specified constant gas viscosity
    amrex::Real m_mu_g0;

    // Average molecular weight of gas
    amrex::Real m_mw_avg;

    // Flag to solve fluid enthalpy equation
    int m_solve_enthalpy;

    // Specified constant gas temperature
    amrex::Real m_T_g0;

    // Specified constant gas reference temperature
    amrex::Real m_T_ref;

    // Specified constant gas phase thermal conductivity coefficient
    amrex::Real m_k_g0;

    // Specified gas phase thermodynamic pressure
    amrex::Real m_thermodynamic_pressure;
    int m_p_therm_defined;

    // Flag to solve species fluid equations
    int m_solve_species;

    // Fluid phase species names
    amrex::Vector<std::string> m_species_names;

    // Species unique identifying code
    amrex::Gpu::HostVector<int> m_species_IDs;
    amrex::Gpu::DeviceVector<int> m_d_species_IDs;

    // Total number of fluid species
    int m_nspecies;

    // Specified constant gas phase species molecular weight
    amrex::Gpu::HostVector<amrex::Real> m_MW_gk0;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_MW_gk0;

    // Specified constant gas phase species diffusion coefficients
    amrex::Real m_D_g0;

    // Flag to understand if fluid is a mixture of fluid species
    int m_is_a_mixture;

    // Specified constant species enthalpy of formation
    amrex::Gpu::HostVector<amrex::Real> m_H_fk0;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_H_fk0;

    // Specified constant gas phase species specific heat
    amrex::Gpu::HostVector<amrex::Real> m_cp_gk0;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_cp_gk0;

    // Tables for stoichiometric coefficients for each species in each reaction
    amrex::Gpu::HostVector<amrex::Real> m_stoich_coeffs;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_stoich_coeffs;

    // Fluid parameters for GPU access
    MFIXFluidParms* m_parameters;

    int m_is_initialized;

}; // class MFIXFluidPhase


class MFIXFluidParms
{
  public:
    AMREX_GPU_HOST_DEVICE
    MFIXFluidParms (const amrex::Real T_ref);

    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_MW_gk (const int n) const
    {
      if (run_on == amrex::RunOn::Gpu)
        return m_d_MW_gk[n];
      else if (run_on == amrex::RunOn::Cpu)
        return m_h_MW_gk[n];
      else
        amrex::Abort("Unrecognized RunOn option");

      return 0;
    }

    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_MW_g () const
    {
      return get_MW_gk<run_on>(0);
    }

    // Compute fluid molecular viscosity
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_mu_g (const amrex::Real /*T_g*/) const
    {
      return m_mu_g;
    }

    // Compute fluid thermal conductivity
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_k_g (const amrex::Real /*T_g*/) const
    {
      return m_k_g;
    }

    // Compute fluid thermal conductivity
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_partial_k_g (const amrex::Real /*T_g*/) const
    {
      return 0;
    }

    // Compute fluid species diffusivity
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_D_g () const
    {
      return m_D_g;
    }

    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_cp_g (const amrex::Real T_g) const
    {
      return calc_cp_gk<run_on>(T_g, 0);
    }

    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_h_g (const amrex::Real T_g,
                          const int cell_is_covered = 0) const
    {
      if (cell_is_covered)
        return 0;

      return calc_h_gk<run_on>(T_g, 0);
    }

    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_partial_h_g (const amrex::Real T_g) const
    {
      return calc_partial_h_gk<run_on>(T_g, 0);
    }

    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_cp_gk (const amrex::Real T_g,
                            const int n) const
    {
      amrex::Real result(0);

      if (m_specific_heat_model == MFIXSpecies::SpecificHeatModel::Constant) {

        if (run_on == amrex::RunOn::Gpu) {

          result = m_d_cp_gk[n];

        } else if (run_on == amrex::RunOn::Cpu) {

          result = m_h_cp_gk[n];
        }

      } else if (m_specific_heat_model == MFIXSpecies::SpecificHeatModel::NASA7Polynomials) {

        amrex::Real coeff(1);

        for (int m(0); m < 5; ++m) {

          int idx = T_g < 1000 ? n*12 + m : n*12 + m+6;

          if (run_on == amrex::RunOn::Gpu) {

            result += m_d_cp_gk[idx]*coeff;

          } else if (run_on == amrex::RunOn::Cpu) {

            result += m_h_cp_gk[idx]*coeff;

          }

          coeff *= T_g;

        }
      }

      return result;
    }

    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_h_gk (const amrex::Real T_g,
                           const int n,
                           const int cell_is_covered = 0) const
    {
      if (cell_is_covered)
        return 0;

      amrex::Real result(0);

      if (m_specific_heat_model == MFIXSpecies::SpecificHeatModel::Constant) {

        if (run_on == amrex::RunOn::Gpu) {

          result = m_d_H_fk[n] + m_d_cp_gk[n]*(T_g - m_T_ref);

        } else if (run_on == amrex::RunOn::Cpu) {

          result = m_h_H_fk[n] + m_h_cp_gk[n]*(T_g - m_T_ref);
        }

      } else if (m_specific_heat_model == MFIXSpecies::SpecificHeatModel::NASA7Polynomials) {

        amrex::Real coeff(T_g);

        for (int m(0); m < 5; ++m) {

          int idx = T_g < 1000 ? n*12 + m : n*12 + m+6;

          if (run_on == amrex::RunOn::Gpu) {

            result += (1./(1.+m))*m_d_cp_gk[idx]*coeff;

          } else if (run_on == amrex::RunOn::Cpu) {

            result += (1./(1.+m))*m_h_cp_gk[idx]*coeff;

          }

          coeff *= T_g;
        }

        int idx = T_g < 1000 ? n*12 + 5 : n*12 + 11;

        if (run_on == amrex::RunOn::Gpu) {

          result += m_d_cp_gk[idx];

        } else if (run_on == amrex::RunOn::Cpu) {

          result += m_h_cp_gk[idx];

        }

      }

      return result;
    }

    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_partial_h_gk (const amrex::Real T_g,
                                   const int n) const
    {
      return calc_cp_gk<run_on>(T_g, n);
    }

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

  private:
    // Friendship so MFIXFluidPhase can access private members
    friend class MFIXFluidPhase;

    amrex::Real  m_T_ref;
    amrex::Real  m_mu_g;
    amrex::Real  m_k_g;
    int          m_nspecies;
    int*         m_h_species_id;
    int*         m_d_species_id;

    amrex::Real* m_h_MW_gk;
    amrex::Real* m_d_MW_gk;
    amrex::Real m_D_g;
    amrex::Real* m_h_cp_gk;
    amrex::Real* m_d_cp_gk;
    amrex::Real* m_h_H_fk;
    amrex::Real* m_d_H_fk;
    int          m_nreactions;
    amrex::Real* m_h_stoich_coeffs;
    amrex::Real* m_d_stoich_coeffs;

    // List of recognized fluid specific heat models
    int m_specific_heat_model;

};


struct FLUID_t {

  FLUID_t(MFIXFluidPhase* fluid_in) :
    fluid(fluid_in),
    volfrac(-1.0),
    pressure(-1.0),
    pressure_defined(0),
    constant_density(true),
    density(-1.0),
    density_defined(0),
    species(0),
    constant_velocity(true),
    flow_thru_eb(false),
    eb_has_velocity(false),
    eb_vel_is_mag(false),
    eb_has_volflow(false),
    constant_volflow(true),
    volflow(0.0),
    constant_temperature(true),
    temperature(0.0),
    temperature_defined(0)
  {}

  MFIXFluidPhase* fluid;

  amrex::Real volfrac;
  amrex::Real pressure;
  int pressure_defined; // Flag to check if pressure is defined

  bool constant_density;
  amrex::Real density;
  int density_defined; // Flag to check if density is defined
  amrex::Vector<amrex::Vector<amrex::Real>> density_table;

  amrex::Real get_density (amrex::Real time = 0) const;

  amrex::Vector< int > constant_species;
  amrex::Vector<SPECIES_t> species;
  amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> species_table;

  amrex::Real get_species (int n, amrex::Real time = 0) const;

  // Flags for velocity
  bool constant_velocity;
  amrex::Vector<amrex::Real> velocity;
  amrex::Vector<amrex::Vector<amrex::Real>> vel_table;

  amrex::RealVect get_velocity (amrex::Real time = 0) const;

  bool flow_thru_eb;
  bool eb_has_velocity;
  bool eb_vel_is_mag;

  amrex::Real get_velocity_mag (amrex::Real time = 0) const;

  bool eb_has_volflow;
  bool constant_volflow;
  amrex::Real volflow;
  amrex::Vector<amrex::Vector<amrex::Real>> volflow_table;

  amrex::Real get_volflow (amrex::Real time = 0) const;

  bool constant_temperature;
  amrex::Real temperature;
  int temperature_defined;
  amrex::Vector<amrex::Vector<amrex::Real>> tg_table;

  amrex::Real get_temperature (amrex::Real time = 0) const;
};


#endif
