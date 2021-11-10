#ifndef MFIX_FLUID_PARMS_H_
#define MFIX_FLUID_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuControl.H>

#include <mfix_species_parms.H>

class FluidParms {

  public:
    static constexpr amrex::Real R =  8.31446261815324;

    const amrex::Real  m_T_ref;
    const amrex::Real  m_mu_g;
    const amrex::Real  m_k_g;
    const int          m_nspecies;
    const int*         m_h_species_id;
    const int*         m_d_species_id;

  private:
    const amrex::Real* m_h_MW_gk;
    const amrex::Real* m_d_MW_gk;
    const amrex::Real m_D_g;
//    const int m_ncoefficients;
    const amrex::Real* m_h_cp_gk;
    const amrex::Real* m_d_cp_gk;
    const amrex::Real* m_h_H_fk;
    const amrex::Real* m_d_H_fk;

    // List of recognized fluid specific heat models
    struct SPECIFICHEATMODEL { enum {Invalid=-1, Constant, NASA7Polynomials, NASA9Polynomials}; };
    int m_specific_heat_model;

  public:
    AMREX_GPU_HOST_DEVICE
    FluidParms ()
      : m_T_ref(0.)
      , m_mu_g(0.)
      , m_k_g(0.)
      , m_nspecies(0)
      , m_h_species_id(nullptr)
      , m_d_species_id(nullptr)
      , m_h_MW_gk(nullptr)
      , m_d_MW_gk(nullptr)
      , m_D_g(0)
//      , m_ncoefficients(0)
      , m_h_cp_gk(nullptr)
      , m_d_cp_gk(nullptr)
      , m_h_H_fk(nullptr)
      , m_d_H_fk(nullptr)
      , m_specific_heat_model(SPECIFICHEATMODEL::Invalid)
    {}

    AMREX_GPU_HOST_DEVICE
    FluidParms (const amrex::Real T_ref,
                const amrex::Real mu_g,
                const amrex::Real k_g,
                const int nspecies,
                const int* h_species_id,
                const int* d_species_id,
                const amrex::Real* h_MW_gk,
                const amrex::Real* d_MW_gk,
                const amrex::Real D_g,
                const int /*ncoefficients*/,
                const amrex::Real* h_cp_gk,
                const amrex::Real* d_cp_gk,
                const amrex::Real* h_H_fk,
                const amrex::Real* d_H_fk,
                const int specific_heat_model)
      : m_T_ref(T_ref)
      , m_mu_g(mu_g)
      , m_k_g(k_g)
      , m_nspecies(nspecies)
      , m_h_species_id(h_species_id)
      , m_d_species_id(d_species_id)
      , m_h_MW_gk(h_MW_gk)
      , m_d_MW_gk(d_MW_gk)
      , m_D_g(D_g)
//      , m_ncoefficients(ncoefficients)
      , m_h_cp_gk(h_cp_gk)
      , m_d_cp_gk(d_cp_gk)
      , m_h_H_fk(h_H_fk)
      , m_d_H_fk(d_H_fk)
      , m_specific_heat_model(specific_heat_model)
    {}


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

      if (m_specific_heat_model == SPECIFICHEATMODEL::Constant) {

        if (run_on == amrex::RunOn::Gpu) {
          
          result = m_d_cp_gk[n];

        } else if (run_on == amrex::RunOn::Cpu) {

          result = m_h_cp_gk[n];
        }

      } else if (m_specific_heat_model == SPECIFICHEATMODEL::NASA7Polynomials) {

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

      if (m_specific_heat_model == SPECIFICHEATMODEL::Constant) {

        if (run_on == amrex::RunOn::Gpu) {
          
          result = m_d_H_fk[n] + m_d_cp_gk[n]*(T_g - m_T_ref);

        } else if (run_on == amrex::RunOn::Cpu) {

          result = m_h_H_fk[n] + m_h_cp_gk[n]*(T_g - m_T_ref);
        }

      } else if (m_specific_heat_model == SPECIFICHEATMODEL::NASA7Polynomials) {

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
};


class FluidPhase {

  public:
    // Constructor
    FluidPhase();

    // Destructor
    ~FluidPhase();

    // List of recognized fluid density models
    struct DENSITYMODEL { enum {Invalid=-1, Constant, IdealGas}; };

    int DensityModel;

    // List of recognized fluid viscosity models
    struct VISCOSITYMODEL { enum {Invalid=-1, Constant, Sutherland}; };

    int ViscosityModel;

    // List of recognized species diffusivity models
    struct DIFFUSIVITYMODEL { enum {Invalid=-1, Constant}; };

    int DiffusivityModel;

    // List of recognized fluid specific heat models
    struct SPECIFICHEATMODEL { enum {Invalid=-1, Constant, NASA7Polynomials, NASA9Polynomials}; };

    int SpecificHeatModel;

    // List of recognized fluid thermal conductivity models
    struct THERMALCONDUCTIVITYMODEL{ enum {Invalid=-1, Constant}; };

    int ThermalConductivityModel;

    // Ideal gas constant
    static constexpr amrex::Real R = 8.31446261815324;

    // Flag to solve fluid equations
    int solve;

    // Flag to solve fluid density equation
    int solve_density;

    // Specified constant gas density
    amrex::Real ro_g0;

    // Flag to solve fluid tracer equation
    int solve_tracer;

    // Specified constant tracer value
    amrex::Real trac_0;

    // Specified constant gas viscosity
    amrex::Real mu_g0;

    // Average molecular weight of gas
    amrex::Real mw_avg;

    // Flag to solve fluid enthalpy equation
    int solve_enthalpy;

    // Specified constant gas temperature
    amrex::Real T_g0;

    // Specified constant gas reference temperature
    amrex::Real T_ref;

    // Specified constant gas phase thermal conductivity coefficient
    amrex::Real k_g0;

    // Flag to solve species fluid equations
    int solve_species;

    // Fluid phase species names
    amrex::Vector<std::string> species;

    // Species unique identifying code
    amrex::Gpu::HostVector<int> species_id;
    amrex::Gpu::DeviceVector<int> d_species_id;

    // Total number of fluid species
    int nspecies;

    // Specified constant gas phase species molecular weight
    amrex::Gpu::HostVector<amrex::Real> MW_gk0;
    amrex::Gpu::DeviceVector<amrex::Real> d_MW_gk0;

    // Specified constant gas phase species diffusion coefficients
    amrex::Real D_g0;

    // Flag to understand if fluid is a mixture of fluid species
    int is_a_mixture;

    // Specified constant species enthalpy of formation
    amrex::Gpu::HostVector<amrex::Real> H_fk0;
    amrex::Gpu::DeviceVector<amrex::Real> d_H_fk0;

    // Specified constant gas phase species specific heat
    amrex::Gpu::HostVector<amrex::Real> cp_gk0;
    amrex::Gpu::DeviceVector<amrex::Real> d_cp_gk0;

    // Name to later reference when building input regions
    std::string name;

    // Fluid parameters for GPU access
    FluidParms* parameters;

    struct FLUID_t {

      amrex::Real volfrac;
      amrex::Real pressure;
      amrex::Real pressure_defined; // Flag to check if pressure is defined


      amrex::Vector< int > constant_species;
      amrex::Vector<SPECIES::SPECIES_t> species;
      amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> species_table;

      amrex::Real get_species (int n, amrex::Real time=0.0) {
        //amrex::Print() << "\n\nLooking for the right BC fluid temperature!\n\n";
        if (constant_species[n]){
          // amrex::Print() << "* Use the constant temperature " << temperature << "\n";
          return species[n].mass_fraction;
        } else {
          const int size = species_table[n].size();
          if( size == 1 || time <= species_table[n][0][0] ) {
            // There is only one entry (use it) OR the current time
            // is less than the first entry so return the first entry.
            // amrex::Print() << "Defaulting to the first entry " << species_table[0][1] << "\n";
            return species_table[n][0][1];

          } else {

            constexpr amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();
            //amrex::Print() << "Interpolating from species table (maybe?)\n";
            for (int lc=1; lc<size; lc++) {
              // Time of current entry
              const amrex::Real entry_time = species_table[n][lc][0];

              if ( amrex::Math::abs(time - entry_time) <= tolerance ) {
                // Too close for any calculations
                // amrex::Print() << "Too close to tell using value form entry "
                //                << lc << "  " << species_table[n][lc][1] << "\n";
                return species_table[n][lc][1];

                // The entries are sorted based on time, so if the current time is
                // less than this entries time, then we need to interpolate between
                // the previous time and this time.
              } else if (time < entry_time) {
                const amrex::Real prv_time = species_table[n][lc-1][0];
                const amrex::Real delta_time = entry_time - prv_time;
                if (amrex::Math::abs(delta_time) <= tolerance) {
                  // amrex::Print() << "Delta time too small. Using value form entry "
                  //                << lc << "  " << species_table[n][lc][1] << "\n";
                  return species_table[n][lc][1];
                } else {
                  const amrex::Real whi = (time - prv_time) / delta_time;
                  const amrex::Real wlo = 1.0 - whi;
                  // amrex::Print() << "Interpolating from table entry " << lc << "\n";
                  return (whi*species_table[n][lc][1] + wlo*species_table[n][lc-1][1]);
                }
              }
            }
            // If we made it here, time is larger than all the entries in the
            // table, therefore, use the last value.
            // amrex::Print() << "Using last table entry " << species_table[n][size-1][1];
            return species_table[n][size-1][1];
          }
        }
      }

      // Flags for constant
      bool constant_velocity;
      amrex::Vector<amrex::Real> velocity;
      amrex::Vector<amrex::Vector<amrex::Real>> vel_table;

      amrex::RealVect get_velocity (amrex::Real time=0.0) {
        //amrex::Print() << "\n\nLooking for the right BC fluid velocity!\n\n";
        if (constant_velocity){
          //amrex::Print() << "* Use the constant velocity!";
          return amrex::RealVect(velocity[0], velocity[1], velocity[2]);
        } else {
          const int size = vel_table.size();
          if( size == 1 || time <= vel_table[0][0] ) {
            // There is only one entry (use it) OR the current time
            // is less than the first entry so return the first entry.
            //amrex::Print() << "Defaulting to the first entry!\n";
            return amrex::RealVect(vel_table[0][1], vel_table[0][2], vel_table[0][3]);

          } else {

            constexpr amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();
            //amrex::Print() << "Interpolating from the velocity table (maybe?)\n";
            for (int lc=1; lc<size; lc++) {
              // Time of current entry
              const amrex::Real entry_time = vel_table[lc][0];

              if ( amrex::Math::abs(time - entry_time) <= tolerance ) {
                // Too close for any calculations
                // amrex::Print() << "Too close to tell using value form entry " << lc << "\n";
                return amrex::RealVect(vel_table[lc][1], vel_table[lc][2], vel_table[lc][3]);

                // The entries are sorted based on time, so if the current time is
                // less than this entries time, then we need to interpolate between
                // the previous time and this time.
              } else if (time < entry_time) {
                const amrex::Real prv_time = vel_table[lc-1][0];
                const amrex::Real delta_time = entry_time - prv_time;
                if (amrex::Math::abs(delta_time) <= tolerance) {
                  //amrex::Print() << "Delta time too small. Using value form entry " << lc << "\n";
                  return amrex::RealVect(vel_table[lc][1], vel_table[lc][2], vel_table[lc][3]);
                } else {
                  //amrex::Print() << "Interpolating from values at entry " << lc << "\n";
                  const amrex::Real whi = (time - prv_time) / delta_time;
                  const amrex::Real wlo = 1.0 - whi;
                  return amrex::RealVect( whi*vel_table[lc][1] + wlo*vel_table[lc-1][1],
                                          whi*vel_table[lc][2] + wlo*vel_table[lc-1][2],
                                          whi*vel_table[lc][3] + wlo*vel_table[lc-1][3]);
                }
              }
            }
            // If we made it here, time is larger than all the entries in the
            // table, therefore, use the last value.
            //amrex::Print() << "Using last table entry\n";
            return amrex::RealVect(vel_table[size-1][1], vel_table[size-1][2], vel_table[size-1][3]);
          }
        }
      }

      bool constant_temperature;
      amrex::Real temperature;
      amrex::Vector<amrex::Vector<amrex::Real>> tg_table;

      amrex::Real get_temperature (amrex::Real time=0.0) {
        //amrex::Print() << "\n\nLooking for the right BC fluid temperature!\n\n";
        if (constant_temperature){
          // amrex::Print() << "* Use the constant temperature " << temperature << "\n";
          return temperature;
        } else {
          const int size = tg_table.size();
          if( size == 1 || time <= tg_table[0][0] ) {
            // There is only one entry (use it) OR the current time
            // is less than the first entry so return the first entry.
            // amrex::Print() << "Defaulting to the first entry " << tg_table[0][1] << "\n";
            return tg_table[0][1];

          } else {

            constexpr amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();
            //amrex::Print() << "Interpolating from the temperature table (maybe?)\n";
            for (int lc=1; lc<size; lc++) {
              // Time of current entry
              const amrex::Real entry_time = tg_table[lc][0];

              if ( amrex::Math::abs(time - entry_time) <= tolerance ) {
                // Too close for any calculations
                // amrex::Print() << "Too close to tell using value form entry "
                //                << lc << "  " << tg_table[lc][1] << "\n";
                return tg_table[lc][1];

                // The entries are sorted based on time, so if the current time is
                // less than this entries time, then we need to interpolate between
                // the previous time and this time.
              } else if (time < entry_time) {
                const amrex::Real prv_time = tg_table[lc-1][0];
                const amrex::Real delta_time = entry_time - prv_time;
                if (amrex::Math::abs(delta_time) <= tolerance) {
                  // amrex::Print() << "Delta time too small. Using value form entry "
                  //                << lc << "  " << tg_table[lc][1] << "\n";
                  return tg_table[lc][1];
                } else {
                  const amrex::Real whi = (time - prv_time) / delta_time;
                  const amrex::Real wlo = 1.0 - whi;
                  // amrex::Print() << "Interpolating from table entry " << lc << "\n";
                  return (whi*tg_table[lc][1] + wlo*tg_table[lc-1][1]);
                }
              }
            }
            // If we made it here, time is larger than all the entries in the
            // table, therefore, use the last value.
            // amrex::Print() << "Using last table entry " << tg_table[size-1][1];
            return tg_table[size-1][1];
          }
        }
      }

      FLUID_t() :
        volfrac(-1.0),
        pressure(-1.0),
        species(0),
        constant_velocity(true),
        constant_temperature(true),
        temperature(0.0)
      {}

    };

    void Initialize ();

}; // class FluidPhase

#endif
