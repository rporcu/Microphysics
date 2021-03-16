#ifndef MFIX_FLUID_PARMS_H_
#define MFIX_FLUID_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_GpuContainers.H>

#include <mfix_species_parms.H>

namespace FLUID
{
  // List of recognized fluid density models
  struct DENSITYMODEL { enum {Invalid=-1, Constant, IdealGas}; };

  extern int DensityModel;

  // List of recognized fluid molecular weight models
  struct MOLECULARWEIGHTMODEL { enum {Invalid=-1, Constant, Mixture}; };

  extern int MolecularWeightModel;

  // List of recognized fluid viscosity models
  struct VISCOSITYMODEL { enum {Invalid=-1, Constant, Sutherland}; };

  extern int ViscosityModel;

  // List of recognized fluid specific heat models
  struct SPECIFICHEATMODEL { enum {Invalid=-1, Constant, NASA9Polynomials}; };

  extern int SpecificHeatModel;

  // List of recognized fluid thermal conductivity models
  struct THERMALCONDUCTIVITYMODEL{ enum {Invalid=-1, Constant}; };

  extern int ThermalConductivityModel;

  static constexpr amrex::Real R = 8.31446261815324;

  // Flag to solve fluid equations
  extern int solve;

  // Specified constant gas density
  extern amrex::Real ro_g0;

  // Specified constant molecular weight of gas
  extern amrex::Real MW_g0;

  // Specified constant tracer value
  extern amrex::Real trac_0;

  // Specified constant gas viscosity
  extern amrex::Real mu_g0;

  // Average molecular weight of gas
  extern amrex::Real mw_avg;

  // Flag to solve enthalpy fluid equation
  extern int solve_enthalpy;

  // Specified constant gas temperature
  extern amrex::Real T_g0;

  // Specified constant gas specific heat
  extern amrex::Real cp_g0;

  // Specified constant gas phase thermal conductivity coefficient
  extern amrex::Real k_g0;

  // Flag to solve species fluid equations
  extern int solve_species;

  // Fluid phase species names
  extern amrex::Vector<std::string> species;

  // Species unique identifying code
  extern amrex::Gpu::HostVector<int> species_id;

  // Total number of fluid species
  extern int nspecies;

  // Specified constant gas phase species molecular weight
  extern amrex::Gpu::HostVector<amrex::Real> MW_gk0;

  // Specified constant gas phase species diffusion coefficients
  extern amrex::Gpu::HostVector<amrex::Real> D_gk0;

  // Flag to understand if fluid is a mixture of fluid species
  extern int is_a_mixture;

  // Specified constant gas phase species specific heat
  extern amrex::Gpu::HostVector<amrex::Real> cp_gk0;

  // Name to later reference when building input regions
  extern std::string name;

  struct FLUID_t {

    amrex::Real volfrac;
    amrex::Real pressure;
    amrex::Real pressure_defined; // Flag to check if pressure is defined

    amrex::Vector<SPECIES::SPECIES_t> species;

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

}

#endif
