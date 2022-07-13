#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>

#include <mfix_fluid.H>
#include <mfix_species.H>
#include <mfix_algorithm.H>


using namespace amrex;


AMREX_GPU_HOST_DEVICE
MFIXFluidParms::MFIXFluidParms (const amrex::Real T_ref)
  : m_T_ref(T_ref)
  , m_mu_g(0.)
  , m_k_g(0.)
  , m_nspecies(0)
  , m_h_species_id(nullptr)
  , m_d_species_id(nullptr)
  , m_h_MW_gk(nullptr)
  , m_d_MW_gk(nullptr)
  , m_D_g(0)
  , m_h_cp_gk(nullptr)
  , m_d_cp_gk(nullptr)
  , m_h_H_fk(nullptr)
  , m_d_H_fk(nullptr)
  , m_nreactions(0)
  , m_h_stoich_coeffs(nullptr)
  , m_d_stoich_coeffs(nullptr)
  , m_specific_heat_model(MFIXSpecies::SpecificHeatModel::Invalid)
{}


MFIXFluidPhase::MFIXFluidPhase ()
  : m_viscosity_model(ViscosityModel::Invalid)
  , m_specific_heat_model(MFIXSpecies::SpecificHeatModel::Invalid)
  , m_thermal_conductivity_model(ThermalConductivityModel::Invalid)
  , m_names(0)
  , m_solve(0)
  , m_solve_density(0)
  , m_solve_tracer(0)
  , m_trac_0(0)
  , m_mu_g0(0)
  , m_mw_avg(0)
  , m_solve_enthalpy(0)
  , m_T_g0(273.15)
  , m_T_ref(0)
  , m_k_g0(0)
  , m_thermodynamic_pressure(-1.0)
  , m_p_therm_defined(0)
  , m_solve_species(0)
  , m_species_names(0)
  , m_species_IDs(0)
  , m_d_species_IDs(0)
  , m_nspecies(0)
  , m_MW_gk0(0)
  , m_d_MW_gk0(0)
  , m_D_g0(0)
  , m_is_a_mixture(0)
  , m_H_fk0(0)
  , m_d_H_fk0(0)
  , m_cp_gk0(0)
  , m_d_cp_gk0(0)
  , m_stoich_coeffs(0)
  , m_d_stoich_coeffs(0)
  , m_parameters(nullptr)
  , m_is_initialized(0)
{}

MFIXFluidPhase::~MFIXFluidPhase()
{
  if (m_parameters != nullptr)
    delete m_parameters;
}

void
MFIXFluidPhase::Initialize (const MFIXSpecies& species,
                            const MFIXReactions& reactions)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species.isInitialized(),
      "Species not initialized. Can't initialize fluid phase before species initialization");

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(reactions.isInitialized(),
      "MFIXReactions not initialized. Can't initialize fluid phase before reactions initialization");

  // Flag for initialization
  m_is_initialized = 1;

  amrex::ParmParse pp("fluid");

  pp.queryarr("solve", m_names);
  m_ntypes = m_names.size();

  // Abort if more than one fluid type because it's not yet implemented
  if (m_ntypes > 1) {
    amrex::Abort("Handling more than one fluid types has not been implemented yet");
  }

  // Disable the fluid solver if the fluid is defined as "" or "None"
  if ((m_ntypes > 0) && (amrex::toLower(m_names[0]).compare("none") != 0)) {

    // Set the fluid solving flag to true
    m_solve = 1;

    // Get mfix global inputs ------------------------------------------------//
    amrex::ParmParse ppMFIX("mfix");

    ppMFIX.query("advect_density", m_solve_density);

    ppMFIX.query("advect_enthalpy", m_solve_enthalpy);

    ppMFIX.query("advect_tracer", m_solve_tracer);

    ppMFIX.query("solve_species", m_solve_species);

    // Query molecular weight average value
    pp.query("mw_avg", m_mw_avg);

    // Get viscosity inputs ----------------------------------//
    {
      std::string model;
      pp.get("viscosity", model);

      if (amrex::toLower(model).compare("constant") == 0) {
        m_viscosity_model = ViscosityModel::Constant;
        pp.get("viscosity.constant", m_mu_g0);

      } else if (amrex::toLower(model).compare("sutherland") == 0) {
        m_viscosity_model = ViscosityModel::Sutherland;
        amrex::Abort("Not yet implemented.");

      } else {
        amrex::Abort("Unknown fluid viscosity model!");
      }
    }

    if (m_solve_tracer) {
      // Query fluid tracer initial value
      pp.query("trac0", m_trac_0);
    }

    if (m_solve_species) {
      // Query fluid species
      pp.queryarr("species", m_species_names);
    }

    if ((!m_solve_species) || (m_species_names.size() == 0)) {
      m_species_names.clear();
      m_solve_species = 0;
      m_nspecies = 0;
    } else if (amrex::toLower(m_species_names[0]).compare("none") == 0) {
      m_species_names.clear();
      m_solve_species = 0;
      m_nspecies = 0;
    } else {
      m_solve_species = 1;
      m_nspecies = m_species_names.size();

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_species_names.size() > 0, 
                                       "No input provided for fluid.names_names");

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_nspecies <= species.nspecies(),
          "Fluid species number is higher than total species number");
    }

    // Flag to determine if we want to solve the fluid as a mixture
    m_is_a_mixture = static_cast<int>(m_nspecies > 1);

    if (m_solve_enthalpy) {

      //  Query fluid clod-flow temperature
      pp.query("T_g0", m_T_g0);

      // Query the reference temperature
      pp.query("reference_temperature", m_T_ref);

      // Get thermal conductivity inputs -----------------------------//
      {
        std::string model;
        pp.get("thermal_conductivity", model);

        if (amrex::toLower(model).compare("constant") == 0) {
          m_thermal_conductivity_model = ThermalConductivityModel::Constant;
          pp.get("thermal_conductivity.constant", m_k_g0);

        } else {
          amrex::Abort("Unknown fluid thermal conductivity model!");
        }
      }

      if (m_solve_species) {
        // Specific heat model from inputs has to be a mixture
        {
          std::string model;
          pp.get("specific_heat", model);

          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(amrex::toLower(model).compare("mixture") == 0,
              "When solving fluid enthalpy and species, fluid specific heat model must be a mixture");
        }

        m_specific_heat_model = species.specific_heat_model();

        if (m_specific_heat_model == MFIXSpecies::SpecificHeatModel::Constant) {
          m_cp_gk0.resize(m_nspecies);

        } else if (m_specific_heat_model == MFIXSpecies::SpecificHeatModel::NASA7Polynomials) {
          m_cp_gk0.resize(m_nspecies*12);
        }

        m_H_fk0.resize(m_nspecies);
      
        for (int n(0); n < m_nspecies; n++) {
          const auto& names = species.names();
          auto it = std::find(names.begin(), names.end(), m_species_names[n]);

          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != names.end(),
                                           "Fluid species missing in input");

          const auto pos = std::distance(names.begin(), it);

          if (m_specific_heat_model == MFIXSpecies::SpecificHeatModel::Constant) {
            m_cp_gk0[n] = species.cp_k(pos);

          } else if (m_specific_heat_model == MFIXSpecies::SpecificHeatModel::NASA7Polynomials) {
            const auto& cp_k = species.cp_k();
            std::copy(&cp_k[pos*12], &cp_k[pos*12] + 12, &m_cp_gk0[n*12]);
          }

          m_H_fk0[n] = species.H_fk(pos);
        }
      
      } else {

        m_H_fk0.resize(1);

        // Get specific heat inputs ------------------------------------//
        {
          std::string model;
          pp.get("specific_heat", model);

          if (amrex::toLower(model).compare("constant") == 0) {
            m_specific_heat_model = MFIXSpecies::SpecificHeatModel::Constant;
            m_cp_gk0.resize(1);
            pp.get("specific_heat.constant", m_cp_gk0[0]);

            // Query the enthalpy_of_formation
            pp.query("enthalpy_of_formation", m_H_fk0[0]);

          } else if (amrex::toLower(model).compare("nasa7-poly") == 0) {

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_MW_gk0.size() > 0, "No fluid molecular weight input");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_MW_gk0[0] > 0., "Wrong fluid molecular weight");
            const amrex::Real coeff = R / m_MW_gk0[0];

            m_specific_heat_model = MFIXSpecies::SpecificHeatModel::NASA7Polynomials;
            m_cp_gk0.resize(12);

            for (int i(0); i < 6; ++i) {
              amrex::Vector<amrex::Real> aa(0);
              std::string field = "specific_heat.NASA7.a"+std::to_string(i);
              pp.getarr(field.c_str(), aa);

              if (aa.size() == 1) {

                m_cp_gk0[i] = aa[0] * coeff;
                m_cp_gk0[i+6] = 0.;

              } else if (aa.size() == 2) {

                m_cp_gk0[i] = aa[0] * coeff;
                m_cp_gk0[i+6] = aa[1] * coeff;

              } else {

                amrex::Abort("Input error");
              }
            }

          } else {
            amrex::Abort("Unknown fluid specific heat model!");
          }
        }
      }
    }

    // Get fluid species m_parameters from species class
    if (m_solve_species) {

      m_species_IDs.resize(m_nspecies);
      m_MW_gk0.resize(m_nspecies);
      m_D_g0 = species.diffusivity();

      for (int n(0); n < m_nspecies; n++) {
        const auto& names = species.names();
        auto it = std::find(names.begin(), names.end(), m_species_names[n]);

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != names.end(),
                                         "Fluid species missing in input");

        const auto pos = std::distance(names.begin(), it);

        m_species_IDs[n] = species.IDs(pos);
        m_MW_gk0[n] = species.MW_k(pos);
      }

    } else {

      if (pp.contains("molecular_weight")) {
        m_MW_gk0.resize(1);
        pp.get("molecular_weight", m_MW_gk0[0]);
      }
    }

    // Create the stoichiometric table for the fluid phase, that is associate
    // the total stoichiometric coefficient for each fluid species in each
    // reaction
    if (reactions.solve()) {

      const int m_nreactions = reactions.nreactions();

      constexpr int Fluid = MFIXChemicalReaction::ChemicalPhase::Fluid;
      constexpr int Heterogeneous = MFIXChemicalReaction::ReactionType::Heterogeneous;

      // Allocate space for necessary data
      m_stoich_coeffs.resize(m_nspecies*m_nreactions, 0.);

      for (int n_g(0); n_g < m_nspecies; n_g++) {
        // Get the ID of the current species n_g
        const int species_id = m_species_IDs[n_g];

        // Loop over reactions to compute each contribution
        for (int q(0); q < m_nreactions; q++) {

          MFIXChemicalReaction* chem_reaction = reactions.get(q);
          const auto& phases = chem_reaction->get_phases();

          const auto& reactants_IDs = chem_reaction->get_reactants_ids();
          const auto& reactants_phases = chem_reaction->get_reactants_phases();
          const auto& reactants_coeffs = chem_reaction->get_reactants_coeffs();

          const auto& products_IDs = chem_reaction->get_products_ids();
          const auto& products_phases = chem_reaction->get_products_phases();
          const auto& products_coeffs = chem_reaction->get_products_coeffs();

          // Do something only if reaction is heterogeneous and contains
          // a solid compound
          if (chem_reaction->get_type() == Heterogeneous &&
              std::find(phases.begin(), phases.end(), Fluid) != phases.end()) {

            // Add reactant contribution (if any)
            for (int pos(0); pos < reactants_IDs.size(); ++pos) {
              if (species_id == reactants_IDs[pos] &&
                  reactants_phases[pos] == Fluid) {
                m_stoich_coeffs[n_g*m_nreactions+q] += reactants_coeffs[pos];
              }
            }

            // Add products contribution (if any)
            for (int pos(0); pos < products_IDs.size(); ++pos) {
              if (species_id == products_IDs[pos] &&
                  products_phases[pos] == Fluid) {
                m_stoich_coeffs[n_g*m_nreactions+q] += products_coeffs[pos];
              }
            }
          }
        }
      }
    }

    // Constraint type
    {
      ParmParse ppMFIX("mfix");

      std::string constraint_str = "IncompressibleFluid";
      ppMFIX.query("constraint_type", constraint_str);
      constraint_str = amrex::toLower(constraint_str);

      if (constraint_str.compare("incompressiblefluid") == 0) {
        m_constraint_type = ConstraintType::IncompressibleFluid;
      }
      else if (constraint_str.compare("idealgasopensystem") == 0) {
        m_constraint_type = ConstraintType::IdealGasOpenSystem;
      }
      else if (constraint_str.compare("idealgasclosedsystem") == 0) {
        m_constraint_type = ConstraintType::IdealGasClosedSystem;
      }
      else {
        amrex::Abort("Don't know this constraint type!");
      }
    }

    m_p_therm_defined = pp.query("thermodynamic_pressure", m_thermodynamic_pressure);

    if (m_constraint_type == ConstraintType::IncompressibleFluid && m_p_therm_defined) {

      amrex::Warning("When the incompressible fluid constraint is selected, "
          "the fluid thermodynamic pressure input will be ignored");
    }

    if (m_constraint_type == ConstraintType::IdealGasClosedSystem && (!m_p_therm_defined)) {

      amrex::Abort("When the idealgas closedsystem constraint is selected, "
          "the fluid thermodynamic pressure input must be provided");
    }

    // Check on inputs in case of Ideal Gas EOS
    if (m_constraint_type == ConstraintType::IdealGasOpenSystem ||
        m_constraint_type == ConstraintType::IdealGasClosedSystem) {
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_MW_gk0.size() > 0, "Inputs error: fluid molecular_weight not provided");

      for (size_t i(0); i < m_MW_gk0.size(); ++i) {

        if (m_MW_gk0[i] < 1.e-15) {

          Print() << "Invalid molecular weight for species " << m_species_names[i] << "\n";
          amrex::Abort("Inputs error");
        }
      }
    }

    // Allocate m_parameters
    m_parameters = new MFIXFluidParms(m_T_ref);

    m_parameters->m_mu_g = m_mu_g0;
    m_parameters->m_k_g = m_k_g0;
    m_parameters->m_nspecies = m_nspecies;

    if (m_species_IDs.size() > 0) {
      m_d_species_IDs.resize(m_species_IDs.size());
      Gpu::copyAsync(Gpu::hostToDevice, m_species_IDs.begin(), m_species_IDs.end(), m_d_species_IDs.begin());
      m_parameters->m_h_species_id = m_species_IDs.data();
      m_parameters->m_d_species_id = m_d_species_IDs.data();
    }

    if (m_MW_gk0.size() > 0) {
      m_d_MW_gk0.resize(m_MW_gk0.size());
      Gpu::copyAsync(Gpu::hostToDevice, m_MW_gk0.begin(), m_MW_gk0.end(), m_d_MW_gk0.begin());
      m_parameters->m_h_MW_gk = m_MW_gk0.data();
      m_parameters->m_d_MW_gk = m_d_MW_gk0.data();
    }

    m_parameters->m_D_g = m_D_g0;

    if (m_cp_gk0.size() > 0) {
      m_d_cp_gk0.resize(m_cp_gk0.size());
      Gpu::copyAsync(Gpu::hostToDevice, m_cp_gk0.begin(), m_cp_gk0.end(), m_d_cp_gk0.begin());
      m_parameters->m_h_cp_gk = m_cp_gk0.data();
      m_parameters->m_d_cp_gk = m_d_cp_gk0.data();
    }

    if (m_H_fk0.size() > 0) {
      m_d_H_fk0.resize(m_H_fk0.size());
      Gpu::copyAsync(Gpu::hostToDevice, m_H_fk0.begin(), m_H_fk0.end(), m_d_H_fk0.begin());
      m_parameters->m_h_H_fk = m_H_fk0.data();
      m_parameters->m_d_H_fk = m_d_H_fk0.data();
    }

    m_parameters->m_nreactions = reactions.nreactions();

    if (m_stoich_coeffs.size() > 0) {
      m_d_stoich_coeffs.resize(m_stoich_coeffs.size());
      Gpu::copyAsync(Gpu::hostToDevice, m_stoich_coeffs.begin(), m_stoich_coeffs.end(), m_d_stoich_coeffs.begin());
      m_parameters->m_h_stoich_coeffs = m_stoich_coeffs.data();
      m_parameters->m_d_stoich_coeffs = m_d_stoich_coeffs.data();
    }

    m_parameters->m_specific_heat_model = m_specific_heat_model;
  }
}


Real
FLUID_t::get_density (amrex::Real time) const
{
  if (density_defined) {
    //amrex::Print() << "\n\nLooking for the right BC fluid temperature!\n\n";
    if (constant_density) {
      // amrex::Print() << "* Use the constant temperature " << temperature << "\n";
      return density;
    } else {
      const int size = density_table.size();
      if( size == 1 || time <= density_table[0][0] ) {
        // There is only one entry (use it) OR the current time
        // is less than the first entry so return the first entry.
        // amrex::Print() << "Defaulting to the first entry " << density_table[0][1] << "\n";
        return density_table[0][1];

      } else {

        constexpr amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();
        //amrex::Print() << "Interpolating from the temperature table (maybe?)\n";
        for (int lc=1; lc<size; lc++) {
          // Time of current entry
          const amrex::Real entry_time = density_table[lc][0];

          if ( amrex::Math::abs(time - entry_time) <= tolerance ) {
            // Too close for any calculations
            // amrex::Print() << "Too close to tell using value form entry "
            //                << lc << "  " << density_table[lc][1] << "\n";
            return density_table[lc][1];

            // The entries are sorted based on time, so if the current time is
            // less than this entries time, then we need to interpolate between
            // the previous time and this time.
          } else if (time < entry_time) {
            const amrex::Real prv_time = density_table[lc-1][0];
            const amrex::Real delta_time = entry_time - prv_time;
            if (amrex::Math::abs(delta_time) <= tolerance) {
              // amrex::Print() << "Delta time too small. Using value form entry "
              //                << lc << "  " << density_table[lc][1] << "\n";
              return density_table[lc][1];
            } else {
              const amrex::Real whi = (time - prv_time) / delta_time;
              const amrex::Real wlo = 1.0 - whi;
              // amrex::Print() << "Interpolating from table entry " << lc << "\n";
              return (whi*density_table[lc][1] + wlo*density_table[lc-1][1]);
            }
          }
        }
        // If we made it here, time is larger than all the entries in the
        // table, therefore, use the last value.
        // amrex::Print() << "Using last table entry " << density_table[size-1][1];
        return density_table[size-1][1];
      }
    }
  } else { // density not defined

    if (fluid.constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasOpenSystem ||
        fluid.constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasClosedSystem) {

      const auto& fluid_parms = fluid.parameters();

      amrex::Real MW_g(0.);

      if (fluid.isMixture()) {
        for (int n(0); n < fluid.nspecies(); n++) {
          MW_g += this->get_species(n, time) / fluid_parms.get_MW_gk<amrex::RunOn::Host>(n);
        }

        MW_g = 1./MW_g;
      } else {
        MW_g = fluid_parms.get_MW_g<amrex::RunOn::Host>();
      }

      return (fluid.thermodynamic_pressure() * MW_g) /
             (MFIXFluidPhase::R * this->get_temperature(time));

    } else {

      amrex::Abort("Error. How did we even arrive here?");
      return 0;
    }
  }
}


Real
FLUID_t::get_species (int n, amrex::Real time) const
{
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


RealVect
FLUID_t::get_velocity (amrex::Real time) const
{
  //amrex::Print() << "\n\nLooking for the right BC fluid velocity!\n\n";
  if (constant_velocity) {
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


Real
FLUID_t::get_velocity_mag (amrex::Real time) const
{
  AMREX_ASSERT(flow_thru_eb);
  //amrex::Print() << "\n\nLooking for the right BC fluid velocity!\n\n";
  if (constant_velocity){
    //amrex::Print() << "* Use the constant velocity!";
    return velocity[0];
  } else {
    const int size = vel_table.size();
    if( size == 1 || time <= vel_table[0][0] ) {
      // There is only one entry (use it) OR the current time
      // is less than the first entry so return the first entry.
      //amrex::Print() << "Defaulting to the first entry!\n";
      return vel_table[0][1];

    } else {

      constexpr amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();
      //amrex::Print() << "Interpolating from the velocity table (maybe?)\n";
      for (int lc=1; lc<size; lc++) {
        // Time of current entry
        const amrex::Real entry_time = vel_table[lc][0];

        if ( amrex::Math::abs(time - entry_time) <= tolerance ) {
          // Too close for any calculations
          // amrex::Print() << "Too close to tell using value form entry " << lc << "\n";
          return vel_table[lc][1];

          // The entries are sorted based on time, so if the current time is
          // less than this entries time, then we need to interpolate between
          // the previous time and this time.
        } else if (time < entry_time) {
          const amrex::Real prv_time = vel_table[lc-1][0];
          const amrex::Real delta_time = entry_time - prv_time;
          if (amrex::Math::abs(delta_time) <= tolerance) {
            //amrex::Print() << "Delta time too small. Using value form entry " << lc << "\n";
            return vel_table[lc][1];
          } else {
            //amrex::Print() << "Interpolating from values at entry " << lc << "\n";
            const amrex::Real whi = (time - prv_time) / delta_time;
            const amrex::Real wlo = 1.0 - whi;
            return whi*vel_table[lc][1] + wlo*vel_table[lc-1][1];
          }
        }
      }
      // If we made it here, time is larger than all the entries in the
      // table, therefore, use the last value.
      //amrex::Print() << "Using last table entry\n";
      return vel_table[size-1][1];
    }
  }
}


Real
FLUID_t::get_volflow (amrex::Real time) const
{
  AMREX_ASSERT(flow_thru_eb);
  //amrex::Print() << "\n\nLooking for the right BC fluid velocity!\n\n";
  if (constant_volflow){
    //amrex::Print() << "* Use the constant velocity!";
    return volflow;
  } else {
    const int size = volflow_table.size();
    if( size == 1 || time <= volflow_table[0][0] ) {
      // There is only one entry (use it) OR the current time
      // is less than the first entry so return the first entry.
      //amrex::Print() << "Defaulting to the first entry!\n";
      return volflow_table[0][1];

    } else {

      constexpr amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();
      //amrex::Print() << "Interpolating from the velocity table (maybe?)\n";
      for (int lc=1; lc<size; lc++) {
        // Time of current entry
        const amrex::Real entry_time = volflow_table[lc][0];

        if ( amrex::Math::abs(time - entry_time) <= tolerance ) {
          // Too close for any calculations
          // amrex::Print() << "Too close to tell using value form entry " << lc << "\n";
          return volflow_table[lc][1];

          // The entries are sorted based on time, so if the current time is
          // less than this entries time, then we need to interpolate between
          // the previous time and this time.
        } else if (time < entry_time) {
          const amrex::Real prv_time = volflow_table[lc-1][0];
          const amrex::Real delta_time = entry_time - prv_time;
          if (amrex::Math::abs(delta_time) <= tolerance) {
            //amrex::Print() << "Delta time too small. Using value form entry " << lc << "\n";
            return volflow_table[lc][1];
          } else {
            //amrex::Print() << "Interpolating from values at entry " << lc << "\n";
            const amrex::Real whi = (time - prv_time) / delta_time;
            const amrex::Real wlo = 1.0 - whi;
            return whi*volflow_table[lc][1] + wlo*volflow_table[lc-1][1];
          }
        }
      }
      // If we made it here, time is larger than all the entries in the
      // table, therefore, use the last value.
      //amrex::Print() << "Using last table entry\n";
      return volflow_table[size-1][1];
    }
  }
}


Real
FLUID_t::get_temperature (amrex::Real time) const
{
  // Sanity check
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(temperature_defined, "Error: temperature BCs were not defined");

  //amrex::Print() << "\n\nLooking for the right BC fluid temperature!\n\n";
  if (constant_temperature) {
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
