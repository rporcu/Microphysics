#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>

#include <mfix_solids.H>
#include <mfix_species.H>
#include <mfix_fluid.H>
#include <mfix_algorithm.H>


using namespace amrex;

MFIXSolidsPhase::MFIXSolidsPhase()
  : m_ntypes(0)
  , m_specific_heat_model(MFIXSpecies::SpecificHeatModel::Invalid)
  , m_thermal_conductivity_model(ThermalConductivityModel::Invalid)
  , m_names(0)
  , m_update_mass(1)
  , m_solve_species(0)
  , m_update_momentum(1)
  , m_T_ref(0)
  , m_enthalpy_source(0)
  , m_update_enthalpy(1)
  , m_solve_enthalpy(0)
  , m_species_names(0)
  , m_species_IDs(0)
  , m_d_species_IDs(0)
  , m_nspecies(0)
  , m_MW_sn0(0)
  , m_d_MW_sn0(0)
  , m_is_a_mixture(0)
  , m_cp_sn0(0)
  , m_d_cp_sn0(0)
  , m_H_fn0(0)
  , m_d_H_fn0(0)
  , m_stoich_coeffs(0)
  , m_d_stoich_coeffs(0)
  , m_kp_sn0(0)
  , m_d_kp_sn0(0)
  , m_flpc(0.4)
  , m_rough(2.E-8)
  , m_do_conduction(0)
  , m_parameters(nullptr)
  , m_is_initialized(0)
{}


MFIXSolidsPhase::~MFIXSolidsPhase()
{
  if (m_parameters != nullptr)
    delete m_parameters;
}


int
MFIXSolidsPhase::name_to_phase (const std::string& name) const
{
  if (m_ntypes == 0) {
    Print() << "Can't look for solid type. No solids types found\n";
    amrex::Abort("Error. Fix inputs");
  }

  for (int n(0); n < m_ntypes; ++n) {
    if (m_names[n].compare(name) == 0) {
      return m_phases[n];
    }
  }

  amrex::Abort("Solid name provided does not exist");

  return -1;
}


void
MFIXSolidsPhase::Initialize (const MFIXSpecies& species,
                             const MFIXReactions& reactions)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species.isInitialized(),
      "Species not initialized. Can't initialize solids phase before species initialization");

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(reactions.isInitialized(),
      "MFIXReactions not initialized. Can't initialize solids phase before reactions initialization");

  // Flag for initialization
  m_is_initialized = 1;

  // Flag for pfp conduction
  int kp_type_read(0), kp_read(0), flpc_read(0), rough_read(0);

  amrex::ParmParse pp("solids");

  pp.queryarr("types", m_names);
  m_ntypes = m_names.size();

  // Disable the solids solver if the solids is defined as "" or "None"
  if ((m_ntypes > 0) && (amrex::toLower(m_names[0]).compare("none") != 0)) {

    // Resize phases and assign integer IDs to each phase as {1,2,...,ntypes+1}
    m_phases.resize(m_ntypes);
    for (int n(0); n < m_ntypes; ++n)
      m_phases[n] = n+1;

    // Get mfix global inputs ------------------------------------------------//
    amrex::ParmParse ppMFIX("mfix");

    ppMFIX.query("advect_enthalpy", m_solve_enthalpy);

    ppMFIX.query("solve_species", m_solve_species);

    // Query updating flags for particles mass, momentum and enthalpy
    {
      amrex::ParmParse pp_solids("mfix.particles");

      pp_solids.query("update_mass", m_update_mass);
      m_update_mass = m_update_mass && m_solve_species;

      pp_solids.query("update_momentum", m_update_momentum);

      pp_solids.query("enthalpy_source", m_enthalpy_source);

      pp_solids.query("update_enthalpy", m_update_enthalpy);
      m_update_enthalpy = m_update_enthalpy && m_solve_enthalpy;
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
                                       "No input provided for solids.species");

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_nspecies <= species.nspecies(),
          "Solids species number is higher than total species number");
    }

    // Flag to determine if we want to solve the solid as a mixture
    m_is_a_mixture = static_cast<int>(m_nspecies > 1);

    if (m_solve_enthalpy) {

      // Query the reference temperature
      pp.query("reference_temperature", m_T_ref);

      if (m_solve_species) {

        {
          std::string model;
          pp.get("specific_heat", model);

          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(amrex::toLower(model).compare("mixture") == 0,
              "When solving solids enthalpy and species, solids specific heat model must be a mixture");
        }

        m_specific_heat_model = species.specific_heat_model();

        if (m_specific_heat_model == MFIXSpecies::SpecificHeatModel::Constant) {
          m_cp_sn0.resize(m_nspecies);

        } else if (m_specific_heat_model == MFIXSpecies::SpecificHeatModel::NASA7Polynomials) {
          m_cp_sn0.resize(m_nspecies*12);
        }

        m_H_fn0.resize(m_nspecies);

        for (int n(0); n < m_nspecies; n++) {
          const auto& names = species.names();
          auto it = std::find(names.begin(), names.end(), m_species_names[n]);

          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != names.end(),
                                           "Solid species missing in input");

          const auto pos = std::distance(names.begin(), it);

          if (m_specific_heat_model == MFIXSpecies::SpecificHeatModel::Constant) {
            m_cp_sn0[n] = species.cp_k(pos);

          } else if (m_specific_heat_model == MFIXSpecies::SpecificHeatModel::NASA7Polynomials) {
            const auto& cp_k = species.cp_k();
            std::copy(&cp_k[pos*12], &cp_k[pos*12] + 12, &m_cp_sn0[n*12]);
          }

          m_H_fn0[n] = species.H_fk(pos);
        }

      } else {

        m_H_fn0.resize(1);

        // Get specific heat model input ------------------------//
        {
          std::string model;
          pp.get("specific_heat", model);

          if (amrex::toLower(model).compare("constant") == 0) {
            m_specific_heat_model = MFIXSpecies::SpecificHeatModel::Constant;
            m_cp_sn0.resize(1);
            pp.get("specific_heat.constant", m_cp_sn0[0]);

            // Get enthalpy of formation model input ------------------------//
            pp.query("enthalpy_of_formation", m_H_fn0[0]);

          } else if (amrex::toLower(model).compare("nasa7-poly") == 0) {

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_MW_sn0.size() > 0, "No solids molecular weight input");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_MW_sn0[0] > 0., "Wrong solids molecular weight");
            const amrex::Real coeff = MFIXFluidPhase::R / m_MW_sn0[0];

            m_specific_heat_model = MFIXSpecies::SpecificHeatModel::NASA7Polynomials;
            m_cp_sn0.resize(12);

            for (int i(0); i < 6; ++i) {
              amrex::Vector<amrex::Real> aa(0);
              std::string field = "specific_heat.NASA7.a"+std::to_string(i);
              pp.getarr(field.c_str(), aa);

              if (aa.size() == 1) {

                m_cp_sn0[i] = aa[0] * coeff;
                m_cp_sn0[i+6] = 0.;

              } else if (aa.size() == 2) {

                m_cp_sn0[i] = aa[0] * coeff;
                m_cp_sn0[i+6] = aa[1] * coeff;

              } else {

                amrex::Abort("Input error");
              }
            }

          } else {
            amrex::Abort("Don't know this specific heat model!");
          }
        }
      }
    }

    // Get solids species parameters from species class
    if (m_solve_species) {

      m_species_IDs.resize(m_nspecies);
      m_MW_sn0.resize(m_nspecies);

      for (int n(0); n < m_nspecies; n++) {
        const auto& names = species.names();
        auto it = std::find(names.begin(), names.end(), m_species_names[n]);

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != names.end(),
                                         "Solids species missing in input");

        const auto pos = std::distance(names.begin(), it);

        m_species_IDs[n] = species.IDs(pos);
        m_MW_sn0[n] = species.MW_k(pos);
      }

    } else {

      if (pp.contains("molecular_weight")) {
        m_MW_sn0.resize(1);
        pp.get("molecular_weight", m_MW_sn0[0]);
      }
    }

    // Create the stoichiometric table for the fluid phase, that is associate
    // the total stoichiometric coefficient for each fluid species in each
    // reaction
    if (reactions.solve()) {

      const int nreactions = reactions.nreactions();

      constexpr int Solid = MFIXChemicalReaction::ChemicalPhase::Solid;
      constexpr int Heterogeneous = MFIXChemicalReaction::ReactionType::Heterogeneous;

      // Allocate space for necessary data
      m_stoich_coeffs.resize(m_nspecies*nreactions, 0.);

      for (int n_s(0); n_s < m_nspecies; n_s++) {
        // Get the ID of the current species n_s
        const int species_id = m_species_IDs[n_s];

        // Loop over reactions to compute each contribution
        for (int q(0); q < nreactions; q++) {

          MFIXChemicalReaction* chem_reaction = reactions.get(q);
          const auto& chem_phases = chem_reaction->get_phases();

          const auto& reactants_IDs = chem_reaction->get_reactants_ids();
          const auto& reactants_chem_phases = chem_reaction->get_reactants_phases();
          const auto& reactants_coeffs = chem_reaction->get_reactants_coeffs();

          const auto& products_IDs = chem_reaction->get_products_ids();
          const auto& products_chem_phases = chem_reaction->get_products_phases();
          const auto& products_coeffs = chem_reaction->get_products_coeffs();

          // Do something only if reaction is heterogeneous and contains
          // a solid compound
          if (chem_reaction->get_type() == Heterogeneous &&
              std::find(chem_phases.begin(), chem_phases.end(), Solid) != chem_phases.end()) {

            // Add reactant contribution (if any)
            {
              for (int pos(0); pos < reactants_IDs.size(); ++pos) {
                if (species_id == reactants_IDs[pos] && reactants_chem_phases[pos] == Solid) {
                  m_stoich_coeffs[n_s*nreactions+q] += reactants_coeffs[pos];
                }
              }
            }

            // Add products contribution (if any)
            {
              for (int pos(0); pos < products_IDs.size(); ++pos) {
                if (species_id == products_IDs[pos] && products_chem_phases[pos] == Solid) {
                  m_stoich_coeffs[n_s*nreactions+q] += products_coeffs[pos];
                }
              }
            }
          }
        }
      }
    }

    // Do not support multiple species for solids conductivities!
    // Always reads in the conductivity model from first solids phase
    std::string model;
    kp_type_read = pp.query("thermal_conductivity", model);

    if (amrex::toLower(model).compare("constant") == 0) {
        m_thermal_conductivity_model = ThermalConductivityModel::Constant;
        m_kp_sn0.resize(1);
        kp_read = pp.query("thermal_conductivity.constant", m_kp_sn0[0]);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_kp_sn0[0] > 0, "Invalid DEM thermal conductivity.");
    } /*else {
        amrex::Abort("Unknown particle conductivity model!");
    }*/


    // FLPC and roughness are always dimension 1
    if (m_solve_enthalpy) {
        flpc_read  = pp.query("flpc", m_flpc);
        rough_read = pp.query("min_conduction_dist", m_rough);
    }

    if(kp_type_read && kp_read && flpc_read && rough_read)
      m_do_conduction = 1;

    // Allocate parameters
    m_parameters = new MFIXSolidsParms(m_T_ref);

    // Species IDs
    if (m_species_IDs.size() > 0) {
      m_d_species_IDs.resize(m_species_IDs.size());
      Gpu::copyAsync(Gpu::hostToDevice, m_species_IDs.begin(), m_species_IDs.end(), m_d_species_IDs.begin());
      m_parameters->m_h_species_id = m_species_IDs.data();
      m_parameters->m_d_species_id = m_d_species_IDs.data();
    }

    // Molecular weights
    if (m_MW_sn0.size() > 0) {
      m_d_MW_sn0.resize(m_MW_sn0.size());
      Gpu::copyAsync(Gpu::hostToDevice, m_MW_sn0.begin(), m_MW_sn0.end(), m_d_MW_sn0.begin());
      m_parameters->m_h_MW_sn = m_MW_sn0.data();
      m_parameters->m_d_MW_sn = m_d_MW_sn0.data();
    }

    // Specific heat coefficients
    if (m_cp_sn0.size() > 0) {
      m_d_cp_sn0.resize(m_cp_sn0.size());
      Gpu::copyAsync(Gpu::hostToDevice, m_cp_sn0.begin(), m_cp_sn0.end(), m_d_cp_sn0.begin());
      m_parameters->m_h_cp_sn = m_cp_sn0.data();
      m_parameters->m_d_cp_sn = m_d_cp_sn0.data();
    }

    // Enthalpies of formation
    if (m_H_fn0.size() > 0) {
      m_d_H_fn0.resize(m_H_fn0.size());
      Gpu::copyAsync(Gpu::hostToDevice, m_H_fn0.begin(), m_H_fn0.end(), m_d_H_fn0.begin());
      m_parameters->m_h_H_fn = m_H_fn0.data();
      m_parameters->m_d_H_fn = m_d_H_fn0.data();
    }

    m_parameters->m_nreactions = reactions.nreactions();

    // Stoichiometric coefficients
    if (m_stoich_coeffs.size() > 0) {
      m_d_stoich_coeffs.resize(m_stoich_coeffs.size());
      Gpu::copyAsync(Gpu::hostToDevice, m_stoich_coeffs.begin(), m_stoich_coeffs.end(), m_d_stoich_coeffs.begin());
      m_parameters->m_h_stoich_coeffs = m_stoich_coeffs.data();
      m_parameters->m_d_stoich_coeffs = m_d_stoich_coeffs.data();
    }

    //
    if (m_kp_sn0.size() > 0) {
      m_d_kp_sn0.resize(m_kp_sn0.size());
      Gpu::copyAsync(Gpu::hostToDevice, m_kp_sn0.begin(), m_kp_sn0.end(), m_d_kp_sn0.begin());
      m_parameters->m_h_kp_sn = m_kp_sn0.data();
      m_parameters->m_d_kp_sn = m_d_kp_sn0.data();
    }

    m_parameters->m_flpc = m_flpc;
    m_parameters->m_rough = m_rough;
    m_parameters->m_do_pfp_cond = m_do_conduction;

    m_parameters->m_specific_heat_model = m_specific_heat_model;
    m_parameters->m_thermal_conductivity_model = m_thermal_conductivity_model;
  }
}


int
INPUT_DIST_t::set (const std::string a_field,
                   const std::string a_pprop)
{
  amrex::ParmParse pp(a_field.c_str());

  std::string field_pprop = a_field  + "." + a_pprop;
  amrex::ParmParse ppDist(field_pprop.c_str());

  // Get the distribution type.
  std::string distribution;
  pp.get(a_pprop.c_str(), distribution);

  int err = 0;
  int fnd_mean(1), fnd_std(1), fnd_min(1), fnd_max(1);

  if( amrex::toLower(distribution).compare("constant") == 0) {
    m_is_constant = 1;
    fnd_mean = ppDist.query("constant", m_mean);
    m_stddev = 0.;
    m_max = m_mean;
    m_min = m_mean;

  } else if( amrex::toLower(distribution).compare("normal") == 0) {
    m_is_normal = 1;

    fnd_mean = ppDist.query("mean", m_mean);
    fnd_std  = ppDist.query("std",  m_stddev);
    fnd_min  = ppDist.query("min",  m_min);
    fnd_max  = ppDist.query("max",  m_max);

    if(m_max < m_min || m_mean < m_min || m_max < m_mean) {
      amrex::Print() << "Invalid uniform distribution parameters!\n";
      err = 1;
    }

  } else if( amrex::toLower(distribution).compare("uniform") == 0) {
    m_is_uniform = 1;
    fnd_min  = ppDist.query("min", m_min);
    fnd_max  = ppDist.query("max", m_max);
    m_mean = m_min + 0.5*(m_max-m_min);
    m_stddev = 0.;

    if(m_max < m_min) {
      amrex::Print() << "Invalid normal distribution parameters!\n";
      err = 1;
    }

  } else {
    amrex::Print() << "Invalid distribution type: " + distribution + "!\n";
    return 1;
  }

  // This lets the calling routine know that something is missing (which
  // was already reported) so that information about the IC/BC region
  // and solid can be printed and the run can be aborted.
  if( !fnd_mean || !fnd_std || !fnd_min || !fnd_max ) {
    err = 1;
  }
  return err;
}
