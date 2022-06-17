#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>

#include <mfix_solids_parms.H>
#include <mfix_species_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_algorithm.H>


using namespace amrex;

SolidsPhase::SolidsPhase()
  : NTYPES(0)
  , specific_heat_model(SPECIFICHEATMODEL::Invalid)
  , thermal_conductivity_model(SOLIDSTHERMALCONDUCTIVITYMODEL::Invalid)
  , names(0)
  , solve_mass(1)
  , solve_momentum(1)
  , T_ref(0)
  , enthalpy_source(0)
  , solve_enthalpy(1)
  , solve_species(0)
  , species_names(0)
  , species_IDs(0)
  , d_species_IDs(0)
  , nspecies(0)
  , MW_sn0(0)
  , d_MW_sn0(0)
  , is_a_mixture(0)
  , cp_sn0(0)
  , d_cp_sn0(0)
  , H_fn0(0)
  , d_H_fn0(0)
  , stoich_coeffs(0)
  , d_stoich_coeffs(0)
  , kp_sn0(0)
  , d_kp_sn0(0)
  , flpc(0.4)
  , rough(2.E-8)
  , do_pfp_cond(0)
  , parameters(nullptr)
  , is_initialized(0)
{}


SolidsPhase::~SolidsPhase()
{
  if (parameters != nullptr)
    delete parameters;
}


int
SolidsPhase::name_to_phase (const std::string& name) const
{
  if (NTYPES == 0) {
    Print() << "Can't look for solid type. No solids types found\n";
    amrex::Abort("Error. Fix inputs");
  }

  for (int n(0); n < NTYPES; ++n) {
    if (names[n].compare(name) == 0) {
      return phases[n];
    }
  }

  amrex::Abort("Solid name provided does not exist");

  return -1;
}


void
SolidsPhase::Initialize (const Species& species,
                         const Reactions& reactions)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species.is_initialized,
      "Species not initialized. Can't initialize solids phase before species initialization");

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(reactions.is_initialized,
      "Reactions not initialized. Can't initialize solids phase before reactions initialization");

  // Flag for initialization
  is_initialized = 1;

  int solve_energy = 0;

  {
    amrex::ParmParse pp_mfix("mfix");
    pp_mfix.query("advect_enthalpy", solve_energy);

    {
      amrex::ParmParse pp_mfix_solids("mfix.particles");
      pp_mfix_solids.query("enthalpy_source", enthalpy_source);
      pp_mfix_solids.query("update_mass", solve_mass);
      pp_mfix_solids.query("update_momentum", solve_momentum);
      pp_mfix_solids.query("update_enthalpy", solve_enthalpy);
    }
  }

  // Flag for pfp conduction
  int kp_type_read(0), kp_read(0), flpc_read(0), rough_read(0);
  int do_conduction(0);

  amrex::ParmParse pp("solids");

  if(solve_energy) {
    // Query the reference temperature
    pp.query("reference_temperature", T_ref);
  }

  parameters = new SolidsParms(T_ref);

  pp.queryarr("types", names);

  if (names.size() > 0 && amrex::toLower(names[0]).compare("none") != 0) {

    NTYPES = names.size();

    phases.resize(NTYPES);
    for (int n(0); n < NTYPES; ++n)
      phases[n] = n+1;

    // Species
    solve_species = pp.queryarr("species", species_names);

    if (!solve_species) {
      species_names.clear();
      nspecies = 0;
    } else {
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species_names.size() > 0,
                                       "No input provided for solids.species");
    }

    // Solids species inputs
    if (solve_species) {

      // Disable the species solver if the species are defined as "None"
      // (caseinsensitive) or 0
      if (amrex::toLower(species_names[0]).compare("none") == 0) {
        solve_species = 0;
        nspecies = 0;

      } else {
        solve_species = 1;
        nspecies = species_names.size();

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nspecies <= species.nspecies,
            "Solids species number is higher than total species number");

        species_IDs.resize(nspecies);
        MW_sn0.resize(nspecies);

        if (solve_energy) {

          if (species.SpecificHeatModel == SPECIFICHEATMODEL::Constant) {
            specific_heat_model = SPECIFICHEATMODEL::Constant;
            cp_sn0.resize(nspecies);
          } else if (species.SpecificHeatModel == SPECIFICHEATMODEL::NASA7Polynomials) {
            specific_heat_model = SPECIFICHEATMODEL::NASA7Polynomials;
            cp_sn0.resize(nspecies*12);
          }

          H_fn0.resize(nspecies);
        }

        for (int n(0); n < nspecies; n++) {
          auto it = std::find(species.names.begin(), species.names.end(), species_names[n]);

          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != species.names.end(),
                                           "Solid species missing in input");

          const auto pos = std::distance(species.names.begin(), it);

          species_IDs[n] = species.IDs[pos];
          MW_sn0[n] = species.MW_k0[pos];

          if (solve_energy) {
            if (specific_heat_model == SPECIFICHEATMODEL::Constant) {
              cp_sn0[n] = species.cp_k0[pos];
            } else if (specific_heat_model == SPECIFICHEATMODEL::NASA7Polynomials) {
              std::copy(&species.cp_k0[pos*12], &species.cp_k0[pos*12] + 12, &cp_sn0[n*12]);
            }

            H_fn0[n] = species.H_fk0[pos];
          }
        }
      }
    } else {

      species_IDs.resize(1, 0);
      MW_sn0.resize(1);

      pp.query("molecular_weight", MW_sn0[0]);

      if (solve_energy) {
        H_fn0.resize(1);

        // Get specific heat model input ------------------------//
        std::string model;
        pp.query("specific_heat", model);

        if (amrex::toLower(model).compare("constant") == 0) {
          specific_heat_model = SPECIFICHEATMODEL::Constant;
          cp_sn0.resize(1);
          pp.get("specific_heat.constant", cp_sn0[0]);

          // Get enthalpy of formation model input ------------------------//
          pp.query("enthalpy_of_formation", H_fn0[0]);

        } else if (amrex::toLower(model).compare("nasa7-poly") == 0) {
          specific_heat_model = SPECIFICHEATMODEL::NASA7Polynomials;
          cp_sn0.resize(12);

          for (int i(0); i < 6; ++i) {

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(MW_sn0[0] > 0., "Wrong solids molecular weight");
            const amrex::Real coeff = FluidPhase::R / MW_sn0[0];

            amrex::Vector<amrex::Real> aa(0);
            std::string field = "specific_heat.NASA7.a"+std::to_string(i);
            pp.getarr(field.c_str(), aa);

            if (aa.size() == 1) {

              cp_sn0[i] = aa[0] * coeff;
              cp_sn0[i+6] = 0.;

            } else if (aa.size() == 2) {

              cp_sn0[i] = aa[0] * coeff;
              cp_sn0[i+6] = aa[1] * coeff;

            } else {

              amrex::Abort("Input error");
            }
          }

        } else {
          amrex::Abort("Don't know this specific heat model!");
        }
      }
    }

    // Flag to determine if we want to solve the solid as a mixture
    is_a_mixture = static_cast<int>(nspecies > 1);

    // Create the stoichiometric table for the fluid phase, that is associate
    // the total stoichiometric coefficient for each fluid species in each
    // reaction
    if (reactions.solve) {

      const int nreactions = reactions.nreactions;

      parameters->m_nreactions = nreactions;

      constexpr int Solid = CHEMICALPHASE::Solid;
      constexpr int Heterogeneous = REACTIONTYPE::Heterogeneous;

      // Allocate space for necessary data
      stoich_coeffs.resize(nspecies*nreactions, 0.);

      for (int n_s(0); n_s < nspecies; n_s++) {
        // Get the ID of the current species n_s
        const int species_id = species_IDs[n_s];

        // Loop over reactions to compute each contribution
        for (int q(0); q < nreactions; q++) {

          ChemicalReaction* chem_reaction = reactions.get(q);
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
                  stoich_coeffs[n_s*nreactions+q] += reactants_coeffs[pos];
                }
              }
            }

            // Add products contribution (if any)
            {
              for (int pos(0); pos < products_IDs.size(); ++pos) {
                if (species_id == products_IDs[pos] && products_chem_phases[pos] == Solid) {
                  stoich_coeffs[n_s*nreactions+q] += products_coeffs[pos];
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
        thermal_conductivity_model = SOLIDSTHERMALCONDUCTIVITYMODEL::Constant;
        kp_sn0.resize(1);
        kp_read = pp.query("thermal_conductivity.constant", kp_sn0[0]);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(kp_sn0[0] > 0, "Invalid DEM thermal conductivity.");
    } /*else {
        amrex::Abort("Unknown particle conductivity model!");
    }*/


    // FLPC and roughness are always dimension 1
    if (solve_energy) {
        flpc_read  = pp.query("flpc", flpc);
        rough_read = pp.query("min_conduction_dist", rough);
    }

    if(kp_type_read && kp_read && flpc_read && rough_read)
      do_conduction = 1;

    // Species IDs
    if (species_IDs.size() > 0) {
      d_species_IDs.resize(species_IDs.size());
      Gpu::copyAsync(Gpu::hostToDevice, species_IDs.begin(), species_IDs.end(), d_species_IDs.begin());
      parameters->m_h_species_id = species_IDs.data();
      parameters->m_d_species_id = d_species_IDs.data();
    }

    // Molecular weights
    if (MW_sn0.size() > 0) {
      d_MW_sn0.resize(MW_sn0.size());
      Gpu::copyAsync(Gpu::hostToDevice, MW_sn0.begin(), MW_sn0.end(), d_MW_sn0.begin());
      parameters->m_h_MW_sn = MW_sn0.data();
      parameters->m_d_MW_sn = d_MW_sn0.data();
    }

    // Specific heat coefficients
    if (cp_sn0.size() > 0) {
      d_cp_sn0.resize(cp_sn0.size());
      Gpu::copyAsync(Gpu::hostToDevice, cp_sn0.begin(), cp_sn0.end(), d_cp_sn0.begin());
      parameters->m_h_cp_sn = cp_sn0.data();
      parameters->m_d_cp_sn = d_cp_sn0.data();
    }

    // Enthalpies of formation
    if (H_fn0.size() > 0) {
      d_H_fn0.resize(H_fn0.size());
      Gpu::copyAsync(Gpu::hostToDevice, H_fn0.begin(), H_fn0.end(), d_H_fn0.begin());
      parameters->m_h_H_fn = H_fn0.data();
      parameters->m_d_H_fn = d_H_fn0.data();
    }

    // Stoichiometric coefficients
    if (stoich_coeffs.size() > 0) {
      d_stoich_coeffs.resize(stoich_coeffs.size());
      Gpu::copyAsync(Gpu::hostToDevice, stoich_coeffs.begin(), stoich_coeffs.end(), d_stoich_coeffs.begin());
      parameters->m_h_stoich_coeffs = stoich_coeffs.data();
      parameters->m_d_stoich_coeffs = d_stoich_coeffs.data();
    }

    //
    if (kp_sn0.size() > 0) {
      d_kp_sn0.resize(kp_sn0.size());
      Gpu::copyAsync(Gpu::hostToDevice, kp_sn0.begin(), kp_sn0.end(), d_kp_sn0.begin());
      parameters->m_h_kp_sn = kp_sn0.data();
      parameters->m_d_kp_sn = d_kp_sn0.data();
    }

    parameters->m_flpc = flpc;
    parameters->m_rough = rough;
    parameters->m_do_pfp_cond = do_conduction;

    parameters->m_specific_heat_model = specific_heat_model;
    parameters->m_thermal_conductivity_model = thermal_conductivity_model;
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
