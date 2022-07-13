#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>

#include <mfix_species.H>
#include <mfix_fluid.H>


using namespace amrex;


MFIXSpecies::MFIXSpecies()
  : m_diffusivity_model(DiffusivityModel::Invalid)
  , m_specific_heat_model(MFIXSpecies::SpecificHeatModel::Invalid)
  , m_solve(0)
  , m_nspecies(0)
  , m_names(0)
  , m_IDs(0)
  , m_MW_k(0)
  , m_D(0)
  , m_cp_k(0)
  , m_H_fk(0)
  , m_is_initialized(0)
{}


void
MFIXSpecies::Initialize ()
{
  // Set the initialization flag
  m_is_initialized = 1;

  int solve_enthalpy(0);

  amrex::ParmParse ppMFIX("mfix");
  ppMFIX.query("advect_enthalpy", solve_enthalpy);

  amrex::ParmParse pp("species");

  if (pp.contains("solve")) {

    pp.getarr("solve", m_names);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_names.size() > 0,
                                     "No input provided for species.solve");

    // Disable the species solver if the species are defined as "None" (case
    // insensitive) or 0
    if (amrex::toLower(m_names[0]).compare("none") == 0) {
      m_solve = 0;

    } else {
      m_solve = 1;
      m_nspecies = m_names.size();

      m_IDs.resize(m_nspecies);
      for (int n(0); n < m_nspecies; n++) {
        m_IDs[n] = n;
      }

      m_MW_k.resize(m_nspecies);

      // Get species temperature inputs -----------------------------//
      if (solve_enthalpy) {

        std::string model;
        pp.query("specific_heat", model);

        if (amrex::toLower(model).compare("constant") == 0) {
          m_specific_heat_model = MFIXSpecies::SpecificHeatModel::Constant;
          m_cp_k.resize(m_nspecies);
        } else if (amrex::toLower(model).compare("nasa7-poly") == 0) {
          m_specific_heat_model = MFIXSpecies::SpecificHeatModel::NASA7Polynomials;
          m_cp_k.resize(m_nspecies*12);
        } else {
          amrex::Abort("Don't know this specific heat model!");
        }

        m_H_fk.resize(m_nspecies);
      }
    }

    if (m_solve) {
      // Get molecular weights input --------------------------------//
      for (int n(0); n < m_nspecies; n++) {
        std::string name = "species." + m_names[n];
        amrex::ParmParse ppSpecies(name.c_str());

        if (!ppSpecies.query("molecular_weight", m_MW_k[n])) {

          if (amrex::ParallelDescriptor::IOProcessor()) {
            std::string message = "Input not provided. Assuming MW_" + m_names[n] + " = 0";
            amrex::Warning(message.c_str());
          }
        }
      }

      // Get diffusivity model input --------------------------------//
      {
        std::string model;
        pp.get("diffusivity", model);

        if (amrex::toLower(model).compare("constant") == 0) {
          m_diffusivity_model = DiffusivityModel::Constant;

          pp.get("diffusivity.constant", m_D);

        } else {
          amrex::Abort("Unknown species mass diffusivity model!");
        }
      }

      if (solve_enthalpy) {
        // Get specific heat model input ------------------------//
        if (m_specific_heat_model == MFIXSpecies::SpecificHeatModel::Constant) {

          for (int n(0); n < m_nspecies; n++) {
            std::string name = "species." + m_names[n];
            amrex::ParmParse ppSpecies(name.c_str());
            ppSpecies.get("specific_heat.constant", m_cp_k[n]);

            if(!ppSpecies.query("enthalpy_of_formation", m_H_fk[n])) {

              if (amrex::ParallelDescriptor::IOProcessor()) {
                std::string message = "Input not provided. Assuming Hf_" + m_names[n] + " = 0";
                amrex::Warning(message.c_str());
              }
            }
          }

        } else if (m_specific_heat_model == MFIXSpecies::SpecificHeatModel::NASA7Polynomials) {

          for (int n(0); n < m_nspecies; n++) {
            // Non-Normalization coefficient
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_MW_k[n] > 0., "Wrong molecular weight inputs");
            const amrex::Real coeff = MFIXFluidPhase::R / m_MW_k[n];

            std::string name = "species." + m_names[n];
            amrex::ParmParse ppSpecies(name.c_str());

            for (int i(0); i < 6; ++i) {
              amrex::Vector<amrex::Real> aa(0);
              std::string field = "specific_heat.NASA7.a"+std::to_string(i);
              ppSpecies.getarr(field.c_str(), aa);

              if (aa.size() == 1) {

                m_cp_k[n*12 + i] = aa[0] * coeff;
                m_cp_k[n*12 + i+6] = 0.;

              } else if (aa.size() == 2) {

                m_cp_k[n*12 + i] = aa[0] * coeff;
                m_cp_k[n*12 + i+6] = aa[1] * coeff;

              } else {

                amrex::Abort("Input error");
              }
            }

          }

        } else {
          amrex::Abort("Unknown specific heat model");

          // Get enthalpy of formation model input ------------------------//
          for (int n(0); n < m_nspecies; n++) {
            std::string name = "species." + m_names[n];
            amrex::ParmParse ppSpecies(name.c_str());

            if(!ppSpecies.query("enthalpy_of_formation", m_H_fk[n])) {

              if (amrex::ParallelDescriptor::IOProcessor()) {
                std::string message = "Input not provided. Assuming Hf_" + m_names[n] + " = 0";
                amrex::Warning(message.c_str());
              }
            }
          }
        }

      }
    }
  }
}
