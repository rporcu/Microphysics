#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_ParmParse.H>

#include <mfix_reactions_parms.H>
#include <mfix_ic_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_regions_parms.H>
#include <mfix_species_parms.H>

#include <string>
#include <sstream>
#include <cctype>
#include <iterator>
#include <regex>


namespace stoichiometry_aux {

std::string ltrim(const std::string& s) {
  return std::regex_replace(s, std::regex("^\\s+"), std::string(""));
}


std::string rtrim(const std::string& s) {
  return std::regex_replace(s, std::regex("\\s+$"), std::string(""));
}


std::string trim(const std::string& s) {
  return ltrim(rtrim(s));
}


std::string get_reactants(const std::string& formula)
{
  {
    std::size_t pos = formula.find("-->");
    if(pos != std::string::npos)
      return trim(formula.substr(0,pos));
  }
  {
    std::size_t pos = formula.find("<--");
    if(pos != std::string::npos)
      return trim(formula.substr(pos+3, formula.size()-(pos+3)));
  }
  {
    std::size_t pos = formula.find("<=>");
    if(pos != std::string::npos) {
      //return trim(formula.substr(0,pos)) + "+" +
      //  trim(formula.substr(pos+3, formula.size()-(pos+3)));
      amrex::Abort("Not yet implemented");
      return std::string();
    }
  }

  amrex::Abort("Error: wrong format of chemical equation");
  return formula;
}


std::string get_products(const std::string& formula)
{
  {
    std::size_t pos = formula.find("<--");
    if(pos != std::string::npos)
      return trim(formula.substr(0,pos));
  }
  {
    std::size_t pos = formula.find("-->");
    if(pos != std::string::npos)
      return trim(formula.substr(pos+3, formula.size()-(pos+3)));
  }
  {
    std::size_t pos = formula.find("<=>");
    if(pos != std::string::npos) {
      //return trim(formula.substr(0,pos)) + "+" +
      //  trim(formula.substr(pos+3, formula.size()-(pos+3)));
      amrex::Abort("Not yet implemented");
      return std::string();
    }
  }

  amrex::Abort("Error: wrong format of chemical equation");
  return formula;
}


void
get_stoichiometric_data(const std::string& s,
                        std::vector<std::string>& compounds,
                        amrex::Gpu::ManagedVector<int>& compounds_id,
                        amrex::Gpu::ManagedVector<amrex::Real>& coefficients,
                        amrex::Gpu::ManagedVector<int>& phases)
{
  std::string formula(trim(s));
  std::replace(formula.begin(), formula.end(), '+', ' ');
  
  std::istringstream iss(formula);
  
  std::vector<std::string> stoichiometry((std::istream_iterator<std::string>(iss)),
                                          std::istream_iterator<std::string>());

  // Number of compounds
  const int nc = stoichiometry.size();

  compounds_id.clear();
  compounds.clear();
  coefficients.clear();
  phases.clear();

  compounds_id.resize(nc);
  compounds.resize(nc);
  coefficients.resize(nc);
  phases.resize(nc);

  for(int n(0); n < nc; n++)
  {
    std::string single_compound = stoichiometry[n];

    std::size_t pos(0);
    while (std::isdigit(single_compound.at(pos)) or single_compound.at(pos) == '.')
      pos++;

    if(pos == 0)
      coefficients[n] = 1.;
    else {
      coefficients[n] = std::stod(single_compound.substr(0, pos));
      single_compound.erase(0, pos);
    }

    pos = std::min(single_compound.find("(g)"), single_compound.find("(s)"));
    if (pos != std::string::npos) {
      std::string loc_phase = single_compound.substr(pos, pos+3);
      if (loc_phase.compare("(g)") == 0)
        phases[n] = CHEMICALPHASE::Fluid;
      else if (loc_phase.compare("(l)") == 0)
        phases[n] = CHEMICALPHASE::Fluid;
      else if (loc_phase.compare("(s)") == 0)
        phases[n] = CHEMICALPHASE::Solid;
      else
        amrex::Abort("Error: not recognized phase in stoichiometric equation");

      single_compound.erase(pos, pos+3);
    }
    else {
      amrex::Abort("Error: wrong format for chemical equation");
    }
    
    compounds[n] = trim(single_compound);

    auto it = std::find(SPECIES::species.begin(), SPECIES::species.end(),
        compounds[n]);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != SPECIES::species.end(),
        "Error: species " + compounds[n] + " is present in a reaction but " +
        "it is not present in the list of declared species");

    const auto it_pos = std::distance(SPECIES::species.begin(), it);

    compounds_id[n] = SPECIES::species_id[it_pos];
  }
}


void
parse_chemical_reaction(const std::string& equation,
                        amrex::Gpu::ManagedVector<int>& phases,
                        int& reaction_type,
                        std::vector<std::string>& reactants,
                        amrex::Gpu::ManagedVector<int>& reactants_id,
                        amrex::Gpu::ManagedVector<amrex::Real>& reactants_coeffs,
                        amrex::Gpu::ManagedVector<int>& reactants_phases,
                        std::vector<std::string>& products,
                        amrex::Gpu::ManagedVector<int>& products_id,
                        amrex::Gpu::ManagedVector<amrex::Real>& products_coeffs,
                        amrex::Gpu::ManagedVector<int>& products_phases)
{
  // Clear reactants containers
  reactants.clear();
  reactants_coeffs.clear();
  reactants_phases.clear();

  // Clear products containers
  products.clear();
  products_coeffs.clear();
  products_phases.clear();

  // Get the reaction part of the equation
  std::string reaction_part = get_reactants(equation);
  
  // Get the products part of the equation
  std::string production_part = get_products(equation);

  // Get the reaction part stoichiometric coefficoents, elements and phases
  get_stoichiometric_data(reaction_part, reactants, reactants_id,
      reactants_coeffs, reactants_phases);

  // Multiply reaction coefficients by -1
  for (auto it = reactants_coeffs.begin(); it != reactants_coeffs.end(); it++)
    *it *= -1;

  // Get the production part stoichiometric coefficoents, elements and phases
  get_stoichiometric_data(production_part, products, products_id,
      products_coeffs, products_phases);

  std::for_each(reactants_phases.begin(), reactants_phases.end(),
    [&phases] (int phase)
  {
    if (std::find(phases.begin(), phases.end(), phase) == phases.end())
      phases.push_back(phase);
  });

  std::for_each(products_phases.begin(), products_phases.end(),
    [&phases] (int phase)
  {
    if (std::find(phases.begin(), phases.end(), phase) == phases.end())
      phases.push_back(phase);
  });

  if (phases.size() == 1)
    reaction_type = REACTIONTYPE::Homogeneous;
  else if (phases.size() >= 2)
    reaction_type = REACTIONTYPE::Heterogeneous;
  else
    amrex::Abort("Error: unrecognized reaction type");

  return;
}

} // end namespace stoichiometry_aux


using namespace stoichiometry_aux;


namespace REACTIONS
{
  // Switch for turning on/off chemical reactions modeling
  bool solve = false;

  // Names of chemical reactions allowed by the model
  std::vector<std::string> reactions;

  // Number of chemical reactions allowed by the model
  int nreactions = 0;

  // Vector of allowed chemical reactions
  amrex::Vector<ChemicalReaction> chemical_reactions;

  // Initialization: read input parameters and set up reactions
  void Initialize ()
  {
    amrex::ParmParse pp("chemistry");

    if (pp.contains("solve"))
    {
      // Get the list of reactions names to define chem equations
      pp.getarr("solve", reactions);

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(reactions.size() > 0,
          "No input provided for chemistry.solve");

      // Disable the species solver if the species are defined as "None" (case
      // insensitive) or 0
      if (amrex::toLower(reactions[0]).compare("none") == 0 or
          (reactions[0]).compare("0") == 0) {
        solve = false;
        reactions.clear();
        nreactions = 0;
      }
      else {
        solve = true;
        nreactions = reactions.size();
      }

      if (solve) {
        for (int n(0); n < nreactions; n++) {
          // Declare the string where we store the equation
          std::string reaction_equation;

          // Get the reation equation relative to given reaction name
          pp.get((reactions[n]+".reaction").c_str(), reaction_equation);

          ChemicalReaction reaction(reaction_equation);

          // Add new chemical reaction to the back of vector
          chemical_reactions.push_back(reaction);
        }
      }
    }
  }

  // REACTION_T Class Constructor
  ChemicalReaction::ChemicalReaction (const std::string& reaction)
    : m_reaction(reaction)
    , m_phases(0)
    , m_reaction_type(REACTIONTYPE::Invalid)
    , m_reactants(0)
    , m_reactants_id(0)
    , m_reactants_coeffs(0)
    , m_reactants_phases(0)
    , m_products(0)
    , m_products_id(0)
    , m_products_coeffs(0)
    , m_products_phases(0)
  {
    parse_chemical_reaction(m_reaction, m_phases, m_reaction_type, m_reactants,
        m_reactants_id, m_reactants_coeffs, m_reactants_phases, m_products,
        m_products_id, m_products_coeffs, m_products_phases);
  }
}
