#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_ParmParse.H>

#include <mfix.H>
#include <mfix_reactions_parms.H>
#include <mfix_ic_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_regions_parms.H>
#include <mfix_species_parms.H>
#include <mfix_des_rrates_K.H>

#include <string>
#include <sstream>
#include <cctype>
#include <iterator>
#include <regex>
#include <algorithm>


namespace chemistry_aux {

std::string ltrim(const std::string& s) {
  return std::regex_replace(s, std::regex("^\\s+"), std::string(""));
}


std::string rtrim(const std::string& s) {
  return std::regex_replace(s, std::regex("\\s+$"), std::string(""));
}


std::string trim(const std::string& s) {
  return ltrim(rtrim(s));
}

} // end namespace chemistry_aux


Reactions::Reactions()
  : solve(0)
  , nreactions(0)
  , reactions(0)
  , reaction_equations(0)
  , m_chemical_reactions(0)
  , heterogeneous(0)
  , is_initialized(0)
{}


Reactions::~Reactions()
{
  if (solve)
    for (int n(0); n < nreactions; n++)
      delete m_chemical_reactions[n];
}


// Initialization: read input parameters and set up reactions
void Reactions::Initialize (const Species& species)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species.is_initialized,
      "Species not initialized. Can't initialize reactions before species initialization");

  // Set the initialization flag
  is_initialized = 1;

  amrex::ParmParse pp("chemistry");

  if (pp.contains("solve")) {

    // Get the list of reactions names to define chem equations
    pp.getarr("solve", reactions);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(reactions.size() > 0,
        "No input provided for chemistry.solve");

    // Disable the species solver if the species are defined as "None" (case
    // insensitive) or 0
    if (amrex::toLower(reactions[0]).compare("none") == 0 || (reactions[0]).compare("0") == 0) {
      solve = 0;
      reactions.clear();
      nreactions = 0;
      reaction_equations.clear();
    } else {
      solve = 1;
      nreactions = reactions.size();
      reaction_equations.clear();
      reaction_equations.resize(nreactions);
    }

    if (solve) {
      for (int n(0); n < nreactions; n++) {
        // Get the reation equation relative to given reaction name
        Vector<std::string> equation(0);
        pp.getarr((reactions[n]+".reaction").c_str(), equation);
        for (auto elem: equation) reaction_equations[n] += elem;
      }
    }
  }

  if (solve) {

    m_chemical_reactions.resize(nreactions, nullptr);
    
    for (int n(0); n < nreactions; n++) {

      const std::string& equation = reaction_equations[n];

      m_chemical_reactions[n] = new ChemicalReaction(equation, species);
    }
  }

  parameters = new ReactionsParms(nreactions);
}


// REACTION_T Class Constructor
ChemicalReaction::ChemicalReaction (const std::string& reaction,
                                    const Species& species)
  : m_type(REACTIONTYPE::Invalid)
  , m_phases(0)
  , m_formula(reaction)
  , m_reactants(0)
  , m_reactants_id(0)
  , m_reactants_coeffs(0)
  , m_reactants_phases(0)
  , m_products(0)
  , m_products_id(0)
  , m_products_coeffs(0)
  , m_products_phases(0)
{
  parse_reaction(species);
}


std::string 
ChemicalReaction::get_reactants(const std::string& formula)
{
  {
    std::size_t pos = formula.find("-->");
    if(pos != std::string::npos)
      return chemistry_aux::trim(formula.substr(0,pos));
  }
  {
    std::size_t pos = formula.find("<--");
    if(pos != std::string::npos)
      return chemistry_aux::trim(formula.substr(pos+3, formula.size()-(pos+3)));
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


std::string 
ChemicalReaction::get_products(const std::string& formula)
{
  {
    std::size_t pos = formula.find("<--");
    if(pos != std::string::npos)
      return chemistry_aux::trim(formula.substr(0,pos));
  }
  {
    std::size_t pos = formula.find("-->");
    if(pos != std::string::npos)
      return chemistry_aux::trim(formula.substr(pos+3, formula.size()-(pos+3)));
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
ChemicalReaction::get_stoichiometric_data(const std::string& s,
                                          amrex::Vector<std::string>& compounds,
                                          amrex::Vector<int>& compounds_id,
                                          amrex::Vector<amrex::Real>& coefficients,
                                          amrex::Vector<int>& phases,
                                          const Species& species)
{
  std::string formula(chemistry_aux::trim(s));
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

  compounds_id.resize(nc, -1);
  compounds.resize(nc, std::string());
  coefficients.resize(nc, 0.);
  phases.resize(nc, -1);

  for(int n(0); n < nc; n++)
  {
    std::string single_compound = stoichiometry[n];

    std::size_t pos(0);
    while (std::isdigit(single_compound.at(pos)) || single_compound.at(pos) == '.')
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
    
    compounds[n] = chemistry_aux::trim(single_compound);

    auto it = std::find(species.names.begin(), species.names.end(),
        compounds[n]);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != species.names.end(),
        "Error: species " + compounds[n] + " is present in a reaction but " +
        "it is not present in the list of declared species");

    const auto it_pos = std::distance(species.names.begin(), it);

    compounds_id[n] = species.IDs[it_pos];
  }
}


void
ChemicalReaction::parse_reaction(const Species& species)
{
  // remove spaces from reaction reaction_formula string
  std::string formula(m_formula);
  formula.erase(std::remove(formula.begin(), formula.end(), ' '), formula.end());

  // Clear reactants containers
  m_reactants.clear();
  m_reactants_coeffs.clear();
  m_reactants_phases.clear();

  // Clear products containers
  m_products.clear();
  m_products_coeffs.clear();
  m_products_phases.clear();

  // Get the reaction part of the formula
  std::string reaction_part = get_reactants(formula);
  
  // Get the products part of the formula
  std::string production_part = get_products(formula);

  // Get the reaction part stoichiometric coefficients, elements and phases
  get_stoichiometric_data(reaction_part, m_reactants, m_reactants_id,
      m_reactants_coeffs, m_reactants_phases, species);

  // Multiply reaction coefficients by -1
  for (size_t i(0); i < m_reactants_coeffs.size(); ++i)
    m_reactants_coeffs[i] *= -1;

  // Get the production part stoichiometric coefficients, elements and phases
  get_stoichiometric_data(production_part, m_products, m_products_id,
      m_products_coeffs, m_products_phases, species);

  // Fill m_phases with all the phases found in reactants
  for (const int phase: m_reactants_phases)
    if (std::find(m_phases.begin(), m_phases.end(), phase) == m_phases.end())
      m_phases.push_back(phase);

  // Add to m_phases all the phases found in products
  for (const int phase: m_products_phases)
    if (std::find(m_phases.begin(), m_phases.end(), phase) == m_phases.end())
      m_phases.push_back(phase);

  // Set the type of reaction, wether homogeneous or heterogeneous
  if (m_phases.size() == 1)
    m_type = REACTIONTYPE::Homogeneous;
  else if (m_phases.size() >= 2)
    m_type = REACTIONTYPE::Heterogeneous;
  else
    amrex::Abort("Error: unrecognized reaction type");

  return;
}
