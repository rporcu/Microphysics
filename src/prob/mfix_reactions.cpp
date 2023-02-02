#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_ParmParse.H>

#include <mfix.H>
#include <mfix_reactions.H>
#include <mfix_ic.H>
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_regions.H>
#include <mfix_species.H>
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


MFIXReactions::MFIXReactions()
  : m_solve(0)
  , m_nreactions(0)
  , m_reactions(0)
  , m_reaction_equations(0)
  , m_is_initialized(0)
  , m_parameters(nullptr)
  , m_chemical_reactions(0)
{}


MFIXReactions::~MFIXReactions()
{
  if (m_solve)
    for (int n(0); n < m_nreactions; n++)
      delete m_chemical_reactions[n];
}


// Initialization: read input parameters and set up reactions
void MFIXReactions::Initialize (const MFIXSpecies& species)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species.isInitialized(),
      "Species not initialized. Can't initialize reactions before species initialization");

  // Set the initialization flag
  m_is_initialized = 1;

  amrex::ParmParse pp("chemistry");

  if (pp.contains("solve")) {

    // Get the list of reactions names to define chem equations
    pp.getarr("solve", m_reactions);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_reactions.size() > 0,
        "No input provided for chemistry.solve");

    // Disable the species solver if the species are defined as "None" (case
    // insensitive) or 0
    if (amrex::toLower(m_reactions[0]).compare("none") == 0 || (m_reactions[0]).compare("0") == 0) {
      m_solve = 0;
      m_reactions.clear();
      m_nreactions = 0;
      m_reaction_equations.clear();
    } else {
      m_solve = 1;
      m_nreactions = m_reactions.size();
      m_reaction_equations.clear();
      m_reaction_equations.resize(m_nreactions);
    }

    if (m_solve) {
      for (int n(0); n < m_nreactions; n++) {
        // Get the reation equation relative to given reaction name
        Vector<std::string> equation(0);
        pp.getarr((m_reactions[n]+".reaction").c_str(), equation);
        for (auto elem: equation) m_reaction_equations[n] += elem;
      }
    }
  }

  if (m_solve) {

    m_chemical_reactions.resize(m_nreactions, nullptr);
    
    for (int n(0); n < m_nreactions; n++) {

      const std::string& equation = m_reaction_equations[n];

      m_chemical_reactions[n] = new MFIXChemicalReaction(equation, species);
    }
  }

  m_parameters = new MFIXReactionsParms(m_nreactions);
}


// REACTION_T Class Constructor
MFIXChemicalReaction::MFIXChemicalReaction (const std::string& reaction,
                                            const MFIXSpecies& species)
  : m_type(ReactionType::Invalid)
  , m_formula(reaction)
  , m_phases(0)
  , m_reactants(0)
  , m_reactants_IDs(0)
  , m_reactants_coeffs(0)
  , m_reactants_phases(0)
  , m_products(0)
  , m_products_IDs(0)
  , m_products_coeffs(0)
  , m_products_phases(0)
  , m_mass_balance_tolerance(1.e-12)
{
  ParmParse pp_reactions("chemistry");
  pp_reactions.query("mass_balance_tolerance", m_mass_balance_tolerance);

  parse_reaction(species);
}


std::string 
MFIXChemicalReaction::parse_reactants(const std::string& formula)
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
MFIXChemicalReaction::parse_products(const std::string& formula)
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
MFIXChemicalReaction::parse_stoichiometric_data(const std::string& s,
                                                amrex::Vector<std::string>& compounds,
                                                amrex::Vector<int>& compounds_id,
                                                amrex::Vector<amrex::Real>& coefficients,
                                                amrex::Vector<int>& phases,
                                                const MFIXSpecies& species)
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
        phases[n] = ChemicalPhase::Fluid;
      else if (loc_phase.compare("(l)") == 0)
        phases[n] = ChemicalPhase::Fluid;
      else if (loc_phase.compare("(s)") == 0)
        phases[n] = ChemicalPhase::Solid;
      else
        amrex::Abort("Error: not recognized phase in stoichiometric equation");

      single_compound.erase(pos, pos+3);
    }
    else {
      amrex::Abort("Error: wrong format for chemical equation");
    }
    
    compounds[n] = chemistry_aux::trim(single_compound);

    const auto& names = species.names();
    auto it = std::find(names.begin(), names.end(), compounds[n]);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != names.end(),
        "Error: species " + compounds[n] + " is present in a reaction but " +
        "it is not present in the list of declared species");

    const auto it_pos = std::distance(names.begin(), it);

    compounds_id[n] = species.IDs(it_pos);
  }
}


void
MFIXChemicalReaction::check_mass_balance (const MFIXSpecies& species)
{
  AMREX_ASSERT(m_reactants_IDs.size() == m_reactants_coeffs.size());
  AMREX_ASSERT(m_products_IDs.size() == m_products_coeffs.size());

  Real reactants_mass(0.);
  Real products_mass(0.);

  for (int n(0); n < m_reactants_IDs.size(); ++n) {
    int ID = m_reactants_IDs[n];
    reactants_mass += (-1.*m_reactants_coeffs[n])*species.MW_k(ID);
  }

  for (int n(0); n < m_products_IDs.size(); ++n) {
    int ID = m_products_IDs[n];
    products_mass += m_products_coeffs[n]*species.MW_k(ID);
  }

  Real diff_mass = std::abs(reactants_mass-products_mass);

  if(diff_mass > m_mass_balance_tolerance) {

    Print() << "\nUnbalanced reaction: " << m_formula << "\n\n";

    Print() << "Reactants mass balance: " << reactants_mass << "\n"
            << "Products mass balance: " << products_mass << "\n\n";

    Print() << "Mass balance difference: " << diff_mass << "\n"
            << "Mass balance tolerance: " << m_mass_balance_tolerance << "\n\n";

    amrex::Abort("Fix inputs either in species mass fractions or chemical reactions");
  }

  return;
}


void
MFIXChemicalReaction::parse_reaction(const MFIXSpecies& species)
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
  std::string reaction_part = parse_reactants(formula);
  
  // Get the products part of the formula
  std::string production_part = parse_products(formula);

  // Get the reaction part stoichiometric coefficients, elements and phases
  parse_stoichiometric_data(reaction_part, m_reactants, m_reactants_IDs,
                            m_reactants_coeffs, m_reactants_phases, species);

  // Multiply reaction coefficients by -1
  for (int i(0); i < m_reactants_coeffs.size(); ++i)
    m_reactants_coeffs[i] *= -1;

  // Get the production part stoichiometric coefficients, elements and phases
  parse_stoichiometric_data(production_part, m_products, m_products_IDs,
                            m_products_coeffs, m_products_phases, species);

  check_mass_balance(species);

  // Fill m_phases with all the phases found in reactants
  for (const int phase: m_reactants_phases)
    if (std::find(m_phases.begin(), m_phases.end(), phase) == m_phases.end())
      m_phases.push_back(phase);

  // Add to m_phases all the phases found in products
  for (const int phase: m_products_phases)
    if (std::find(m_phases.begin(), m_phases.end(), phase) == m_phases.end())
      m_phases.push_back(phase);

  // Set the type of reaction, whether homogeneous or heterogeneous
  if (m_phases.size() == 1)
    m_type = ReactionType::Homogeneous;
  else if (m_phases.size() >= 2)
    m_type = ReactionType::Heterogeneous;
  else
    amrex::Abort("Error: unrecognized reaction type");

  return;
}
