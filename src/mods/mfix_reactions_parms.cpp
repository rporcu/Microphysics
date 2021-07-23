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
#include <mfix_reactions_rates_K.H>

#include <string>
#include <sstream>
#include <cctype>
#include <iterator>
#include <regex>


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
{}


Reactions::~Reactions()
{
  if (solve)
    for (int n(0); n < nreactions; n++)
      delete m_chemical_reactions[n];
}


// Initialization: read input parameters and set up reactions
void Reactions::Initialize () {

  amrex::ParmParse pp("chemistry");

  if (pp.contains("solve")) {

    // Get the list of reactions names to define chem equations
    pp.getarr("solve", reactions);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(reactions.size() > 0,
        "No input provided for chemistry.solve");

    // Disable the species solver if the species are defined as "None" (case
    // insensitive) or 0
    if (amrex::toLower(reactions[0]).compare("none") == 0 || (reactions[0]).compare("0") == 0) {
      solve = false;
      reactions.clear();
      nreactions = 0;
      reaction_equations.clear();
    } else {
      solve = true;
      nreactions = reactions.size();
      reaction_equations.clear();
      reaction_equations.resize(nreactions);
    }

    if (solve) {
      for (int n(0); n < nreactions; n++) {
        // Get the reation equation relative to given reaction name
        pp.get((reactions[n]+".reaction").c_str(), reaction_equations[n]);
      }
    }
  }

  if (solve) {

    m_chemical_reactions.resize(nreactions, nullptr);
    
    for (int n(0); n < nreactions; n++) {

      const std::string& equation = reaction_equations[n];

      m_chemical_reactions[n] = new ChemicalReaction(equation);
    }
  }

  // Comment
  Gpu::HostVector<int> h_types(nreactions);
  std::transform(m_chemical_reactions.begin(), m_chemical_reactions.end(),
                 h_types.begin(), [](auto reaction){ return reaction->m_reaction_type; });
  d_types.resize(nreactions);
  Gpu::copyAsync(Gpu::hostToDevice, h_types.begin(), h_types.end(), d_types.begin());
  //const int* p_types = d_types.data();

  // Comment
  Gpu::HostVector<int> h_nphases(nreactions);
  std::transform(m_chemical_reactions.begin(), m_chemical_reactions.end(),
                 h_nphases.begin(), [](auto reaction){ return reaction->m_phases.size(); });
  d_nphases.resize(nreactions);
  Gpu::copyAsync(Gpu::hostToDevice, h_nphases.begin(), h_nphases.end(), d_nphases.begin());
  //const int* p_nphases = d_nphases.data();

  // Comment
  Gpu::HostVector<const int*> h_phases(nreactions);
  std::transform(m_chemical_reactions.begin(), m_chemical_reactions.end(),
                 h_phases.begin(), [](auto reaction){ return reaction->m_phases.data(); });
  d_phases.resize(nreactions);
  Gpu::copyAsync(Gpu::hostToDevice, h_phases.begin(), h_phases.end(), d_phases.begin());
  //const int** p_phases = d_phases.data();

  // Comment
  Gpu::HostVector<int> h_nreactants(nreactions);
  std::transform(m_chemical_reactions.begin(), m_chemical_reactions.end(),
                 h_nreactants.begin(), [](auto reaction){ return reaction->m_reactants.size(); });
  d_nreactants.resize(nreactions);
  Gpu::copyAsync(Gpu::hostToDevice, h_nreactants.begin(), h_nreactants.end(), d_nreactants.begin());
  const int* p_nreactants = d_nreactants.data();

  // Comment
  Gpu::HostVector<const int*> h_reactants_id(nreactions);
  std::transform(m_chemical_reactions.begin(), m_chemical_reactions.end(),
                 h_reactants_id.begin(), [](auto reaction){ return reaction->m_reactants_id.data(); });
  d_reactants_id.resize(nreactions);
  Gpu::copyAsync(Gpu::hostToDevice, h_reactants_id.begin(), h_reactants_id.end(), d_reactants_id.begin());
  const int** p_reactants_id = d_reactants_id.data();

  // Comment
  Gpu::HostVector<const Real*> h_reactants_coeffs(nreactions);
  std::transform(m_chemical_reactions.begin(), m_chemical_reactions.end(),
                 h_reactants_coeffs.begin(), [](auto reaction){ return reaction->m_reactants_coeffs.data(); });
  d_reactants_coeffs.resize(nreactions);
  Gpu::copyAsync(Gpu::hostToDevice, h_reactants_coeffs.begin(), h_reactants_coeffs.end(), d_reactants_coeffs.begin());
  const Real** p_reactants_coeffs = d_reactants_coeffs.data();

  // Comment
  Gpu::HostVector<const int*> h_reactants_phases(nreactions);
  std::transform(m_chemical_reactions.begin(), m_chemical_reactions.end(),
                 h_reactants_phases.begin(), [](auto reaction){ return reaction->m_reactants_phases.data(); });
  d_reactants_phases.resize(nreactions);
  Gpu::copyAsync(Gpu::hostToDevice, h_reactants_phases.begin(), h_reactants_phases.end(), d_reactants_phases.begin());
  const int** p_reactants_phases = d_reactants_phases.data();

  // Comment
  Gpu::HostVector<int> h_nproducts(nreactions);
  std::transform(m_chemical_reactions.begin(), m_chemical_reactions.end(),
                 h_nproducts.begin(), [](auto reaction){ return reaction->m_products.size(); });
  d_nproducts.resize(nreactions);
  Gpu::copyAsync(Gpu::hostToDevice, h_nproducts.begin(), h_nproducts.end(), d_nproducts.begin());
  const int* p_nproducts = d_nproducts.data();

  // Comment
  Gpu::HostVector<const int*> h_products_id(nreactions);
  std::transform(m_chemical_reactions.begin(), m_chemical_reactions.end(),
                 h_products_id.begin(), [](auto reaction){ return reaction->m_products_id.data(); });
  d_products_id.resize(nreactions);
  Gpu::copyAsync(Gpu::hostToDevice, h_products_id.begin(), h_products_id.end(), d_products_id.begin());
  const int** p_products_id = d_products_id.data();

  // Comment
  Gpu::HostVector<const Real*> h_products_coeffs(nreactions);
  std::transform(m_chemical_reactions.begin(), m_chemical_reactions.end(),
                 h_products_coeffs.begin(), [](auto reaction){ return reaction->m_products_coeffs.data(); });
  d_products_coeffs.resize(nreactions);
  Gpu::copyAsync(Gpu::hostToDevice, h_products_coeffs.begin(), h_products_coeffs.end(), d_products_coeffs.begin());
  const Real** p_products_coeffs = d_products_coeffs.data();

  // Comment
  Gpu::HostVector<const int*> h_products_phases(nreactions);
  std::transform(m_chemical_reactions.begin(), m_chemical_reactions.end(),
                 h_products_phases.begin(), [](auto reaction){ return reaction->m_products_phases.data(); });
  d_products_phases.resize(nreactions);
  Gpu::copyAsync(Gpu::hostToDevice, h_products_phases.begin(), h_products_phases.end(), d_products_phases.begin());
  const int** p_products_phases = d_products_phases.data();

  parameters = new ReactionsParms(nreactions, p_nreactants, p_reactants_id,
                                  p_reactants_coeffs, p_reactants_phases,
                                  p_nproducts, p_products_id, p_products_coeffs,
                                  p_products_phases);

  Gpu::synchronize();
}


// REACTION_T Class Constructor
ChemicalReaction::ChemicalReaction (const std::string& reaction)
  : m_reaction_type(REACTIONTYPE::Invalid)
  , m_phases(0)
  , m_reaction(reaction)
  , m_reactants(0)
  , m_reactants_id(0)
  , m_reactants_coeffs(0)
  , m_reactants_phases(0)
  , m_products(0)
  , m_products_id(0)
  , m_products_coeffs(0)
  , m_products_phases(0)
{
  Gpu::HostVector<int> h_phases;

  Gpu::HostVector<int> h_reactants_id;
  Gpu::HostVector<Real> h_reactants_coeffs;
  Gpu::HostVector<int> h_reactants_phases;

  Gpu::HostVector<int> h_products_id;
  Gpu::HostVector<Real> h_products_coeffs;
  Gpu::HostVector<int> h_products_phases;

  parse_reaction(m_reaction, h_phases, m_reaction_type, m_reactants,
                 h_reactants_id, h_reactants_coeffs, h_reactants_phases, m_products,
                 h_products_id, h_products_coeffs, h_products_phases);

  m_phases.resize(h_phases.size());
  Gpu::copyAsync(Gpu::hostToDevice, h_phases.begin(), h_phases.end(), m_phases.begin());

  m_reactants_id.resize(h_reactants_id.size());
  Gpu::copyAsync(Gpu::hostToDevice, h_reactants_id.begin(), h_reactants_id.end(), m_reactants_id.begin());

  m_reactants_coeffs.resize(h_reactants_coeffs.size());
  Gpu::copyAsync(Gpu::hostToDevice, h_reactants_coeffs.begin(), h_reactants_coeffs.end(), m_reactants_coeffs.begin());

  m_reactants_phases.resize(h_reactants_phases.size());
  Gpu::copyAsync(Gpu::hostToDevice, h_reactants_phases.begin(), h_reactants_phases.end(), m_reactants_phases.begin());

  m_products_id.resize(h_products_id.size());
  Gpu::copyAsync(Gpu::hostToDevice, h_products_id.begin(), h_products_id.end(), m_products_id.begin());

  m_products_coeffs.resize(h_products_coeffs.size());
  Gpu::copyAsync(Gpu::hostToDevice, h_products_coeffs.begin(), h_products_coeffs.end(), m_products_coeffs.begin());

  m_products_phases.resize(h_products_phases.size());
  Gpu::copyAsync(Gpu::hostToDevice, h_products_phases.begin(), h_products_phases.end(), m_products_phases.begin());

  Gpu::synchronize();
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
                                          std::vector<std::string>& compounds,
                                          amrex::Gpu::HostVector<int>& compounds_id,
                                          amrex::Gpu::HostVector<amrex::Real>& coefficients,
                                          amrex::Gpu::HostVector<int>& phases)
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

  compounds_id.resize(nc);
  compounds.resize(nc);
  coefficients.resize(nc);
  phases.resize(nc);

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
ChemicalReaction::parse_reaction(const std::string& equation,
                                 amrex::Gpu::HostVector<int>& phases,
                                 int& reaction_type,
                                 std::vector<std::string>& reactants,
                                 amrex::Gpu::HostVector<int>& reactants_id,
                                 amrex::Gpu::HostVector<amrex::Real>& reactants_coeffs,
                                 amrex::Gpu::HostVector<int>& reactants_phases,
                                 std::vector<std::string>& products,
                                 amrex::Gpu::HostVector<int>& products_id,
                                 amrex::Gpu::HostVector<amrex::Real>& products_coeffs,
                                 amrex::Gpu::HostVector<int>& products_phases)
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
