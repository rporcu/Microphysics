#ifndef MFIX_REACTIONS_PARMS_H_
#define MFIX_REACTIONS_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>

#include <AMReX_Vector.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Gpu.H>


struct REACTIONTYPE { enum {Invalid=-1, Homogeneous, Heterogeneous}; };

struct CHEMICALPHASE { enum {Invalid=-1, Fluid, Solid}; };


// Forward declaration
class ChemicalReaction;


class ReactionsParms {

  public:
    const int m_nreactions;

    AMREX_GPU_HOST_DEVICE
    ReactionsParms ()
      : m_nreactions(0)
    {}

    AMREX_GPU_HOST_DEVICE
    ReactionsParms (const int nreactions)
      : m_nreactions(nreactions)
    {}
};


class Reactions {

  public:
    // Constructor
    Reactions();

    // Destructor
    ~Reactions();

    // Max number of chemical reactions
    static constexpr int NMAX = 10;

    // Switch for turning on/off chemical reactions modeling
    int solve;

    // Number of chemical reactions allowed by the model
    int nreactions;

    // Names of chemical reactions allowed by the model
    std::vector<std::string> reactions;

    // Chemical reactions equations
    std::vector<std::string> reaction_equations;

    int is_initialized;

    // Initialization procedure
    void Initialize(const Species& species);

    // Getter
    ChemicalReaction* get(const int q) const
    {
      AMREX_ALWAYS_ASSERT(q < nreactions);
      return m_chemical_reactions[q];
    }

    ReactionsParms* parameters;

  private:
    // Vector of chemical reactions
    amrex::Vector<ChemicalReaction*> m_chemical_reactions;
};


// Class for chemical reaction details
class ChemicalReaction {

  public:
    // Default constructor
    explicit ChemicalReaction () = default;

    // class constructor
    explicit ChemicalReaction (const std::string& reaction,
                               const Species& species);

    ChemicalReaction (const ChemicalReaction&) = delete;

    ChemicalReaction& operator = (const ChemicalReaction&) = delete;

    // Default destructor
    ~ChemicalReaction () = default;

    int get_type () 
    { return m_type; }

    const amrex::Vector<int>& get_phases () const
    { return m_phases; }

    const amrex::Vector<std::string>& get_reactants () const
    { return m_reactants; }

    const amrex::Vector<int>& get_reactants_ids () const
    { return m_reactants_IDs; }

    const amrex::Vector<amrex::Real>& get_reactants_coeffs () const
    { return m_reactants_coeffs; }

    const amrex::Vector<int>& get_reactants_phases () const
    { return m_reactants_phases; }

    const amrex::Vector<std::string>& get_products () const
    { return m_products; }

    const amrex::Vector<int>& get_products_ids () const
    { return m_products_IDs; }

    const amrex::Vector<amrex::Real>& get_products_coeffs () const
    { return m_products_coeffs; }

    const amrex::Vector<int>& get_products_phases () const
    { return m_products_phases; }

  private:
    // Reaction type (homogeneous, heterogeneous)
    int m_type;

    // chemical formula
    const std::string m_formula;

    // Phases
    amrex::Vector<int> m_phases;

    // Reactants
    amrex::Vector<std::string> m_reactants;
    amrex::Vector<int> m_reactants_IDs;
    amrex::Vector<amrex::Real> m_reactants_coeffs;
    amrex::Vector<int> m_reactants_phases;

    // Products
    amrex::Vector<std::string> m_products;
    amrex::Vector<int> m_products_IDs;
    amrex::Vector<amrex::Real> m_products_coeffs;
    amrex::Vector<int> m_products_phases;

    // Get reactants data
    std::string parse_reactants(const std::string& formula);

    // Get reactants data
    std::string parse_products(const std::string& formula);

    // Extract stoichiometric data
    void parse_stoichiometric_data(const std::string& s,
                                   amrex::Vector<std::string>& compounds,
                                   amrex::Vector<int>& compounds_id,
                                   amrex::Vector<amrex::Real>& coefficients,
                                   amrex::Vector<int>& phases,
                                   const Species& species);

    // Parse reaction
    void parse_reaction(const Species& species);
};

#endif
