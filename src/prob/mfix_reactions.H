#ifndef MFIX_REACTIONS_PARMS_H_
#define MFIX_REACTIONS_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>

#include <AMReX_Vector.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Gpu.H>


// Forward declaration
class MFIXChemicalReaction;
class MFIXReactionsParms;


class MFIXReactions {

  public:
    // Constructor
    MFIXReactions();

    // Destructor
    ~MFIXReactions();

    // Initialization procedure
    void Initialize(const MFIXSpecies& species);

    int solve () const { return m_solve; }

    int nreactions () const { return m_nreactions; }

    const std::vector<std::string>& reactions () const { return m_reactions; }

    const std::vector<std::string>& reaction_equations () const { return m_reaction_equations; }

    int isInitialized () const { return m_is_initialized; }

    const MFIXReactionsParms& parameters () const { return *m_parameters; }

    // Getter
    MFIXChemicalReaction* get(const int q) const
    {
      AMREX_ALWAYS_ASSERT(q < m_nreactions);
      return m_chemical_reactions[q];
    }

  private:
    // Switch for turning on/off chemical reactions modeling
    int m_solve;

    // Number of chemical reactions allowed by the model
    int m_nreactions;

    // Names of chemical reactions allowed by the model
    std::vector<std::string> m_reactions;

    // Chemical reactions equations
    std::vector<std::string> m_reaction_equations;

    int m_is_initialized;

    MFIXReactionsParms* m_parameters;

    // Vector of chemical reactions
    amrex::Vector<MFIXChemicalReaction*> m_chemical_reactions;
};


class MFIXReactionsParms {

  public:
    AMREX_GPU_HOST_DEVICE
    MFIXReactionsParms ()
      : m_nreactions(0)
    {}

    AMREX_GPU_HOST_DEVICE
    MFIXReactionsParms (const int nreactions)
      : m_nreactions(nreactions)
    {}

    int m_nreactions;
};


// Class for chemical reaction details
class MFIXChemicalReaction {

  public:
    struct ReactionType { enum {Invalid=-1, Homogeneous, Heterogeneous}; };

    struct ChemicalPhase { enum {Invalid=-1, Fluid, Solid}; };

    // Default constructor
    explicit MFIXChemicalReaction () = default;

    // class constructor
    explicit MFIXChemicalReaction (const std::string& reaction,
                               const MFIXSpecies& species);

    MFIXChemicalReaction (const MFIXChemicalReaction&) = delete;

    MFIXChemicalReaction& operator = (const MFIXChemicalReaction&) = delete;

    // Default destructor
    ~MFIXChemicalReaction () = default;

    int get_type () { return m_type; }

    const amrex::Vector<int>& get_phases () const { return m_phases; }

    const amrex::Vector<std::string>& get_reactants () const { return m_reactants; }

    const amrex::Vector<int>& get_reactants_ids () const { return m_reactants_IDs; }

    const amrex::Vector<amrex::Real>& get_reactants_coeffs () const { return m_reactants_coeffs; }

    const amrex::Vector<int>& get_reactants_phases () const { return m_reactants_phases; }

    const amrex::Vector<std::string>& get_products () const { return m_products; }

    const amrex::Vector<int>& get_products_ids () const { return m_products_IDs; }

    const amrex::Vector<amrex::Real>& get_products_coeffs () const { return m_products_coeffs; }

    const amrex::Vector<int>& get_products_phases () const { return m_products_phases; }

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
    std::string parse_reactants (const std::string& formula);

    // Get reactants data
    std::string parse_products (const std::string& formula);

    // Extract stoichiometric data
    void parse_stoichiometric_data (const std::string& s,
                                    amrex::Vector<std::string>& compounds,
                                    amrex::Vector<int>& compounds_id,
                                    amrex::Vector<amrex::Real>& coefficients,
                                    amrex::Vector<int>& phases,
                                    const MFIXSpecies& species);

    void check_mass_balance (const MFIXSpecies& species);

    // Parse reaction
    void parse_reaction (const MFIXSpecies& species);

    amrex::Real m_mass_balance_tolerance;
};

#endif
