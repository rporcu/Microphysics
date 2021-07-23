#ifndef MFIX_REACTIONS_PARMS_H_
#define MFIX_REACTIONS_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>

#include <AMReX_Vector.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Gpu.H>

#include <mfix_fluid_parms.H>
#include <mfix_solids_parms.H>


// Forward declaration
class ChemicalReaction;


class ReactionsParms {

  public:
    const int m_nreactions;
    const int* m_nreactants;
    const int** m_reactants_id;
    const amrex::Real** m_reactants_coeffs;
    const int** m_reactants_phases;
    const int* m_nproducts;
    const int** m_products_id;
    const amrex::Real** m_products_coeffs;
    const int** m_products_phases;

  public:
    AMREX_GPU_HOST_DEVICE
    ReactionsParms ()
      : m_nreactions(0)
      , m_nreactants(nullptr)
      , m_reactants_id(nullptr)
      , m_reactants_coeffs(nullptr)
      , m_reactants_phases(nullptr)
      , m_nproducts(nullptr)
      , m_products_id(nullptr)
      , m_products_coeffs(nullptr)
      , m_products_phases(nullptr)
    {}

    AMREX_GPU_HOST_DEVICE
    ReactionsParms (const int nreactions,
                    const int* nreactants,
                    const int** reactants_id,
                    const amrex::Real** reactants_coeffs,
                    const int** reactants_phases,
                    const int* nproducts,
                    const int** products_id,
                    const amrex::Real** products_coeffs,
                    const int** products_phases)
      : m_nreactions(nreactions)
      , m_nreactants(nreactants)
      , m_reactants_id(reactants_id)
      , m_reactants_coeffs(reactants_coeffs)
      , m_reactants_phases(reactants_phases)
      , m_nproducts(nproducts)
      , m_products_id(products_id)
      , m_products_coeffs(products_coeffs)
      , m_products_phases(products_phases)
    {}

    AMREX_GPU_HOST_DEVICE
    ~ReactionsParms() {}

    // Compute heterogeneous reactions
    template <class HeterogeneousRatesFunc>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void heterogeneous_rates(amrex::Real* R_q,
                             const SolidsParms& solids_parms,
                             const amrex::Real* X_sn,
                             const amrex::Real ro_s,
                             const amrex::Real ep_s,
                             const FluidParms& fluid_parms,
                             const amrex::Real* X_gk,
                             const amrex::Real ro_g,
                             const amrex::Real ep_g,
                             const HeterogeneousRatesFunc& heterogeneous_rates_func) const
    { heterogeneous_rates_func(R_q, *this, solids_parms, X_sn, ro_s, ep_s, fluid_parms, X_gk, ro_g, ep_g); }
    
    // Compute heterogeneous reactions
    template <class HomogeneousRatesFunc>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void homogeneous_rates(amrex::Real* R_q,
                           const SolidsParms& solids_parms,
                           const amrex::Real* X_sn,
                           const amrex::Real ro_s,
                           const amrex::Real ep_s,
                           const HomogeneousRatesFunc& homogeneous_rates_func) const
    { homogeneous_rates_func(R_q, *this, solids_parms, X_sn, ro_s, ep_s); }
    
    // Compute heterogeneous reactions
    template <class HomogeneousRatesFunc>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void homogeneous_rates(amrex::Real* R_q,
                           const FluidParms& fluid_parms,
                           const amrex::Real* X_sn,
                           const amrex::Real ro_s,
                           const amrex::Real ep_s,
                           const HomogeneousRatesFunc& homogeneous_rates_func) const
    { homogeneous_rates_func(R_q, *this, fluid_parms, X_sn, ro_s, ep_s); }
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
    bool solve;

    // Number of chemical reactions allowed by the model
    int nreactions;

    // Names of chemical reactions allowed by the model
    std::vector<std::string> reactions;

    // Chemical reactions equations
    std::vector<std::string> reaction_equations;

    amrex::Gpu::DeviceVector<int> d_types;

    amrex::Gpu::DeviceVector<int> d_nphases;

    amrex::Gpu::DeviceVector<const int*> d_phases;

    amrex::Gpu::DeviceVector<int> d_nreactants;
    amrex::Gpu::DeviceVector<const int*> d_reactants_id;
    amrex::Gpu::DeviceVector<const amrex::Real*> d_reactants_coeffs;
    amrex::Gpu::DeviceVector<const int*> d_reactants_phases;

    amrex::Gpu::DeviceVector<int> d_nproducts;
    amrex::Gpu::DeviceVector<const int*> d_products_id;
    amrex::Gpu::DeviceVector<const amrex::Real*> d_products_coeffs;
    amrex::Gpu::DeviceVector<const int*> d_products_phases;

    // Initialization procedure
    void Initialize();

    // Getter
    ChemicalReaction* get(const int q) const
    {
      AMREX_ALWAYS_ASSERT(q < nreactions);
      return m_chemical_reactions[q];
    }

    ReactionsParms* parameters;

  private:
    // Vector of chemical reactions
    amrex::Vector< ChemicalReaction* > m_chemical_reactions;
};


// Class for chemical reaction details
class ChemicalReaction {

  public:
    // Default constructor
    explicit ChemicalReaction () = default;

    // class constructor
    explicit ChemicalReaction (const std::string& reaction);

    ChemicalReaction (const ChemicalReaction&) = delete;

    ChemicalReaction& operator = (const ChemicalReaction&) = delete;

    // Default destructor
    ~ChemicalReaction () = default;

    struct REACTIONTYPE { enum {Invalid=-1, Homogeneous, Heterogeneous}; };

    // Reaction type (homogeneous, heterogeneous)
    int m_reaction_type;

    struct CHEMICALPHASE { enum {Invalid=-1, Fluid, Solid}; };

    // Phases
    amrex::Gpu::DeviceVector<int> m_phases;

    // chemical formula
    const std::string m_reaction;

    // Reactants
    amrex::Vector<std::string> m_reactants;
    amrex::Gpu::DeviceVector<int> m_reactants_id;
    amrex::Gpu::DeviceVector<amrex::Real> m_reactants_coeffs;
    amrex::Gpu::DeviceVector<int> m_reactants_phases;

    // Products
    amrex::Vector<std::string> m_products;
    amrex::Gpu::DeviceVector<int> m_products_id;
    amrex::Gpu::DeviceVector<amrex::Real> m_products_coeffs;
    amrex::Gpu::DeviceVector<int> m_products_phases;

  private:
    // Get reactants data
    std::string get_reactants(const std::string& formula);

    // Get reactants data
    std::string get_products(const std::string& formula);

    // Extract stoichiometric data
    void get_stoichiometric_data(const std::string& s,
                                 std::vector<std::string>& compounds,
                                 amrex::Gpu::HostVector<int>& compounds_id,
                                 amrex::Gpu::HostVector<amrex::Real>& coefficients,
                                 amrex::Gpu::HostVector<int>& phases);

    // Parse reaction
    void parse_reaction(const std::string& equation,
                        amrex::Gpu::HostVector<int>& phases,
                        int& reaction_type,
                        std::vector<std::string>& reactants,
                        amrex::Gpu::HostVector<int>& reactants_id,
                        amrex::Gpu::HostVector<amrex::Real>& reactants_coeffs,
                        amrex::Gpu::HostVector<int>& reactants_phases,
                        std::vector<std::string>& products,
                        amrex::Gpu::HostVector<int>& products_id,
                        amrex::Gpu::HostVector<amrex::Real>& products_coeffs,
                        amrex::Gpu::HostVector<int>& products_phases);
};

#endif
