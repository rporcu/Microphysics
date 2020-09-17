#ifndef MFIX_REACTIONS_PARMS_H_
#define MFIX_REACTIONS_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>

#include <AMReX_Vector.H>
#include <AMReX_GpuContainers.H>

struct CHEMICALPHASE { enum {Invalid=-1, Fluid, Solid}; };

struct REACTIONTYPE { enum {Invalid=-1, Homogeneous, Heterogeneous}; };

namespace REACTIONS
{
  constexpr int NMAX = 10;

  // Switch for turning on/off chemical reactions modeling
  extern bool solve;

  // Names of chemical reactions allowed by the model
  extern std::vector<std::string> reactions;

  // Number of chemical reactions allowed by the model
  extern int nreactions;

  // Chemical reactions equations
  extern std::vector<std::string> reaction_equations;

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

      // chemical formula
      const std::string m_reaction;

      amrex::Gpu::ManagedVector<int> m_phases;

      // Reaction type (homogeneous, heterogeneous)
      int m_reaction_type;

      // Reactants
      amrex::Vector<std::string> m_reactants;
      amrex::Gpu::ManagedVector<int> m_reactants_id;
      amrex::Gpu::ManagedVector<amrex::Real> m_reactants_coeffs;
      amrex::Gpu::ManagedVector<int> m_reactants_phases;

      // Products
      amrex::Vector<std::string> m_products;
      amrex::Gpu::ManagedVector<int> m_products_id;
      amrex::Gpu::ManagedVector<amrex::Real> m_products_coeffs;
      amrex::Gpu::ManagedVector<int> m_products_phases;
  };

  void Initialize();

}

#endif
