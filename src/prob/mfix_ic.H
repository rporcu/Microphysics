#ifndef MFIX_IC_PARMS_H_
#define MFIX_IC_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>

#include <AMReX_Vector.H>

#include <mfix_regions.H>
#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_pic.H>


// Forward declaration
struct IC_t;


class MFIXInitialConditions
{
  public:
    void Initialize (const MFIXRegions& regions,
                     MFIXFluidPhase& fluid,
                     const MFIXSolidsPhase& solids,
                     MFIXDEM& dem,
                     MFIXPIC& pic);

    const amrex::Vector<IC_t>& ic () const { return m_ic; };
    const IC_t& ic (const int icv) const { return m_ic[icv]; };

    enum struct ParticleInitType {Invalid, Auto, AsciiFile};

    int AutoParticleInit () const {
      return (m_particle_init_type == ParticleInitType::Auto) ? 1 : 0;
    };

    // Store the number of particles generated in the IC region.
    void set_particle_count ( const int icv, const long np_in ) {
      AMREX_ASSERT( icv < m_np.size() );
      m_np[icv] = np_in;
    }

    // Return the number of particles generated in the IC region.
    long get_particle_count ( const int icv ) const {
      AMREX_ASSERT( icv < m_np.size() );
      return m_np[icv];
    }

  private:

    amrex::Vector<IC_t> m_ic;

    ParticleInitType m_particle_init_type = ParticleInitType::Invalid;

    amrex::Vector<long> m_np;

};


struct IC_t {

  // Constructor
  IC_t (MFIXFluidPhase& fluid_in)
    : fluid(fluid_in)
  {}

  const amrex::RealBox* region;

  FLUID_t fluid;

  amrex::Vector<SOLIDS_t> solids;
  std::string packing;

  const SOLIDS_t* get_solid (const int phase) const
  {
    for (int n(0); n < solids.size(); ++n)
      if (phase == solids[n].phase)
        return &(solids[n]);

    amrex::Abort("Error: solid type not found");
    return nullptr;
  }
};

#endif
