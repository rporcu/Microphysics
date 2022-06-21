#ifndef MFIX_REGIONS_PARMS_H_
#define MFIX_REGIONS_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>

// Name and extents of regions
class Regions
{
  public:
    // Constructor
    Regions ()
      : m_is_initialized(0)
    {}

    // Initialization function -- read input parameters
    void Initialize ();

    // Returns the index of region with matching name.
    // -1 is returned if the region is not matched.
    int get_index(const std::string& input) const;

    // Return region extents by region name
    const amrex::RealBox* get_region(const std::string& name) const;

    int is_initialized() const
    { return m_is_initialized; }

  private:
    amrex::Vector<std::string> m_names;

    amrex::Vector<amrex::RealBox> m_extents;

    int m_is_initialized;
};

#endif
