#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>
#include <AMReX_RealBox.H>

#include <AMReX_ParmParse.H>
#include <mfix_regions_parms.H>

void
Regions::Initialize ()
{
  m_is_initialized = 1;

  // Read in the region data
  amrex::ParmParse pp("mfix");

  std::vector<std::string> regions_in;
  pp.queryarr("regions", regions_in);

  for(size_t lc=0; lc < regions_in.size(); lc++){

    amrex::RealBox new_box;

    amrex::Vector<amrex::Real> reg_lo(3);
    amrex::Vector<amrex::Real> reg_hi(3);

    {
      std::string tmp_input = "regions."+regions_in[lc];
      amrex::ParmParse ptmp(tmp_input.c_str());

      ptmp.getarr("lo", reg_lo, 0, 3);
      ptmp.getarr("hi", reg_hi, 0, 3);

      new_box.setLo(reg_lo);
      new_box.setHi(reg_hi);

      m_names.push_back(regions_in[lc]);
      m_extents.push_back(new_box);
    }
  }

#if 0

  if(m_names.size() > 0 ) {
    amrex::Print() << std::endl << "Found regions!"  << std::endl;
    for (int lc=0; lc < m_names.size(); ++lc ) {

      amrex::RealBox myBox = m_extents[lc];

      amrex::Print() <<
        std::endl << " Index: " << lc <<
        std::endl << "  Name: " << m_names[lc] <<
        std::endl << "    Lo: " << myBox.lo(0) << "  " << myBox.lo(1) << "  " << myBox.lo(2) <<
        std::endl << "    Hi: " << myBox.hi(0) << "  " << myBox.hi(1) << "  " << myBox.hi(2) <<
        std::endl;
    }

  } else {
    amrex::Print() << std::endl << "No regions found :(" << std::endl;
  }
#endif

} // End Initialize


/*-----------------------------------------------------------------\
|  Given 'name', return the index of the index of the region with  |
|  the matching name.                                              |
\-----------------------------------------------------------------*/
int
Regions::get_index (const std::string& name) const
{
  for (int lc=0; lc < m_names.size(); lc++) {
    if (name.compare(m_names[lc]) == 0) {
      return lc;
    }
  }

  return -1; // Failed to find a match
}


/*-----------------------------------------------------------------\
|  Given 'name', return a pointer to the RealBox of the region     |
|  with the matching name.                                         |
\-----------------------------------------------------------------*/
const amrex::RealBox*
Regions::get_region (const std::string& name) const
{
  for (int lc=0; lc < m_names.size(); lc++) {
    if (name.compare(m_names[lc]) == 0) {
      return &m_extents[lc];
    }
  }

  return nullptr;
}
