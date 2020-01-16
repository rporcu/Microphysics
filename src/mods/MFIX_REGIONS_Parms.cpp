#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>
#include <AMReX_RealBox.H>

#include <AMReX_ParmParse.H>
#include <MFIX_REGIONS_Parms.H>

namespace REGIONS
{

  amrex::Vector<std::string> names;
  amrex::Vector<amrex::RealBox> extents;

  void Initialize ()
  {

    // Read in the region data
    amrex::ParmParse pp("mfix");

    std::vector<std::string> regions_in;
    pp.queryarr("regions", regions_in);

    for(int lc=0; lc<regions_in.size(); lc++){

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

        names.push_back(regions_in[lc]);
        extents.push_back(new_box);

      }


    }

#if 0

    if( REGIONS::names.size() > 0 ) {
      amrex::Print() << std::endl << "Found regions!"  << std::endl;
      for (int lc=0; lc < REGIONS::names.size(); ++lc ) {

        amrex::RealBox myBox = REGIONS::extents[lc];

        amrex::Print() <<
          std::endl << " Index: " << lc <<
          std::endl << "  Name: " << REGIONS::names[lc] <<
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
  int getIndex(std::string name){

    int idx=-1;
    for(int lc=0; lc<names.size(); lc++){
      if( names[lc] == name){
        return lc;
      }
    }
    return -1; // Failed to find a match
  }



  /*-----------------------------------------------------------------\
  |  Given 'name', return a pointer to the RealBox of the region     |
  |  with the matching name.                                         |
  \-----------------------------------------------------------------*/
  const amrex::RealBox* getRegion(std::string name){

    for(int lc=0; lc<names.size(); lc++){
      if( names[lc] == name){
        return &extents[lc];
      }
    }
    return NULL;
  }


}
