#include <mfix_bc_list.H>
#include <AMReX_Geometry.H>

using namespace amrex;


// Constructor
BCList::BCList (const int nlev_in)
  : nlev(nlev_in)
{}


// Destructor
BCList::~BCList()
{
  for (int lev(0); lev < nlev; ++lev) {
    delete bc_ilo[lev];
    delete bc_ihi[lev];
    delete bc_jlo[lev];
    delete bc_jhi[lev];
    delete bc_klo[lev];
    delete bc_khi[lev];
  }
}


void
BCList::MakeBCArrays (int nghost,
                      bool ooo_debug,
                      amrex::Vector<amrex::Geometry>& geom)
{
  for (int lev = 0; lev < bc_ilo.size(); lev++)
  {
    if (bc_ilo[lev] != nullptr) delete bc_ilo[lev];
    if (bc_ihi[lev] != nullptr) delete bc_ihi[lev];
    if (bc_jlo[lev] != nullptr) delete bc_jlo[lev];
    if (bc_jhi[lev] != nullptr) delete bc_jhi[lev];
    if (bc_klo[lev] != nullptr) delete bc_klo[lev];
    if (bc_khi[lev] != nullptr) delete bc_khi[lev];
  }

  if (ooo_debug) amrex::Print() << "MakeBCArrays" << std::endl;

  bc_ilo.clear(); bc_ilo.resize(nlev, nullptr);
  bc_ihi.clear(); bc_ihi.resize(nlev, nullptr);
  bc_jlo.clear(); bc_jlo.resize(nlev, nullptr);
  bc_jhi.clear(); bc_jhi.resize(nlev, nullptr);
  bc_klo.clear(); bc_klo.resize(nlev, nullptr);
  bc_khi.clear(); bc_khi.resize(nlev, nullptr);

  for (int lev = 0; lev < nlev; lev++)
  {
    // Define and allocate the integer MultiFab that is the outside adjacent
    // cells of the problem domain.
    Box domainx(geom[lev].Domain());
    domainx.grow(1,nghost);
    domainx.grow(2,nghost);
    Box box_ilo = amrex::adjCellLo(domainx,0,1);
    Box box_ihi = amrex::adjCellHi(domainx,0,1);

    Box domainy(geom[lev].Domain());
    domainy.grow(0,nghost);
    domainy.grow(2,nghost);
    Box box_jlo = amrex::adjCellLo(domainy,1,1);
    Box box_jhi = amrex::adjCellHi(domainy,1,1);

    Box domainz(geom[lev].Domain());
    domainz.grow(0,nghost);
    domainz.grow(1,nghost);
    Box box_klo = amrex::adjCellLo(domainz,2,1);
    Box box_khi = amrex::adjCellHi(domainz,2,1);

    // Note that each of these is a single IArrayBox so every process has a copy of them
    bc_ilo[lev] = new IArrayBox(box_ilo,2);
    bc_ihi[lev] = new IArrayBox(box_ihi,2);
    bc_jlo[lev] = new IArrayBox(box_jlo,2);
    bc_jhi[lev] = new IArrayBox(box_jhi,2);
    bc_klo[lev] = new IArrayBox(box_klo,2);
    bc_khi[lev] = new IArrayBox(box_khi,2);
  }
}
