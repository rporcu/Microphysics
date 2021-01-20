#include <mfix.H>
#include <mfix_algorithm.H>

#include <MOL.H>
#include <mfix_mol_convective_fluxes_K.H>
#include <mfix_mol_convective_fluxes_eb_K.H>

using namespace aux;

//
// Compute the three components of the convection term
//
void
mfix::mfix_compute_fluxes (int lev,
                           Vector< MultiFab* >& a_fx,
                           Vector< MultiFab* >& a_fy,
                           Vector< MultiFab* >& a_fz,
                           Vector< MultiFab* > const& state_in,
                           const int state_comp, const int ncomp,
                           Vector< MultiFab* > const& xslopes_in,
                           Vector< MultiFab* > const& yslopes_in,
                           Vector< MultiFab* > const& zslopes_in,
                           const int slopes_comp,
                           Vector< MultiFab* > const& ep_u_mac,
                           Vector< MultiFab* > const& ep_v_mac,
                           Vector< MultiFab* > const& ep_w_mac,
                           const int  nghost,
                           const Real covered_val,
                           const GpuArray<int, 3> bc_types,
                           Array4<int const> const& bct_ilo,
                           Array4<int const> const& bct_ihi,
                           Array4<int const> const& bct_jlo,
                           Array4<int const> const& bct_jhi,
                           Array4<int const> const& bct_klo,
                           Array4<int const> const& bct_khi,
                           EBFArrayBoxFactory const* ebfact,
                           Vector<Geometry> geom)
{
  // Get EB geometric info
  Array< const MultiCutFab*,3> areafrac;
  Array< const MultiCutFab*,3> facecent;
  const amrex::MultiFab*    volfrac;
  const amrex::MultiCutFab* bndrycent;

  areafrac  =   ebfact->getAreaFrac();
  facecent  =   ebfact->getFaceCent();
  volfrac   = &(ebfact->getVolFrac());
  bndrycent = &(ebfact->getBndryCent());

  const auto& cellcent = ebfact->getCentroid();

  // We do this here to avoid any confusion about the FAB setVal.
  a_fx[lev]->setVal(covered_val);
  a_fy[lev]->setVal(covered_val);
  a_fz[lev]->setVal(covered_val);


  for (MFIter mfi(*state_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    // Tilebox
    Box bx = mfi.tilebox ();

    const Box ubx = amrex::surroundingNodes(bx,0);
    const Box vbx = amrex::surroundingNodes(bx,1);
    const Box wbx = amrex::surroundingNodes(bx,2);

    // this is to check efficiently if this tile contains any eb stuff
    const EBFArrayBox& state_fab = static_cast<EBFArrayBox const&>((*state_in[lev])[mfi]);
    const EBCellFlagFab&  flags = state_fab.getEBCellFlagFab();

    if (flags.getType(amrex::grow(bx,0)) != FabType::covered )
    {
      // No cut cells in tile + nghost-cell width halo -> use non-eb routine
      if (flags.getType(amrex::grow(bx,nghost)) == FabType::regular )
      {
        mol::mfix_compute_fluxes_on_box(
               lev, bx, ubx, vbx, wbx,
              (*a_fx[lev])[mfi], (*a_fy[lev])[mfi], (*a_fz[lev])[mfi],
              (*state_in[lev])[mfi], state_comp, ncomp,
              (*xslopes_in[lev])[mfi], (*yslopes_in[lev])[mfi], (*zslopes_in[lev])[mfi], slopes_comp,
               (*ep_u_mac[lev])[mfi], (*ep_v_mac[lev])[mfi], (*ep_w_mac[lev])[mfi],
               bc_types, bct_ilo, bct_ihi, bct_jlo, bct_jhi, bct_klo, bct_khi, geom);
      }
      else
      {
        mol::mfix_compute_eb_fluxes_on_box(
               lev, bx, ubx, vbx, wbx,
              (*a_fx[lev])[mfi], (*a_fy[lev])[mfi], (*a_fz[lev])[mfi],
              (*state_in[lev])[mfi], state_comp, ncomp,
              (*xslopes_in[lev])[mfi], (*yslopes_in[lev])[mfi], (*zslopes_in[lev])[mfi], slopes_comp,
              (*ep_u_mac[lev])[mfi], (*ep_v_mac[lev])[mfi], (*ep_w_mac[lev])[mfi],
              (*areafrac[0])[mfi], (*areafrac[1])[mfi], (*areafrac[2])[mfi],
              (*facecent[0])[mfi], (*facecent[1])[mfi], (*facecent[2])[mfi],
               cellcent[mfi], (*volfrac)[mfi], (*bndrycent)[mfi], flags,
               bc_types, bct_ilo, bct_ihi, bct_jlo, bct_jhi, bct_klo, bct_khi, geom);
      }
    }
  } // MFIter
}
