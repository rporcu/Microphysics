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
mfix::mfix_compute_fluxes (const int lev,
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

  auto const& flags = ebfact->getMultiEBCellFlagFab();

  areafrac  =   ebfact->getAreaFrac();
  facecent  =   ebfact->getFaceCent();
  volfrac   = &(ebfact->getVolFrac());
  bndrycent = &(ebfact->getBndryCent());

  const auto& cellcent = ebfact->getCentroid();

  // We do this here to avoid any confusion about the FAB setVal.
  a_fx[lev]->setVal(covered_val);
  a_fy[lev]->setVal(covered_val);
  a_fz[lev]->setVal(covered_val);

  const Box& domain_bx = geom[lev].Domain();

  for (MFIter mfi(*state_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    // Tilebox
    Box bx = mfi.tilebox ();

    const Box xbx = amrex::surroundingNodes(bx,0);
    const Box ybx = amrex::surroundingNodes(bx,1);
    const Box zbx = amrex::surroundingNodes(bx,2);

    Array4<Real const> const& state_fab = state_in[lev]->const_array(mfi);

    // Face-centered MAC velocities (ep*u_mac, ep*v_mac, ep*w_mac)
    Array4<Real const> const& ep_u_mac_fab = ep_u_mac[lev]->const_array(mfi);
    Array4<Real const> const& ep_v_mac_fab = ep_v_mac[lev]->const_array(mfi);
    Array4<Real const> const& ep_w_mac_fab = ep_w_mac[lev]->const_array(mfi);

    // Face-centered fluxes
    Array4<Real      > const& fx_fab = a_fx[lev]->array(mfi);
    Array4<Real      > const& fy_fab = a_fy[lev]->array(mfi);
    Array4<Real      > const& fz_fab = a_fz[lev]->array(mfi);

    // this is to check efficiently if this tile contains any eb stuff
    EBCellFlagFab const& flagfab = flags[mfi];
    Array4<EBCellFlag const> const& flagarr = flagfab.const_array();

    auto const typ = flagfab.getType(amrex::grow(bx,2));

    if (flagfab.getType(amrex::grow(bx,0)) != FabType::covered )
    {
      // No cut cells in tile + nghost-cell width halo -> use non-eb routine
      if (flagfab.getType(amrex::grow(bx,nghost)) == FabType::regular )
      {
        mol::mfix_compute_fluxes_on_box(
               lev, domain_bx, xbx, ybx, zbx, ncomp, state_comp, state_fab,
               fx_fab, fy_fab, fz_fab, ep_u_mac_fab, ep_v_mac_fab, ep_w_mac_fab,
               bc_types, bct_ilo, bct_ihi, bct_jlo, bct_jhi, bct_klo, bct_khi);
      }
      else
      {
        mol::mfix_compute_eb_fluxes_on_box(
               lev, bx, xbx, ybx, zbx,
              (*a_fx[lev])[mfi], (*a_fy[lev])[mfi], (*a_fz[lev])[mfi],
              (*state_in[lev])[mfi], state_comp, ncomp,
              (*xslopes_in[lev])[mfi], (*yslopes_in[lev])[mfi], (*zslopes_in[lev])[mfi], slopes_comp,
              (*ep_u_mac[lev])[mfi], (*ep_v_mac[lev])[mfi], (*ep_w_mac[lev])[mfi],
              (*areafrac[0])[mfi], (*areafrac[1])[mfi], (*areafrac[2])[mfi],
              (*facecent[0])[mfi], (*facecent[1])[mfi], (*facecent[2])[mfi],
               cellcent[mfi], (*volfrac)[mfi], (*bndrycent)[mfi],
               bc_types, bct_ilo, bct_ihi, bct_jlo, bct_jhi, bct_klo, bct_khi, geom);
      }
    }
  } // MFIter
}
