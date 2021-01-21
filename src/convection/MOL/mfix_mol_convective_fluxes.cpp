#include <mfix_algorithm.H>

#include <MOL.H>
#include <mfix_mol_convective_fluxes_K.H>
#include <mfix_mol_convective_fluxes_eb_K.H>

using namespace aux;

//
// Compute the three components of the convection term
//
void
mol::mfix_compute_fluxes (const int lev,
                          Vector< MultiFab* >& a_fx,
                          Vector< MultiFab* >& a_fy,
                          Vector< MultiFab* >& a_fz,
                          Vector< MultiFab* > const& state_in,
                          const int state_comp, const int ncomp,
                          Vector< MultiFab* > const& ep_u_mac,
                          Vector< MultiFab* > const& ep_v_mac,
                          Vector< MultiFab* > const& ep_w_mac,
                          const int  nghost,
                          const Real covered_val,
                          const GpuArray<int, 2> bc_types,
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
  Array< const MultiCutFab*,3> areafrac  =   ebfact->getAreaFrac();
  Array< const MultiCutFab*,3> facecent  =   ebfact->getFaceCent();
  const amrex::MultiFab*       volfrac   = &(ebfact->getVolFrac());
  const amrex::MultiCutFab*    bndrycent = &(ebfact->getBndryCent());

  auto const& flags = ebfact->getMultiEBCellFlagFab();

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
        // Face centroids
        const auto& fcx_fab = facecent[0]->const_array(mfi);
        const auto& fcy_fab = facecent[1]->const_array(mfi);
        const auto& fcz_fab = facecent[2]->const_array(mfi);

        // Cell centroids
        const auto& ccc_fab = cellcent.const_array(mfi);

        mol::mfix_compute_eb_fluxes_on_box(
               lev, domain_bx, xbx, ybx, zbx, ncomp, state_comp, state_fab,
               fx_fab, fy_fab, fz_fab, ep_u_mac_fab, ep_v_mac_fab, ep_w_mac_fab,
               flagarr, fcx_fab, fcy_fab, fcz_fab, ccc_fab,
               bc_types, bct_ilo, bct_ihi, bct_jlo, bct_jhi, bct_klo, bct_khi, geom);
      }
    }
  } // MFIter
}
