#include <MOL.H>

#include <mfix_mol_predict_vels_on_faces_K.H>
#include <mfix_mol_predict_vels_on_faces_eb_K.H>

void
mol::predict_vels_on_faces (int lev,
                            MultiFab& ep_u_mac,  MultiFab& ep_v_mac,  MultiFab& ep_w_mac,
                            const MultiFab& vel_in,
                            std::map<std::string, Gpu::DeviceVector<int>>& bc_types,
                            Array4<int const> const& bc_ilo,
                            Array4<int const> const& bc_ihi,
                            Array4<int const> const& bc_jlo,
                            Array4<int const> const& bc_jhi,
                            Array4<int const> const& bc_klo,
                            Array4<int const> const& bc_khi,
                            EBFArrayBoxFactory const* ebfact,
                            Vector<Geometry> geom)

{
  BL_PROFILE("mol::predict_vels_on_faces");

  auto const& flags = ebfact->getMultiEBCellFlagFab();
  auto const& facecent = ebfact->getFaceCent();
  auto const& cellcent = ebfact->getCentroid();
  auto const& vfrac    = ebfact->getVolFrac();

  // ****************************************************************************
  // Then predict to face centers
  // ****************************************************************************

  for (MFIter mfi(vel_in,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& domain_bx = geom[lev].Domain();

    const Box& ubx = mfi.nodaltilebox(0);
    const Box& vbx = mfi.nodaltilebox(1);
    const Box& wbx = mfi.nodaltilebox(2);

    // Face-centered velocity components
    const auto& umac_arr = ep_u_mac.array(mfi);
    const auto& vmac_arr = ep_v_mac.array(mfi);
    const auto& wmac_arr = ep_w_mac.array(mfi);

    // Cell-centered velocity
    const auto& ccvel_arr = vel_in.const_array(mfi);

    EBCellFlagFab const& flagfab = flags[mfi];
    Array4<EBCellFlag const> const& flagarr = flagfab.const_array();
    const Box& bx = mfi.tilebox();
    auto const typ = flagfab.getType(amrex::grow(bx,2));

    if (typ == FabType::covered )
    {
      amrex::ParallelFor(ubx, vbx, wbx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { umac_arr(i,j,k) = 0.0; },
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { vmac_arr(i,j,k) = 0.0; },
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { wmac_arr(i,j,k) = 0.0; });
    }
    // Cut cells in this FAB

    else if (typ == FabType::singlevalued)
    {

      // Face centroids
      const auto& fcx_arr = facecent[0]->const_array(mfi);
      const auto& fcy_arr = facecent[1]->const_array(mfi);
      const auto& fcz_arr = facecent[2]->const_array(mfi);

      // Cell centroids
      const auto& ccc_arr = cellcent.const_array(mfi);

      const auto& vfrac_arr = vfrac.const_array(mfi);

      mol::predict_vels_on_faces_eb(domain_bx, ubx, vbx, wbx,
                                    umac_arr, vmac_arr, wmac_arr, ccvel_arr,
                                    flagarr, fcx_arr, fcy_arr,  fcz_arr,
                                    ccc_arr, vfrac_arr, bc_types,
                                    bc_ilo, bc_ihi,
                                    bc_jlo, bc_jhi,
                                    bc_klo, bc_khi);
    }
    // No cut cells in this FAB
    else if (typ == FabType::regular )
      {

        mol::predict_vels_on_faces(domain_bx, ubx, vbx, wbx,
                                   umac_arr, vmac_arr, wmac_arr, ccvel_arr,
                                   bc_types,
                                   bc_ilo, bc_ihi,
                                   bc_jlo, bc_jhi,
                                   bc_klo, bc_khi);

    }

  } // MFIter

}
