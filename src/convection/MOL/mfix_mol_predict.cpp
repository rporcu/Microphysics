#include <MOL.H>

#include <mfix_mol_predict_vels_on_faces_K.H>
#include <mfix_mol_predict_vels_on_faces_eb_K.H>

void
mol::predict_vels_on_faces (int lev,
                            MultiFab& ep_u_mac,  MultiFab& ep_v_mac,  MultiFab& ep_w_mac,
                            MultiFab& ep_face_x, MultiFab& ep_face_y, MultiFab& ep_face_z,
                            MultiFab& vel_in,
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
    const auto& umac_fab = ep_u_mac.array(mfi);
    const auto& vmac_fab = ep_v_mac.array(mfi);
    const auto& wmac_fab = ep_w_mac.array(mfi);

    // Cell-centered velocity
    const auto& ccvel_fab = vel_in.const_array(mfi);

    // Face-centered ep
    const auto& epx_fab = ep_face_x.const_array(mfi);
    const auto& epy_fab = ep_face_y.const_array(mfi);
    const auto& epz_fab = ep_face_z.const_array(mfi);

    EBCellFlagFab const& flagfab = flags[mfi];
    Array4<EBCellFlag const> const& flagarr = flagfab.const_array();
    const Box& bx = mfi.tilebox();
    auto const typ = flagfab.getType(amrex::grow(bx,2));

    if (typ == FabType::covered )
    {
      amrex::ParallelFor(ubx, vbx, wbx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { umac_fab(i,j,k) = 0.0; },
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { vmac_fab(i,j,k) = 0.0; },
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { wmac_fab(i,j,k) = 0.0; });
    }
    // Cut cells in this FAB

    else if (typ == FabType::singlevalued)
    {

      // Face centroids
      const auto& fcx_fab = facecent[0]->const_array(mfi);
      const auto& fcy_fab = facecent[1]->const_array(mfi);
      const auto& fcz_fab = facecent[2]->const_array(mfi);

      // Cell centroids
      const auto& ccc_fab = cellcent.const_array(mfi);

      mol::predict_vels_on_faces_eb(domain_bx, ubx, vbx, wbx,
                                    umac_fab, vmac_fab, wmac_fab, ccvel_fab,
                                    epx_fab, epy_fab, epz_fab,
                                    flagarr, fcx_fab, fcy_fab,  fcz_fab,
                                    ccc_fab, bc_types,
                                    bc_ilo, bc_ihi,
                                    bc_jlo, bc_jhi,
                                    bc_klo, bc_khi);
    }
    // No cut cells in this FAB
    else if (typ == FabType::regular )
      {

        mol::predict_vels_on_faces(domain_bx, ubx, vbx, wbx,
                                   umac_fab, vmac_fab, wmac_fab, ccvel_fab,
                                   epx_fab,  epy_fab,  epz_fab,
                                   bc_types,
                                   bc_ilo, bc_ihi,
                                   bc_jlo, bc_jhi,
                                   bc_klo, bc_khi);

    }

  } // MFIter

}
