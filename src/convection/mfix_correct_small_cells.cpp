#include <mfix.H>

void
mfix::mfix_correct_small_cells (Vector<MultiFab*      > const& vel_in,
                                Vector<MultiFab const*> const& ep_u_mac,
                                Vector<MultiFab const*> const& ep_v_mac,
                                Vector<MultiFab const*> const& ep_w_mac)
{
  BL_PROFILE("mfix::mfix_correct_small_cells");

  // Get EB geometric info
  Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;

  for (int lev = 0; lev < nlev; lev++)
  {
    areafrac = ebfactory[lev]->getAreaFrac();

    for (MFIter mfi(*vel_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      // Tilebox
      const Box bx = mfi.tilebox();

      const EBFArrayBox& vel_fab = static_cast<EBFArrayBox const&>((*vel_in[lev])[mfi]);
      const EBCellFlagFab& flags = vel_fab.getEBCellFlagFab();

      // Face-centered velocity components
      const auto& umac_fab = (ep_u_mac[lev])->array(mfi);
      const auto& vmac_fab = (ep_v_mac[lev])->array(mfi);
      const auto& wmac_fab = (ep_w_mac[lev])->array(mfi);

      if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
      {
        // do nothing
      }

      // No cut cells in this FAB
      else if (flags.getType(amrex::grow(bx,1)) == FabType::regular )
      {
        // do nothing
      }

      // Cut cells in this FAB
      else
      {
        const MultiFab* volfrac = &(ebfactory[lev] -> getVolFrac());

        // Face-centered areas
        const auto& apx_fab = areafrac[0]->array(mfi);
        const auto& apy_fab = areafrac[1]->array(mfi);
        const auto& apz_fab = areafrac[2]->array(mfi);

        const auto& vfrac_fab = volfrac->array(mfi);
        const auto& ccvel_fab = vel_in[lev]->array(mfi);

        // This FAB has cut cells -- we define the centroid value in terms
        // of the MAC velocities onfaces
        amrex::ParallelFor(bx,
          [vfrac_fab,apx_fab,apy_fab,apz_fab,ccvel_fab,umac_fab,vmac_fab,wmac_fab]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          Real vfrac = vfrac_fab(i,j,k);
          if (vfrac > 0.0 && vfrac < 1.e-4)
          {
            const Real apx_mns = apx_fab(i,j,k);
            const Real apx_pls = apx_fab(i+1,j,k);

            const Real apy_mns = apy_fab(i,j,k);
            const Real apy_pls = apy_fab(i,j+1,k);

            const Real apz_mns = apz_fab(i,j,k);
            const Real apz_pls = apz_fab(i,j,k+1);

            ccvel_fab(i,j,k,0) =
              (apx_mns*umac_fab(i,j,k)+apx_pls*umac_fab(i+1,j,k))/(apx_mns+apx_pls);
            ccvel_fab(i,j,k,1) =
              (apy_mns*vmac_fab(i,j,k)+apy_pls*vmac_fab(i,j+1,k))/(apy_mns+apy_pls);
            ccvel_fab(i,j,k,2) =
              (apz_mns*wmac_fab(i,j,k)+apz_pls*wmac_fab(i,j,k+1))/(apz_mns+apz_pls);
          }
        });
      }
    }
  }
}
