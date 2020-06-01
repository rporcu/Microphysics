#include <mfix.H>

#include <AMReX_NodalProjector.H>
#include <DiffusionOp.H>
#include <MFIX_BC_Parms.H>

void
mfix::mfix_init_solvers ()
{
    BL_PROFILE("mfix::mfix_init_solvers");

    diffusion_op.reset(new DiffusionOp(this, &ebfactory,
                                       BC::diff_vel_lobc,         BC::diff_vel_hibc,
                                       BC::diff_scal_lobc,        BC::diff_scal_hibc,
                                       BC::diff_temperature_lobc, BC::diff_temperature_hibc,
                                       BC::diff_species_lobc,     BC::diff_species_hibc,
                                       nghost));

    // Initialize MultiFabs to be used in the MLMG solvers in DiffusionOp for
    // mixed inhomogeneous Dirichlet - homogeneous Neumann boundary conditions
    // on some of the EB regions
    if (advect_enthalpy)
    {
      bool eb_temperature_is_dirichlet = false;

      for (int bcv(0); bcv < BC::bc.size(); bcv++) {
        // If at least one BC condition is eb_
        if ( BC::bc[bcv].type == bc_list.get_eb() ) {
          eb_temperature_is_dirichlet = true;
          break;
        }
      }

      if ( eb_temperature_is_dirichlet ) {
        Vector< std::unique_ptr< MultiFab > > T_g_on_eb(finest_level+1);
        Vector< std::unique_ptr< MultiFab > > k_g_on_eb(finest_level+1);

        for(int lev = 0; lev <= max_level; lev++) {
          T_g_on_eb[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));
          k_g_on_eb[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));
        }

        // Assign Dirichlet values to the MultiFabs on the proper EB region
        // where we want to set inhomogeneous Dirichlet BCs
        this->mfix_set_eb_temperature_bcs(T_g_on_eb, k_g_on_eb);

        // Transfer the ownership of the MultiFabs allocated here to the
        // DiffusionOp class
        diffusion_op->setup_eb_temperature(T_g_on_eb, k_g_on_eb);
      }
    }
}

void
mfix::mfix_setup_solvers ()
{
    BL_PROFILE("mfix::mfix_setup_solvers");

    diffusion_op->setup(this, &ebfactory);
}
