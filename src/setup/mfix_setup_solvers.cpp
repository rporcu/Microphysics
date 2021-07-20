#include <mfix.H>

#include <AMReX_NodalProjector.H>
#include <mfix_diffusion_op.H>
#include <mfix_bc_parms.H>

void
mfix::mfix_init_solvers ()
{
    BL_PROFILE("mfix::mfix_init_solvers");

    diffusion_op = std::make_unique<DiffusionOp>(this, amrex::GetVecOfConstPtrs(ebfactory), fluid,
                                       BC::diff_vel_lobc,         BC::diff_vel_hibc,
                                       BC::diff_scal_lobc,        BC::diff_scal_hibc,
                                       BC::diff_temperature_lobc, BC::diff_temperature_hibc,
                                       BC::diff_species_lobc,     BC::diff_species_hibc,
                                       nghost_state());

    macproj = std::make_unique<MacProjector>(Geom(0,finest_level),
                                   MLMG::Location::FaceCentroid,  // Location of mac_vec
                                   MLMG::Location::FaceCentroid,  // Location of beta
                                   MLMG::Location::CellCenter,    // Location of solution variable phi
                                   MLMG::Location::CellCentroid);// Location of MAC RHS
}

void
mfix::mfix_setup_solvers ()
{
    BL_PROFILE("mfix::mfix_setup_solvers");

    diffusion_op->setup(this, amrex::GetVecOfConstPtrs(ebfactory));
}
