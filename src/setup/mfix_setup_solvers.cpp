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
}

void
mfix::mfix_setup_solvers ()
{
    BL_PROFILE("mfix::mfix_setup_solvers");

    diffusion_op->setup(this, &ebfactory);
}
