#include <mfix_proj_F.H>
#include <diffusion_F.H>
#include <mfix.H>

#include <AMReX_NodalProjector.H>
#include <DiffusionOp.H>
#include <MFIX_BC_Parms.H>

void
mfix::mfix_init_solvers ()
{
    BL_PROFILE("mfix::mfix_init_solvers");

    int bc_lo[3], bc_hi[3];
    Box domain(geom[0].Domain());

    //

    diffusion_op.reset(new DiffusionOp(this, &ebfactory,
                                       BC::diff_vel_lobc,  BC::diff_vel_hibc,
                                       BC::diff_scal_lobc, BC::diff_scal_hibc, nghost));
}

void
mfix::mfix_setup_solvers ()
{
    BL_PROFILE("mfix::mfix_setup_solvers");

    int bc_lo[3], bc_hi[3];
    Box domain(geom[0].Domain());

    //
    // Now the diffusion solver
    //
    diffusion_op->setup(this, &ebfactory);
}
