#include <mfix_proj_F.H>
#include <mfix.H>

#include <MFIX_NodalProjection.H>

void
mfix::mfix_setup_nodal_solver ()
{
    BL_PROFILE("mfix::mfix_setup_nodal_solver");

    // Set domain BCs for Poisson's solver
    // The domain BCs refer to level 0 only
    int bc_lo[3], bc_hi[3];
    Box domain(geom[0].Domain());

    set_ppe_bcs(bc_lo, bc_hi,
                domain.loVect(), domain.hiVect(),
                &nghost,
                bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
                bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
                bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

    ppe_lobc = {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]};
    ppe_hibc = {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]};

    nodal_projector.reset(new NodalProjection(this, ppe_lobc, ppe_hibc));

    // //
    // // First define the matrix (operator).
    // //
    // //        (del dot b sigma grad)) phi
    // //
    // LPInfo                       info;
    // info.setMaxCoarseningLevel(nodal_mg_max_coarsening_level);
    // nodal_matrix.reset(new MLNodeLaplacian( geom, grids, dmap, info,
    //                                         GetVecOfConstPtrs(ebfactory)));

    // nodal_matrix->setGaussSeidel(true);
    // nodal_matrix->setHarmonicAverage(false);
    // nodal_matrix->setDomainBC ( ppe_lobc, ppe_hibc);

    // //
    // // Then setup the solver ----------------------
    // //
    // nodal_solver.reset(new MLMG(*nodal_matrix));

    // nodal_solver->setMaxIter    (nodal_mg_maxiter);q

    // nodal_solver->setVerbose   (nodal_mg_verbose);
    // nodal_solver->setCGVerbose (nodal_mg_cg_verbose);
    // nodal_solver->setCGMaxIter (nodal_mg_cg_maxiter);
}
