#include <mfix_diff_F.H>
#include <mfix.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>

// For multigrid
#include <AMReX_MLMG.H>
#include <AMReX_MLTensorOp.H>
#include <AMReX_MLEBTensorOp.H>

//
// Computation of divtau
//
void
mfix::mfix_compute_divtau ( Vector< std::unique_ptr<MultiFab> >& divtau,
                            Vector< std::unique_ptr<MultiFab> >& vel   ,
                            Real time)
{
    BL_PROFILE("mfix::mfix_diffuse_explicit");

    Vector<std::unique_ptr<MultiFab> > divtau_aux(nlev);
    for (int lev = 0; lev < nlev; lev++)
    {
        divtau_aux[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost,
                                           MFInfo(), *ebfactory[lev]));
        divtau_aux[lev]->setVal(0.0);
    }

    // Swap ghost cells and apply BCs to velocity
    int extrap_dir_bcs = 0;
    mfix_set_velocity_bcs (time, vel, extrap_dir_bcs);

    // The boundary conditions need only be set once -- we do this at level 0
    int bc_lo[3], bc_hi[3];

    // Whole domain
    Box domain(geom[0].Domain());

    // Set BCs for Poisson's solver
    set_diff_bc (bc_lo, bc_hi,
                 domain.loVect(), domain.hiVect(),
                 &nghost,
                 bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
                 bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
                 bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

    //
    // First define the operator "ebtensorop"
    //
    //       (alpha * a - beta * (del dot b grad)) sol
    //
    // LPInfo                       info;
    MLEBTensorOp ebtensorop(geom, grids, dmap, LPInfo().setMaxCoarseningLevel(0),
                            amrex::GetVecOfConstPtrs(ebfactory));

    // It is essential that we set MaxOrder of the solver to 2
    // if we want to use the standard sol(i)-sol(i-1) approximation
    // for the gradient at Dirichlet boundaries.
    // The solver's default order is 3 and this uses three points for the
    // gradient at a Dirichlet boundary.
    ebtensorop.setMaxOrder(2);

    // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
    ebtensorop.setDomainBC ( {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]},
                             {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]} );

    // Return div (mu grad)) phi
    ebtensorop.setScalars(0.0, -1.0);

    // Compute the coefficients
    for (int lev = 0; lev < nlev; lev++)
    {
        average_cellcenter_to_face( GetArrOfPtrs(bcoeff[lev]), *mu_g[lev], geom[lev] );

        bcoeff[lev][0] -> FillBoundary(geom[lev].periodicity());
        bcoeff[lev][1] -> FillBoundary(geom[lev].periodicity());
        bcoeff[lev][2] -> FillBoundary(geom[lev].periodicity());

        ebtensorop.setShearViscosity  (lev, GetArrOfConstPtrs(bcoeff[lev]));
        ebtensorop.setEBShearViscosity(lev, (*mu_g[lev]));

        ebtensorop.setLevelBC ( lev, GetVecOfConstPtrs(vel)[lev] );
    }

    MLMG solver(ebtensorop);

    solver.apply(GetVecOfPtrs(divtau_aux), GetVecOfPtrs(vel));

    //
    // Apply redistribution and divide by ro_g*ep_g
    //
    for (int lev = 0; lev < nlev; lev++)
    {
        amrex::single_level_weighted_redistribute(lev, *divtau_aux[lev], *divtau[lev], *ep_g[lev], 0, 3, geom);

        // Divide by (ro_g ep_g)
        for (int n = 0; n < 3; n++)
        {
            MultiFab::Divide( *divtau[lev], *ep_g[lev], 0, n, 1, 0 );
            MultiFab::Divide( *divtau[lev], *ro_g[lev], 0, n, 1, 0 );
        }
    }
}
