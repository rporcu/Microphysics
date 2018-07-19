#include <AMReX_ParmParse.H>

#include <mfix_proj_F.H>
#include <mfix_level.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>

// For multigrid
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>

//
// Explicit diffusion
// 
void
mfix_level::mfix_compute_divtau ( int lev,
			          MultiFab& divtau, 
			          Vector< std::unique_ptr<MultiFab> >& vel) 
{
    BL_PROFILE("mfix_level::mfix_compute_divtau");
    Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*vel[lev],true); mfi.isValid(); ++mfi)
    {
	// Tilebox
	Box bx = mfi.tilebox ();

	compute_divtau (
	    BL_TO_FORTRAN_BOX(bx),  
	    BL_TO_FORTRAN_ANYD(divtau[mfi]),
	    BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
            (*mu_g[lev])[mfi].dataPtr(),
            (*lambda_g[lev])[mfi].dataPtr(),
            BL_TO_FORTRAN_ANYD((*rop_g[lev])[mfi]),
            domain.loVect (), domain.hiVect (),
	    bc_ilo.dataPtr(), bc_ihi.dataPtr(),
	    bc_jlo.dataPtr(), bc_jhi.dataPtr(),
	    bc_klo.dataPtr(), bc_khi.dataPtr(),
	    geom[lev].CellSize(), &nghost);
    }
}

//
// Implicit diffusion
// 
void
mfix_level::mfix_diffuse_velocity ( int lev, MultiFab& vel, amrex::Real dt)

{
    BL_PROFILE("mfix_level::mfix_diffuse_velocity");

    // Whole domain
    Box domain(geom[lev].Domain());
    
#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*vel_g[lev],true); mfi.isValid(); ++mfi) {
	
	// Tilebox
	Box bx = mfi.tilebox();

	compute_intermediate_velocity ( BL_TO_FORTRAN_BOX(bx),  
					BL_TO_FORTRAN_ANYD((*vel_g[lev])[mfi]),
					BL_TO_FORTRAN_ANYD((*f_gds[lev])[mfi]),
					BL_TO_FORTRAN_ANYD((*rop_g[lev])[mfi]),
					&dt );
    }
}
