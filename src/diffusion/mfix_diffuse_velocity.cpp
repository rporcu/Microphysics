#include <AMReX_MultiFabUtil.H>
#include <mfix_diffusion_op.H>
#include <mfix_eb_parms.H>

using namespace amrex;

//
// Implicit tensor solve for velocity diffusion
//
void DiffusionOp::diffuse_velocity (const Vector< MultiFab* >& vel_in,
                                    const Vector< MultiFab* >& ep_ro_in,
                                    const Vector< MultiFab* >& T_g_in,
                                    const int advect_enthalpy,
                                    Real dt,
                                    const amrex::Vector< const amrex::MultiFab* >& eb_flow_vel)
{
    BL_PROFILE("DiffusionOp::diffuse_velocity");

    int finest_level = amrcore->finestLevel();

    Vector< MultiFab* > mu_g(finest_level+1);

    Vector<BCRec> bcs_dummy; // This is just to satisfy the call to EB_interp...
    bcs_dummy.resize(3);

    auto& fluid_parms = *fluid.parameters;

    for(int lev = 0; lev <= finest_level; lev++)
    {
      mu_g[lev] = new MultiFab(ep_ro_in[lev]->boxArray(), ep_ro_in[lev]->DistributionMap(),
                               ep_ro_in[lev]->nComp(), ep_ro_in[lev]->nGrow(), MFInfo(),
                               ep_ro_in[lev]->Factory());

      mu_g[lev]->setVal(0);

      const Real mu_g0 = fluid.mu_g0;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*vel_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        Box const& bx = mfi.growntilebox(vel_in[lev]->nGrowVect());

        if (bx.ok())
        {
          Array4<Real      > const& mu_g_array = mu_g[lev]->array(mfi);
          Array4<Real const> const& T_g_array  = advect_enthalpy ?
            T_g_in[lev]->const_array(mfi) : Array4<const Real>();

          ParallelFor(bx, [mu_g_array,T_g_array,advect_enthalpy,mu_g0,fluid_parms]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            if (advect_enthalpy)
              mu_g_array(i,j,k) = fluid_parms.calc_mu_g(T_g_array(i,j,k));
            else
              mu_g_array(i,j,k) = mu_g0;
          });
        }
      }

//      EB_set_covered(*mu_g[lev], 0, mu_g[lev]->nComp(), mu_g[lev]->nGrow(), covered_val);
      EB_set_covered(*mu_g[lev], 0, mu_g[lev]->nComp(), mu_g[lev]->nGrow(), 1.e40);
    }

    // Update the coefficients of the matrix going into the solve based on the
    // current state of the simulation. Recall that the relevant matrix is
    //
    //      alpha a - beta div ( b grad )   <--->   rho - dt div ( mu grad )
    //
    // So the constants and variable coefficients are:
    //
    //      alpha: 1
    //      beta: dt
    //      a: ro
    //      b: mu

    // Set alpha and beta
    vel_matrix->setScalars(1.0, dt);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // Compute the spatially varying b coefficients (on faces) to equal the
        // apparent viscosity
        // average_cellcenter_to_face(GetArrOfPtrs(b[lev]), *mu_g[lev], geom[lev]);
        EB_interp_CellCentroid_to_FaceCentroid (*mu_g[lev], GetArrOfPtrs(b[lev]), 0, 0, 1, geom[lev], bcs_dummy);

        // This sets the coefficients
        vel_matrix->setACoeffs(lev, (*ep_ro_in[lev]));
        vel_matrix->setShearViscosity  (lev, GetArrOfConstPtrs(b[lev]), MLMG::Location::FaceCentroid);
        vel_matrix->setEBShearViscosity(lev, (*mu_g[lev]));

        if (EB::has_flow) {
            vel_matrix->setEBShearViscosityWithInflow(lev, (*mu_g[lev]), *eb_flow_vel[lev]);
        }
    }

    if(verbose > 0)
        amrex::Print() << "Diffusing velocity components all together..." << std::endl;

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // Set the right hand side to equal rho
        MultiFab::Copy((*rhs[lev]),(*vel_in[lev]), 0, 0, 3, 0);

        // Multiply rhs by rho to get momentum
        // Note that vel holds the updated velocity:
        //
        //      u_old + dt ( - u grad u + div ( mu (grad u)^T ) / rho - grad p / rho + gravity )
        //
        for (int i = 0; i < 3; i++)
           MultiFab::Multiply((*rhs[lev]), (*ep_ro_in[lev]), 0, i, 1, 0);

        // By this point we must have filled the Dirichlet values of phi stored in ghost cells
        MultiFab::Copy(*phi[lev],*vel_in[lev], 0, 0, 3, 1);
        vel_matrix->setLevelBC(lev, GetVecOfConstPtrs(phi)[lev]);

        // matrix->setEBHomogDirichlet(lev, *mu_g[lev]);
    }

    MLMG solver(*vel_matrix);
    setSolverSettings(solver);

    // This ensures that ghost cells of sol are correctly filled when returned from the solver
    solver.setFinalFillBC(true);

    solver.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        MultiFab::Copy(*vel_in[lev], *phi[lev], 0, 0, AMREX_SPACEDIM, 1);
    }

    if(verbose > 0)
        amrex::Print() << " Done diffusing all velocity components" << std::endl;

    for (int lev = 0; lev <= finest_level; lev++)
    {
      delete mu_g[lev];
    }
}
