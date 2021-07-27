#include <mfix_solvers.H>


using namespace amrex;
using namespace Solvers;


void DampedNewton::solve(const amrex::Vector<amrex::MultiFab*>& solution,
                         ResidueMF& R,
                         GradientMF& partial_R,
                         NormMF& norm,
                         DumpingFactor dumping_factor,
                         const amrex::Real /*abs_tol*/,
                         const amrex::Real rel_tol,
                         const int max_iterations)
{
  const int finest_level = solution.size()-1;

  // Allocate auxiliary data
  amrex::Vector<amrex::MultiFab*> residue(finest_level+1, nullptr);
  //amrex::Vector<amrex::MultiFab*> residue_old(finest_level+1, nullptr);
  amrex::Vector<amrex::MultiFab*> gradient(finest_level+1, nullptr);
  amrex::Vector<amrex::MultiFab*> update(finest_level+1, nullptr);
  amrex::Vector<amrex::MultiFab*> solution_old_old(finest_level+1, nullptr);
  amrex::Vector<amrex::MultiFab*> original_solution(finest_level+1, nullptr);

  for (int lev(0); lev <= finest_level; ++lev) {
    residue[lev] = new amrex::MultiFab(solution[lev]->boxArray(),
                                       solution[lev]->DistributionMap(),
                                       solution[lev]->nComp(),
                                       solution[lev]->nGrow(),
                                       amrex::MFInfo(), solution[lev]->Factory());
    residue[lev]->setVal(0.);


    //residue_old[lev] = new amrex::MultiFab(solution[lev]->boxArray(),
    //                                       solution[lev]->DistributionMap(),
    //                                       solution[lev]->nComp(),
    //                                       solution[lev]->nGrow(),
    //                                       amrex::MFInfo(), solution[lev]->Factory());
    //residue_old[lev]->setVal(0.);


    gradient[lev] = new amrex::MultiFab(solution[lev]->boxArray(),
                                        solution[lev]->DistributionMap(),
                                        solution[lev]->nComp(),
                                        solution[lev]->nGrow(),
                                        amrex::MFInfo(), solution[lev]->Factory());
    gradient[lev]->setVal(0.);


    update[lev] = new amrex::MultiFab(solution[lev]->boxArray(),
                                      solution[lev]->DistributionMap(),
                                      solution[lev]->nComp(), 0, amrex::MFInfo(),
                                      solution[lev]->Factory());
    update[lev]->setVal(0.);
    amrex::MultiFab::Copy(*update[lev], *solution[lev], 0, 0, 1, 0);
    amrex::EB_set_covered(*update[lev], 0, 1, 0, 0.);


    solution_old_old[lev] = new amrex::MultiFab(solution[lev]->boxArray(),
                                           solution[lev]->DistributionMap(),
                                           solution[lev]->nComp(), 0,
                                           amrex::MFInfo(), solution[lev]->Factory());
    solution_old_old[lev]->setVal(0.);
    amrex::MultiFab::Copy(*solution_old_old[lev], *solution[lev], 0, 0, 1, 0);
    amrex::EB_set_covered(*solution_old_old[lev], 0, 1, 0, 0.);


    original_solution[lev] = new amrex::MultiFab(solution[lev]->boxArray(),
                                                 solution[lev]->DistributionMap(),
                                                 solution[lev]->nComp(),
                                                 solution[lev]->nGrow(),
                                                 amrex::MFInfo(), solution[lev]->Factory());
    amrex::MultiFab::Copy(*original_solution[lev], *solution[lev], 0, 0,
                          solution[lev]->nComp(), solution[lev]->nGrow());
  }

  R(residue, solution);

  int iter(0);
  const amrex::Real update_rel_tol = rel_tol*norm(update);
  //const amrex::Real residue_rel_tol = rel_tol*norm(residue);

  for (int lev(0); lev <= finest_level; ++lev) {
    update[lev]->setVal(0.);
  }

  do {
    partial_R(gradient, solution);

    for (int lev(0); lev <= finest_level; ++lev) {
      const auto& factory = dynamic_cast<amrex::EBFArrayBoxFactory const&>(solution[lev]->Factory());
      const auto& flags = factory.getMultiEBCellFlagFab();
      const auto& volfrac = factory.getVolFrac();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
      for (amrex::MFIter mfi(*solution[lev]); mfi.isValid(); ++mfi)
      {
        amrex::Box const& bx = mfi.tilebox();

        if (bx.ok()) {
          amrex::Array4<amrex::Real      > const& update_arr      = update[lev]->array(mfi);
          amrex::Array4<amrex::Real      > const& solution_arr    = solution[lev]->array(mfi);
          amrex::Array4<amrex::Real      > const& solution_oo_arr = solution_old_old[lev]->array(mfi);
          //amrex::Array4<amrex::Real const> const& epg_arr       = ep_g[lev]->const_array(mfi);
          amrex::Array4<amrex::Real const> const& residue_arr     = residue[lev]->const_array(mfi);
          //amrex::Array4<amrex::Real      > const& residue_o_arr   = residue_old[lev]->array(mfi);
          amrex::Array4<amrex::Real const> const& gradient_arr    = gradient[lev]->const_array(mfi);

          auto const& flags_arr = flags.const_array(mfi);
          auto const& volfrac_arr = volfrac.const_array(mfi);

          amrex::ParallelFor(bx, [update_arr,residue_arr,/*residue_o_arr,*/
              solution_arr,/*epg_arr,*/gradient_arr,flags_arr,solution_oo_arr,iter,
              volfrac_arr,dumping_factor]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            if (!flags_arr(i,j,k).isCovered()) {
              //const amrex::Real epg = epg_arr(i,j,k);
              const amrex::Real vfrac = volfrac_arr(i,j,k);

              const amrex::Real sol_o = solution_arr(i,j,k);
              const amrex::Real sol_oo = solution_oo_arr(i,j,k);

              const amrex::Real res = residue_arr(i,j,k);
              //const amrex::Real res_o = residue_o_arr(i,j,k);

              const amrex::Real learning_rate = dumping_factor(1., vfrac);

              amrex::Real sol = sol_o - learning_rate*(res/gradient_arr(i,j,k));

              // If the method is diverging, try to recover convergence with
              // a midpoint step
              if (iter != 0) {
                if (std::abs(sol - sol_oo) < std::abs(sol - sol_o)) {
                  sol = (sol_o + sol_oo) / 2.;
                }
                else if (std::abs(sol_o - sol_oo) < std::abs(sol - sol_o)) {
                  sol = (sol + sol_oo) / 2.;
                }

                //residue_o_arr(i,j,k) = res;
              }

              solution_arr(i,j,k) = sol;
              solution_oo_arr(i,j,k) = sol_o;

              update_arr(i,j,k) = std::abs(sol - sol_o);
            }
          });
        }
      }
    }

    R(residue, solution);

    ++iter;

    if(iter > max_iterations) {
      // Reset the solution
      for (int lev(0); lev <= finest_level; ++lev) {
        amrex::MultiFab::Copy(*solution[lev], *original_solution[lev], 0, 0,
                              solution[lev]->nComp(), solution[lev]->nGrow());
      }

      NonConvergingExc exception(iter, norm(update), norm(residue));
      throw exception;
    }
  } while(//(norm(residue) > residue_rel_tol) ||
          (norm(update) > update_rel_tol));

  amrex::Print() << "Damped-Newton converged after iterations nb. = " << iter << "\n";
  amrex::Print() << "Damped-Newton final update norm = " << norm(update) << "\n";
  amrex::Print() << "Damped-Newton final residue norm = " << norm(residue) << "\n";

  for (int lev(0); lev <= finest_level; ++lev) {
    delete residue[lev];
    //delete residue_old[lev];
    delete gradient[lev];
    delete update[lev];
    delete solution_old_old[lev];
    delete original_solution[lev];
  }
}
