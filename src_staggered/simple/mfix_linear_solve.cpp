#include <mfix.H>
#include <mfix_F.H>

void
mfix::mfix_solve_linear_equation(int eq_id,int lev,MultiFab& sol, MultiFab& matrix, MultiFab& rhs)
{
    int sweep_type, precond_type, max_it;
    Real tol;

    BL_PROFILE("mfix::mfix_solve_linear_equation()");
    get_solver_params (&eq_id,&sweep_type,&precond_type,&max_it,&tol);

    bool debug_solve = false;
    if (debug_solve)
    {
       BoxArray ba(geom[lev].Domain()); 
       BoxArray this_ba = ba;

       // In this case we are doing pressure solve so use the cell-centered version
       if (eq_id == 1) {

       // In this case we are doing x-vel solve so use x-faces
       } else if (eq_id == 2) {
         this_ba.surroundingNodes(0);

       // In this case we are doing y-vel solve so use y-faces
       } else if (eq_id == 3) {
         this_ba.surroundingNodes(1);

       // In this case we are doing z-vel solve so use z-faces
       } else if (eq_id == 4) {
         this_ba.surroundingNodes(2);
       } 

       DistributionMapping dm(this_ba);

       int ng = sol.nGrow();
       int nc = sol.nComp();

       std::unique_ptr<MultiFab> sol_new(new MultiFab(this_ba,dm,nc,ng));
       sol_new->copy(sol,0,0,nc,ng,ng);
       sol_new->FillBoundary(geom[lev].periodicity());

       ng = rhs.nGrow();
       nc = rhs.nComp();
       std::unique_ptr<MultiFab> rhs_new(new MultiFab(this_ba,dm,nc,ng));
       rhs_new->copy(rhs,0,0,nc,ng,ng);
       rhs_new->FillBoundary(geom[lev].periodicity());

       ng = matrix.nGrow();
       nc = matrix.nComp();
       std::unique_ptr<MultiFab> matrix_new(new MultiFab(this_ba,dm,nc,ng));
       matrix_new->copy(matrix,0,0,nc,ng,ng);
       matrix_new->FillBoundary(geom[lev].periodicity());

       solve_bicgstab(*sol_new, *rhs_new, *matrix_new, sweep_type, precond_type, max_it, tol, lev);

       ng = sol.nGrow();
       sol.copy(*sol_new,0,0,1,ng,ng);
       sol.FillBoundary(geom[lev].periodicity());

    }
    else
    {
       solve_bicgstab(sol, rhs, matrix, sweep_type, precond_type, max_it, tol, lev);
    }
}
