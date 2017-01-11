#include <mfix_F.H>
#include <mfix_level.H>

//
// Do a one-component dot product of r & z using supplied components.
//
static
Real
dotxy (const MultiFab& r,
       const MultiFab& z,
       bool            local = false)
{
    const int ncomp = 1;
    const int nghost = 0;
    return MultiFab::Dot(r,0,z,0,ncomp,nghost,local);
}

static
void
sxay (MultiFab&       ss,
      const MultiFab& xx,
      Real            a,
      const MultiFab& yy,
      int             yycomp)
{
    const int ncomp  = 1;
    const int sscomp = 0;
    const int xxcomp = 0;
    MultiFab::LinComb(ss, 1.0, xx, xxcomp, a, yy, yycomp, sscomp, ncomp, 0);
}

inline
void
sxay (MultiFab&       ss,
      const MultiFab& xx,
      Real            a,
      const MultiFab& yy)
{
    sxay(ss,xx,a,yy,0);
}

int
mfix_level::solve_bicgstab (MultiFab&       sol,
                            const MultiFab& rhs,
                            const MultiFab& A_m,
                            int             sweep_type,
                            int             precond_type,
                            int             maxiter,
                            Real            eps_rel,
                            Real            eps_abs)
{
    // We're not quite ready to use this yet ... just want it all to compile
    int ret = 0, nit = 1;

    const int ncomp  = 1;
    const int nghost = sol.nGrow();

    const BoxArray& ba = sol.boxArray();
    const DistributionMapping& dm = sol.DistributionMap();

    MultiFab ph(ba, ncomp, nghost, dm);
    MultiFab sh(ba, ncomp, nghost, dm);

    MultiFab sorig(ba, ncomp, 0, dm);
    MultiFab p    (ba, ncomp, 0, dm);
    MultiFab r    (ba, ncomp, 0, dm);
    MultiFab s    (ba, ncomp, 0, dm);
    MultiFab rh   (ba, ncomp, 0, dm);
    MultiFab v    (ba, ncomp, 0, dm);
    MultiFab t    (ba, ncomp, 0, dm);

//  Lp.residual(r, rhs, sol, bc_mode);
    Array<int> slo(3);
    Array<int> shi(3);

    std::cout << "Hello from BiCGStab! " << '\n';

    std::cout << "Calc leq_residual " << '\n';

    for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
    {
      const Box& bx = (mfi.validbox()).shift(IntVect(2,2,2));

      const int* sslo = rhs[mfi].loVect();
      const int* sshi = rhs[mfi].hiVect();

      slo[0] = sslo[0]+2;
      slo[1] = sslo[1]+2;
      slo[2] = sslo[2]+2;

      shi[0] = sshi[0]+2;
      shi[1] = sshi[1]+2;
      shi[2] = sshi[2]+2;

      // Compute r = rhs - A_m*sol
      leq_residual(rhs[mfi].dataPtr(), sol[mfi].dataPtr(), A_m[mfi].dataPtr(),
                   r[mfi].dataPtr(), slo.dataPtr(),shi.dataPtr(),bx.loVect(),bx.hiVect());
    }


    std::cout << "copy sol to sorig " << '\n';
    MultiFab::Copy(sorig,sol,0,0,1,0);
    std::cout << "copy r to rh " << '\n';
    MultiFab::Copy(rh,   r,  0,0,1,0);
    std::cout << "done with copy " << '\n';

    sol.setVal(0);

    Real rnorm = r.norm0();
    const Real rnorm0   = rnorm;

    if ( verbose > 0 && ParallelDescriptor::IOProcessor() )
    {
      std::cout << "BiCGStab: Initial error (error0) =        " << rnorm0 << '\n';
    }
    Real rho_1 = 0, alpha = 0, omega = 0;

    if ( rnorm0 == 0 || rnorm0 < eps_abs )
    {
      if ( verbose > 0 && ParallelDescriptor::IOProcessor())
      {
        std::cout << "BiCGStab: niter = 0,"
                  << ", rnorm = " << rnorm
                  << ", eps_abs = " << eps_abs << std::endl;
      }
      return ret;
    }




    for (; nit <= maxiter; ++nit)
    {
      const Real rho = dotxy(rh,r);
      if ( rho == 0 )
        {
          ret = 1; break; }
      if ( nit == 1 )
      {
        MultiFab::Copy(p,r,0,0,1,0);
      }
      else
      {
        const Real beta = (rho/rho_1)*(alpha/omega);
        sxay(p, p, -omega, v);
        sxay(p, r,   beta, p);
      }

      //    if ( use_jacobi_precond )
      if ( precond_type == 0 ) // pc_type == line
      {
        ph.setVal(0);
        // Lp.jacobi_smooth(ph, p, temp_bc_mode);
      }
      else if ( precond_type = 1 ) // pc_type == diag
      {
        for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
        {
          const Box& bx = (mfi.validbox()).shift(IntVect(2,2,2));

          const int* sslo = rhs[mfi].loVect();
          const int* sshi = rhs[mfi].hiVect();

          slo[0] = sslo[0]+2;
          slo[1] = sslo[1]+2;
          slo[2] = sslo[2]+2;

          shi[0] = sshi[0]+2;
          shi[1] = sshi[1]+2;
          shi[2] = sshi[2]+2;

          leq_msolve1(bx.loVect(),bx.hiVect(),p[mfi].dataPtr(), A_m[mfi].dataPtr(),
                      ph[mfi].dataPtr(), &sweep_type);
        }
      }
      else // pc_type ==None
      {
        MultiFab::Copy(ph,p,0,0,1,0);
      }


      for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
      {
        const Box& bx = (mfi.validbox()).shift(IntVect(2,2,2));

        const int* sslo = rhs[mfi].loVect();
        const int* sshi = rhs[mfi].hiVect();

        slo[0] = sslo[0]+2;
        slo[1] = sslo[1]+2;
        slo[2] = sslo[2]+2;

        shi[0] = sshi[0]+2;
        shi[1] = sshi[1]+2;
        shi[2] = sshi[2]+2;

        leq_matvec(ph[mfi].dataPtr(), A_m[mfi].dataPtr(), v[mfi].dataPtr(),
                   slo.dataPtr(),shi.dataPtr(),bx.loVect(),bx.hiVect());
      }


      if ( Real rhTv = dotxy(rh,v) )
      {
        alpha = rho/rhTv;
      }
      else
      {
        ret = 2; break;
      }
      sxay(sol, sol,  alpha, ph);
      sxay(s,     r, -alpha,  v);

      rnorm = s.norm0();

      if ( verbose > 2 && ParallelDescriptor::IOProcessor())
      {
        std::cout << "BiCGStab: Half Iter " << nit << " rel. err. "
                  << rnorm/(rnorm0) << '\n';
      }

      if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) break;

      //    if ( use_jacobi_precond )
      if ( precond_type == 0 ) // pc_type == line
      {
        ph.setVal(0);
        // Lp.jacobi_smooth(ph, p, temp_bc_mode);
      }
      else if ( precond_type = 1 ) // pc_type == diag
      {
        for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
        {
          const Box& bx = (mfi.validbox()).shift(IntVect(2,2,2));

          const int* sslo = rhs[mfi].loVect();
          const int* sshi = rhs[mfi].hiVect();

          slo[0] = sslo[0]+2;
          slo[1] = sslo[1]+2;
          slo[2] = sslo[2]+2;

          shi[0] = sshi[0]+2;
          shi[1] = sshi[1]+2;
          shi[2] = sshi[2]+2;

          leq_msolve1(bx.loVect(),bx.hiVect(),s[mfi].dataPtr(), A_m[mfi].dataPtr(),
                      sh[mfi].dataPtr(), &sweep_type);
        }
      }
      else // pc_type ==None
      {
        MultiFab::Copy(sh,s,0,0,1,0);
      }


      for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
      {
        const Box& bx = (mfi.validbox()).shift(IntVect(2,2,2));

        const int* sslo = rhs[mfi].loVect();
        const int* sshi = rhs[mfi].hiVect();

        slo[0] = sslo[0]+2;
        slo[1] = sslo[1]+2;
        slo[2] = sslo[2]+2;

        shi[0] = sshi[0]+2;
        shi[1] = sshi[1]+2;
        shi[2] = sshi[2]+2;

        leq_matvec(sh[mfi].dataPtr(), A_m[mfi].dataPtr(), t[mfi].dataPtr(),
                   slo.dataPtr(),shi.dataPtr(),bx.loVect(),bx.hiVect());
      }

      // This is a little funky.  I want to elide one of the reductions
      // in the following two dotxy()s.  We do that by calculating the "local"
      // values and then reducing the two local values at the same time.
      Real vals[2] = { dotxy(t,t,true), dotxy(t,s,true) };

      ParallelDescriptor::ReduceRealSum(vals,2);

      if ( vals[0] )
      {
        omega = vals[1]/vals[0];
      }
      else
      {
        ret = 3; break;
      }
      sxay(sol, sol,  omega, sh);
      sxay(r,     s, -omega,  t);

      rnorm = r.norm0();

      if ( verbose > 2 && ParallelDescriptor::IOProcessor())
      {
        std::cout << "BiCGStab: Iteration " << nit << " rel. err. "
                  << rnorm/(rnorm0) << '\n';
      }

      if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) break;

      if ( omega == 0 )
      {
        ret = 4; break;
      }
      rho_1 = rho;
    }

    if ( verbose > 0 && ParallelDescriptor::IOProcessor())
    {
      std::cout << "BiCGStab: Final: Iteration " << nit << " rel. err. "
                << rnorm/(rnorm0) << '\n';
    }

    if ( ret == 0 && rnorm > eps_rel*rnorm0 && rnorm > eps_abs)
    {
      if ( ParallelDescriptor::IOProcessor())
        BoxLib::Error("BiCGStab:: failed to converge!");
      ret = 8;
    }

    if ( ( ret == 0 || ret == 8 ) && (rnorm < rnorm0) )
    {
      sol.plus(sorig, 0, 1, 0);
    }
    else
    {
      sol.setVal(0);
      sol.plus(sorig, 0, 1, 0);
    }

    return ret;
}
