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
    const int nghost = 1;
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
    MultiFab::LinComb(ss, 1.0, xx, xxcomp, a, yy, yycomp, sscomp, ncomp, 1);
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
                            Real            eps_rel)
{
    // We're not quite ready to use this yet ... just want it all to compile
    int ret = 0, nit = 1;

    const int ncomp  = 1;
    const int nghost = sol.nGrow();

    const BoxArray& ba = sol.boxArray();
    const DistributionMapping& dm = sol.DistributionMap();

    MultiFab ph(ba, ncomp, nghost, dm);
    MultiFab sh(ba, ncomp, nghost, dm);

    MultiFab sorig(ba, ncomp, nghost, dm);
    MultiFab p    (ba, ncomp, nghost, dm);
    MultiFab r    (ba, ncomp, nghost, dm);
    MultiFab s    (ba, ncomp, nghost, dm);
    MultiFab rh   (ba, ncomp, nghost, dm);
    MultiFab v    (ba, ncomp, nghost, dm);
    MultiFab t    (ba, ncomp, nghost, dm);

//  Lp.residual(r, rhs, sol, bc_mode);

    const Real eps_abs=1.0E-12;

    // Unit scaling
    //-------------------------------------------------------------------
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
    {
        const Box& rbx = rhs[mfi].box();
        const Box& abx = A_m[mfi].box();

        leq_scale(rhs[mfi].dataPtr(), rbx.loVect(), rbx.hiVect(),
                  A_m[mfi].dataPtr(), abx.loVect(), abx.hiVect());
    }


    // Compute initial residual r = rhs - A*sol
    //-------------------------------------------------------------------
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
    {
      const Box& hbx = rhs[mfi].box();
      const Box& rbx =   r[mfi].box();
      const Box& abx = A_m[mfi].box();
      const Box& sbx = sol[mfi].box();

      // Compute r = rhs - A_m*sol
      leq_residual(rhs[mfi].dataPtr(), hbx.loVect(), hbx.hiVect(),
                   sol[mfi].dataPtr(), sbx.loVect(), sbx.hiVect(),
                   A_m[mfi].dataPtr(), abx.loVect(), abx.hiVect(),
                     r[mfi].dataPtr(), rbx.loVect(), rbx.hiVect());
    }

    MultiFab::Copy(sorig,sol,0,0,1,nghost);
    MultiFab::Copy(rh,   r,  0,0,1,nghost);

    //Real rnorm = r.norm0(0,1);  HACK Below gives same value a MFIX
    Real rnorm = dotxy(r,r,true);
    ParallelDescriptor::ReduceRealSum(rnorm);
    const Real rnorm0   = sqrt(rnorm);

    if ( verbose > 0 && ParallelDescriptor::IOProcessor() )
    {
      std::cout << "BiCGStab: Initial error (error0) = " << rnorm0 << '\n';
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


    // Main loop
    //-------------------------------------------------------------------
    for (; nit <= maxiter; ++nit)
    {

      Real rho = dotxy(rh,r,true);
      if ( rho == 0 )
      {
        ret = 1;
        break;
      }


      if ( nit == 1 )
      {
        MultiFab::Copy(p,r,0,0,1,nghost);
      }
      else
      {
        const Real beta = (rho/rho_1)*(alpha/omega);
        sxay(p, p, -omega, v);
        sxay(p, r,   beta, p);
      }

      //  A*ph = p
      //  v = A*Ph
      //-----------------------------------------------------------------
      if ( precond_type == 0 ) // pc_type == line
      {
        std::cout << " DNE should not be here" << '\n';
        ph.setVal(0);
      }
      else if ( precond_type == 1) // pc_type == diag
      {
        for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
        {
          const Box& pbx =   p[mfi].box();
          const Box& abx = A_m[mfi].box();
          const Box& hbx =  ph[mfi].box();

          leq_msolve1(  p[mfi].dataPtr(), pbx.loVect(), pbx.hiVect(),
                      A_m[mfi].dataPtr(), abx.loVect(), abx.hiVect(),
                       ph[mfi].dataPtr(), hbx.loVect(), hbx.hiVect());
        }
      }
      else // pc_type ==None
      {
        MultiFab::Copy(ph,p,0,0,1,nghost);
      }


      for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
      {
        const Box& pbx =  ph[mfi].box();
        const Box& abx = A_m[mfi].box();
        const Box& vbx =   v[mfi].box();
        leq_matvec(ph[mfi].dataPtr(), pbx.loVect(), pbx.hiVect(),
                  A_m[mfi].dataPtr(), abx.loVect(), abx.hiVect(),
                    v[mfi].dataPtr(), vbx.loVect(), vbx.hiVect());
      }

      Real rhTv = dotxy(rh,v,true);
      ParallelDescriptor::ReduceRealSum(rhTv);

      // Compute alpha
      //----------------------------------------------------------------
      if ( rhTv )
      {
        alpha = rho/rhTv;
      }
      else
      {
        ret = 2; break;
      }
      // Compute s
      //----------------------------------------------------------------
      sxay(s,     r, -alpha,  v);


      // A*sh = s
      // t=A*sh
      //----------------------------------------------------------------
      if ( precond_type == 0 ) // pc_type == line
      {
        ph.setVal(0);
      }
      else if ( precond_type == 1 ) // pc_type == diag
      {
        for (MFIter mfi(A_m); mfi.isValid(); ++mfi)
        {
          const Box& sbx =   s[mfi].box();
          const Box& abx = A_m[mfi].box();
          const Box& hbx =  sh[mfi].box();

         leq_msolve1(  s[mfi].dataPtr(), sbx.loVect(), sbx.hiVect(),
                     A_m[mfi].dataPtr(), abx.loVect(), abx.hiVect(),
                      sh[mfi].dataPtr(), hbx.loVect(), hbx.hiVect());

        }
      }
      else // pc_type ==None
      {
        MultiFab::Copy(sh,s,0,0,1,nghost);
      }

      for (MFIter mfi(A_m); mfi.isValid(); ++mfi)
      {
        const Box& sbx =  sh[mfi].box();
        const Box& abx = A_m[mfi].box();
        const Box& tbx =   t[mfi].box();

        leq_matvec(sh[mfi].dataPtr(), sbx.loVect(), sbx.hiVect(),
                  A_m[mfi].dataPtr(), abx.loVect(), abx.hiVect(),
                    t[mfi].dataPtr(), tbx.loVect(), tbx.hiVect());
      }

      // This is a little funky.  I want to elide one of the reductions
      // in the following two dotxy()s.  We do that by calculating the "local"
      // values and then reducing the two local values at the same time.
      Real vals[2] = { dotxy(t,t,true), dotxy(t,s,true) };
      ParallelDescriptor::ReduceRealSum(vals,2);


      // Compute omega
      // ---------------------------------------------------------------
      if ( vals[0] )
      {
        omega = vals[1]/vals[0];
      }
      else
      {
        ret = 3; break;
      }

      sxay(sol, sol,  alpha, ph);
      sxay(sol, sol,  omega, sh);

      sxay(r,     s, -omega,  t);


      //rnorm = s.norm0();
      rnorm = dotxy(r,r,true);
      ParallelDescriptor::ReduceRealSum(rnorm);
      rnorm = sqrt(rnorm);

      if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) break;

      if ( omega == 0 )
      {
        ret = 4; break;
      }
      rho_1 = rho;
    }


    if ( verbose > 0 && ParallelDescriptor::IOProcessor())
    {
      std::cout << "BiCGStab: final Rnorm " << rnorm << '\n';
      std::cout << "BiCGStab: ratio " << nit << " L-2 "<< rnorm/(rnorm0) << '\n';
    }

    // HACK: MFIX isn't enforcing this. It sets an error flag but does
    // do anything with that information.
    //-------------------------------------------------------------------
    // if ( ret == 0 && rnorm > eps_rel*rnorm0 && rnorm > eps_abs)
    // {
    //   if ( ParallelDescriptor::IOProcessor())
    //     BoxLib::Error("BiCGStab:: failed to converge!");
    //   ret = 8;
    // }

    // if ( ( ret == 0 || ret == 8 ) && (rnorm < rnorm0) )
    // {
    //   sol.plus(sorig, 0, 1, 0);
    // }
    // else
    // {
    //   sol.setVal(0);
    //   sol.plus(sorig, 0, 1, 0);
    // }

    return ret;
}
