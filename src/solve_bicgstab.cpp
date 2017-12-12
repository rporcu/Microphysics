#include <mfix_F.H>
#include <mfix_level.H>

//
// Do a one-component dot product of r & z using supplied components.
//
static
Real
dotxy (const MultiFab& r,
       const MultiFab& z,
       const Periodicity& period,
       bool            local = false)
{
    const int ncomp = 1;
    const int nghost = 0;
    BL_ASSERT(r.boxArray().ixType() == z.boxArray().ixType());

    Real val;

    // If the MultiFab is cell-centered we can use the standard dot product routine
    if (r.boxArray().ixType().cellCentered())
    {
      val = MultiFab::Dot(r,0,z,0,ncomp,nghost,local);
      ParallelDescriptor::ReduceRealSum(val);
    }

    // If the MultiFab is not cell-centered we have to create a special mask to make sure we
    //    don't double (or more) count the faces, edges or corners
    else
    {
       MultiFab tmpmf(r.boxArray(), r.DistributionMap(), ncomp, nghost);
       MultiFab::Copy(tmpmf, r, 0, 0, ncomp, nghost);

       auto mask = r.OverlapMask(period);
       MultiFab::Divide(tmpmf, *mask, 0, 0, ncomp, nghost);

       val = MultiFab::Dot(z, 0, tmpmf, 0, ncomp, nghost);
    }

    // Note that the MultiFab::Dot has already done the ParallelDescriptor::ReduceRealSum()
    //      so we don't need to do that here
    return val;
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
                            Real            eps_rel, int lev)
{
    BL_PROFILE("solve_bicgstab");
    Real strt_time = ParallelDescriptor::second();

    int bicg_verbose = 0;
    int ret = 0, nit = 1;

    const int ncomp  = 1;
    const int nghost = sol.nGrow();

    const BoxArray& ba = sol.boxArray();
    const DistributionMapping& dm = sol.DistributionMap();

    MultiFab ph(ba, dm, ncomp, nghost);
    MultiFab sh(ba, dm, ncomp, nghost);

    MultiFab sorig(ba, dm, ncomp, nghost);
    MultiFab rh   (ba, dm, ncomp, nghost);
    MultiFab p    (ba, dm, ncomp, nghost);
    MultiFab r    (ba, dm, ncomp, nghost);
    MultiFab s    (ba, dm, ncomp, nghost);
    MultiFab v    (ba, dm, ncomp, nghost);
    MultiFab t    (ba, dm, ncomp, nghost);

    // Initialize these to zero so valgrind doesn't complain -- in future we should look
    // at whether some of these can get away without ghost cells
    sorig.setVal(0.0);
       rh.setVal(0.0);
        p.setVal(0.0);
        r.setVal(0.0);
        s.setVal(0.0);
        t.setVal(0.0);
        v.setVal(0.0);
       ph.setVal(0.0);
       sh.setVal(0.0);

    const Real eps_abs=1.0E-12;

    // Unit scaling
    //-------------------------------------------------------------------
#ifdef _OPENMP
#pragma omp parallel
#endif 
    for (MFIter mfi(rhs, true); mfi.isValid(); ++mfi)
    {
        Real wt = ParallelDescriptor::second();

        const Box&  bx = mfi.tilebox();
        const Box& rbx = rhs[mfi].box();
        const Box& abx = A_m[mfi].box();

        leq_scale(bx.loVect(), bx.hiVect(), 
                  rhs[mfi].dataPtr(), rbx.loVect(), rbx.hiVect(),
                  A_m[mfi].dataPtr(), abx.loVect(), abx.hiVect());

       if (costs[lev]) {
	 const Box& tbx = mfi.tilebox(IntVect::TheZeroVector());
	 wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
	 (*costs[lev])[mfi].plus(wt, tbx);
       }
    }

    // Enforce periodicity on sol in case it hasn't been done yet
    sol.FillBoundary(geom[lev].periodicity());

    // Compute initial residual r = rhs - A*sol
    //-------------------------------------------------------------------
#ifdef _OPENMP
#pragma omp parallel
#endif 
    for (MFIter mfi(rhs, true); mfi.isValid(); ++mfi)
    {

      Real wt = ParallelDescriptor::second();

      const Box&  bx = mfi.tilebox();
      const Box& hbx = rhs[mfi].box();
      const Box& rbx =   r[mfi].box();
      const Box& abx = A_m[mfi].box();
      const Box& sbx = sol[mfi].box();

      // Compute r = rhs - A_m*sol
      leq_residual( bx.loVect(), bx.hiVect(), 
                   rhs[mfi].dataPtr(), hbx.loVect(), hbx.hiVect(),
                   sol[mfi].dataPtr(), sbx.loVect(), sbx.hiVect(),
                   A_m[mfi].dataPtr(), abx.loVect(), abx.hiVect(),
                     r[mfi].dataPtr(), rbx.loVect(), rbx.hiVect());

       if (costs[lev]) {
	 const Box& tbx = mfi.tilebox(IntVect::TheZeroVector());
	 wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
	 (*costs[lev])[mfi].plus(wt, tbx);
       }
    }

    r.FillBoundary(geom[lev].periodicity());

    MultiFab::Copy(sorig,sol,0,0,1,nghost);
    MultiFab::Copy(rh,   r,  0,0,1,nghost);

    Real rnorm = dotxy(r,r,geom[lev].periodicity(),true);
    const Real rnorm0   = sqrt(rnorm);

    if ( bicg_verbose > 0 && ParallelDescriptor::IOProcessor() )
    {
      std::cout << "BiCGStab: Initial error (error0) = " << rnorm0 << '\n';
    }
    Real rho_1 = 0, alpha = 0, omega = 0;

    if ( rnorm0 == 0 || rnorm0 < eps_abs )
    {
      if ( bicg_verbose > 0 && ParallelDescriptor::IOProcessor())
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

      Real rho = dotxy(rh,r,geom[lev].periodicity(),true);
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

      // Enforce periodicity on sol in case it hasn't been done yet
      p.FillBoundary(geom[lev].periodicity());

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
#ifdef _OPENMP
#pragma omp parallel
#endif 
        for (MFIter mfi(p, true); mfi.isValid(); ++mfi)
        {

	  Real wt = ParallelDescriptor::second();

          const Box&  bx = mfi.tilebox();
          const Box& pbx =   p[mfi].box();
          const Box& abx = A_m[mfi].box();
          const Box& hbx =  ph[mfi].box();

          leq_msolve1( bx.loVect(), bx.hiVect() ,
                        p[mfi].dataPtr(), pbx.loVect(), pbx.hiVect(),
                      A_m[mfi].dataPtr(), abx.loVect(), abx.hiVect(),
                       ph[mfi].dataPtr(), hbx.loVect(), hbx.hiVect());

	  if (costs[lev]) {
	    const Box& tbx = mfi.tilebox(IntVect::TheZeroVector());
	    wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
	    (*costs[lev])[mfi].plus(wt, tbx);
	  }
        }
      }
      else // pc_type ==None
      {
        MultiFab::Copy(ph,p,0,0,1,nghost);
      }

      ph.FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel
#endif 
      for (MFIter mfi(ph, true); mfi.isValid(); ++mfi)
      {

        Real wt = ParallelDescriptor::second();

        const Box&  bx = mfi.tilebox();
        const Box& pbx =  ph[mfi].box();
        const Box& abx = A_m[mfi].box();
        const Box& vbx =   v[mfi].box();
        leq_matvec( bx.loVect(), bx.hiVect(), 
		    ph[mfi].dataPtr(), pbx.loVect(), pbx.hiVect(),
		    A_m[mfi].dataPtr(), abx.loVect(), abx.hiVect(),
                    v[mfi].dataPtr(), vbx.loVect(), vbx.hiVect());

	if (costs[lev]) {
	  const Box& tbx = mfi.tilebox(IntVect::TheZeroVector());
	  wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
	  (*costs[lev])[mfi].plus(wt, tbx);
	}
      }

      // Enforce periodicity on sol in case it hasn't been done yet
      rh.FillBoundary(geom[lev].periodicity());
       v.FillBoundary(geom[lev].periodicity());

      Real rhTv = dotxy(rh,v,geom[lev].periodicity(),true);

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

      // HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK
      s.FillBoundary(geom[lev].periodicity());

      // A*sh = s
      // t=A*sh
      //----------------------------------------------------------------
      if ( precond_type == 0 ) // pc_type == line
      {
        ph.setVal(0);
      }
      else if ( precond_type == 1 ) // pc_type == diag
      {
#ifdef _OPENMP
#pragma omp parallel
#endif 
        for (MFIter mfi(A_m, true); mfi.isValid(); ++mfi)
        {

	  Real wt = ParallelDescriptor::second();

          const Box&  bx = mfi.tilebox();
          const Box& sbx =   s[mfi].box();
          const Box& abx = A_m[mfi].box();
          const Box& hbx =  sh[mfi].box();

          leq_msolve1( bx.loVect(), bx.hiVect() ,
		       s[mfi].dataPtr(), sbx.loVect(), sbx.hiVect(),
		       A_m[mfi].dataPtr(), abx.loVect(), abx.hiVect(),
		       sh[mfi].dataPtr(), hbx.loVect(), hbx.hiVect());

	  if (costs[lev]) {
	    const Box& tbx = mfi.tilebox(IntVect::TheZeroVector());
	    wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
	    (*costs[lev])[mfi].plus(wt, tbx);
	  }
        }
        sh.FillBoundary(geom[lev].periodicity());
      }
      else // pc_type ==None
      {
        MultiFab::Copy(sh,s,0,0,1,nghost);
      }

#ifdef _OPENMP
#pragma omp parallel
#endif 
      for (MFIter mfi(A_m, true); mfi.isValid(); ++mfi)
      {

        Real wt = ParallelDescriptor::second();

        const Box&  bx = mfi.tilebox();
        const Box& sbx =  sh[mfi].box();
        const Box& abx = A_m[mfi].box();
        const Box& tbx =   t[mfi].box();

        leq_matvec( bx.loVect(), bx.hiVect(), 
		    sh[mfi].dataPtr(), sbx.loVect(), sbx.hiVect(),
		    A_m[mfi].dataPtr(), abx.loVect(), abx.hiVect(),
                    t[mfi].dataPtr(), tbx.loVect(), tbx.hiVect());

	if (costs[lev]) {
	  const Box& tbx = mfi.tilebox(IntVect::TheZeroVector());
	  wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
	  (*costs[lev])[mfi].plus(wt, tbx);
	}
      }

      // Enforce periodicity on sol in case it hasn't been done yet
      sh.FillBoundary(geom[lev].periodicity());

      // This is a little funky.  I want to elide one of the reductions
      // in the following two dotxy()s.  We do that by calculating the "local"
      // values and then reducing the two local values at the same time.
      Real vals[2] = { dotxy(t,t,geom[lev].periodicity(),true),
                       dotxy(t,s,geom[lev].periodicity(),true) };

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

      sol.FillBoundary(geom[lev].periodicity());

      sxay(r,     s, -omega,  t);

      rnorm = dotxy(r,r,geom[lev].periodicity(),true);
      rnorm = sqrt(rnorm);

      if ( bicg_verbose > 0 && ParallelDescriptor::IOProcessor())
        std::cout << "BiCGStab:       Rnorm " << nit << " L-2 "<< rnorm/(rnorm0) << '\n';

      if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) break;

      if ( omega == 0 )
      {
        ret = 4; break;
      }
      rho_1 = rho;
    }

    if ( bicg_verbose > 0 && ParallelDescriptor::IOProcessor())
    {
      std::cout << "BiCGStab: final Rnorm " << rnorm << '\n';
      std::cout << "BiCGStab: ratio " << nit << " L-2 "<< rnorm/(rnorm0) << '\n';
    }

    Real end_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(end_time, ParallelDescriptor::IOProcessorNumber());
    // amrex::Print() << "Time spent in bicgsolve " << end_time << std::endl;

    return ret;
}
