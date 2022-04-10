#ifndef MFIX_UTIL_H
#define MFIX_UTIL_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Box.H>
#include <AMReX_Vector.H>
#include <AMReX_EBFabFactory.H>


struct Utils {

  static amrex::Real volWgtSum (int lev,
                                const amrex::MultiFab & mf,
                                int comp,
                                amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>& ebfactory,
                                bool local=false)
  {
    BL_PROFILE("mfix::volWgtSum()");

    const amrex::MultiFab* volfrac =  &(ebfactory[lev]->getVolFrac());

    amrex::Real sum = amrex::ReduceSum(mf, *volfrac, 0,
        [comp] AMREX_GPU_HOST_DEVICE (amrex::Box const & bx,
                                      amrex::Array4<const amrex::Real> const & rho,
                                      amrex::Array4<const amrex::Real> const & vfrc)
        {
          amrex::Real dm = 0.0;

          amrex::Loop(bx, [rho,vfrc,comp,&dm] (int i, int j, int k) noexcept
              { dm += rho(i,j,k,comp) * vfrc(i,j,k); });

          return dm;
        });

    if (!local)
        amrex::ParallelDescriptor::ReduceRealSum(sum);

    return sum;
  }


  static amrex::Real volEpsWgtSum (int lev,
                                   const amrex::MultiFab & mf,
                                   const amrex::MultiFab & ep_g,
                                   int comp,
                                   amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>& ebfactory,
                                   bool local=false)
  {
    BL_PROFILE("mfix::volEpsWgtSum()");

    const amrex::MultiFab* volfrac =  &(ebfactory[lev]->getVolFrac());

    amrex::Real sum = amrex::ReduceSum(mf, *volfrac, ep_g, 0,
        [comp] AMREX_GPU_HOST_DEVICE (amrex::Box const & bx,
                                      amrex::Array4<const amrex::Real> const & rho,
                                      amrex::Array4<const amrex::Real> const & vfrc,
                                      amrex::Array4<const amrex::Real> const & ep)
        {
          amrex::Real dm = 0.0;

          amrex::Loop(bx, [rho,vfrc,ep,comp,&dm] (int i, int j, int k) noexcept
              { dm += rho(i,j,k,comp) * vfrc(i,j,k) * ep(i,j,k); });

          return dm;
        });

    if (!local)
        amrex::ParallelDescriptor::ReduceRealSum(sum);

    return sum;
  }


  static amrex::Real volWgtSumBox (int lev,
                                   const amrex::MultiFab & mf,
                                   int comp,
                                   amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>& ebfactory,
                                   const amrex::Box a_bx,
                                   bool local=false)
  {
    BL_PROFILE("mfix::volWgtSumBox()");

    const amrex::MultiFab* volfrac =  &(ebfactory[lev]->getVolFrac());

    amrex::Real sum = amrex::ReduceSum(mf, *volfrac, 0, [comp, a_bx]
      AMREX_GPU_HOST_DEVICE (amrex::Box const & bx,
                             amrex::Array4<const amrex::Real> const & rho,
                             amrex::Array4<const amrex::Real> const & vfrc)
        {

          // We want the intersection of this box (bx) and the box
          // provided in the function call.
          const amrex::Box insect_bx = bx&a_bx;

          amrex::Real dm = 0.0;
          amrex::Loop(insect_bx, [rho,vfrc,comp,&dm] (int i, int j, int k) noexcept
              { dm += rho(i,j,k,comp) * vfrc(i,j,k); });

          return dm;
        });

    if (!local)
        amrex::ParallelDescriptor::ReduceRealSum(sum);

    return sum;
  }

};


#endif
