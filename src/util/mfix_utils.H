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

  const amrex::MultiFab& volfrac =  ebfactory[lev]->getVolFrac();

  amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
  amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {

    const amrex::Box& bx = mfi.tilebox();

    amrex::Array4<amrex::Real const> const& mf_arr   = mf.const_array(mfi);
    amrex::Array4<amrex::Real const> const& vfrc_arr = volfrac.const_array(mfi);

    reduce_op.eval(bx, reduce_data, [comp,mf_arr,vfrc_arr]
        AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
    {
      return mf_arr(i,j,k,comp) * vfrc_arr(i,j,k);
    });
  }

  ReduceTuple host_tuple = reduce_data.value(reduce_op);

  amrex::Real sum = amrex::get<0>(host_tuple);

  if (!local)
    amrex::ParallelDescriptor::ReduceRealSum(&sum, 1);

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

  const amrex::MultiFab& volfrac =  ebfactory[lev]->getVolFrac();

  amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
  amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {

    const amrex::Box& bx = mfi.tilebox();

    amrex::Array4<amrex::Real const> const& mf_arr   = mf.const_array(mfi);
    amrex::Array4<amrex::Real const> const& epg_arr  = ep_g.const_array(mfi);
    amrex::Array4<amrex::Real const> const& vfrc_arr = volfrac.const_array(mfi);

    reduce_op.eval(bx, reduce_data, [comp,mf_arr,vfrc_arr,epg_arr]
        AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
    {
      return mf_arr(i,j,k,comp) * vfrc_arr(i,j,k) * epg_arr(i,j,k);
    });
  }

  ReduceTuple host_tuple = reduce_data.value(reduce_op);

  amrex::Real sum = amrex::get<0>(host_tuple);

  if (!local)
    amrex::ParallelDescriptor::ReduceRealSum(&sum, 1);

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

  const amrex::MultiFab& volfrac =  ebfactory[lev]->getVolFrac();

  bool reduction_computed(false);

  amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
  amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {

    const amrex::Box bx = mfi.tilebox() & a_bx;

    if (bx.ok()) {
      reduction_computed = true;

      amrex::Array4<amrex::Real const> const& mf_arr   = mf.const_array(mfi);
      amrex::Array4<amrex::Real const> const& vfrc_arr = volfrac.const_array(mfi);

      reduce_op.eval(bx, reduce_data, [comp,mf_arr,vfrc_arr]
          AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
      {
        return mf_arr(i,j,k,comp) * vfrc_arr(i,j,k);
      });
    }
  }

  ReduceTuple host_tuple = reduce_data.value(reduce_op);

  amrex::Real sum(0.);
  
  if (reduction_computed)
    sum = amrex::get<0>(host_tuple);

  if (!local)
    amrex::ParallelDescriptor::ReduceRealSum(&sum, 1);

  return sum;
}


static amrex::Real volWgtAvg (const amrex::MultiFab& mf,
                              int comp,
                              const amrex::MultiFab& volfrac)
{
  BL_PROFILE("mfix::volWgtSum()");

  //const amrex::MultiFab& volfrac =  ebfactory[lev]->getVolFrac();

  amrex::ReduceOps<amrex::ReduceOpSum,amrex::ReduceOpSum> reduce_op;
  amrex::ReduceData<amrex::Real,amrex::Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {

    const amrex::Box& bx = mfi.tilebox();

    amrex::Array4<amrex::Real const> const& mf_arr   = mf.const_array(mfi);
    amrex::Array4<amrex::Real const> const& vfrc_arr = volfrac.const_array(mfi);

    reduce_op.eval(bx, reduce_data, [comp,mf_arr,vfrc_arr]
        AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
    {
      return {mf_arr(i,j,k,comp) * vfrc_arr(i,j,k), vfrc_arr(i,j,k)};
    });
  }

  ReduceTuple host_tuple = reduce_data.value(reduce_op);

  amrex::Real sum = amrex::get<0>(host_tuple);
  amrex::Real vol = amrex::get<1>(host_tuple);

  amrex::ParallelDescriptor::ReduceRealSum(&sum, 1);
  amrex::ParallelDescriptor::ReduceRealSum(&vol, 1);

  return sum/vol;
}


static amrex::Real volEpsWgtAvg (const amrex::MultiFab& mf,
                                 int comp,
                                 const amrex::MultiFab& ep_g,
                                 const amrex::MultiFab& volfrac)
{
  BL_PROFILE("mfix::volEpsWgtSum()");

//  const amrex::MultiFab& volfrac =  ebfactory[lev]->getVolFrac();

  amrex::ReduceOps<amrex::ReduceOpSum,amrex::ReduceOpSum> reduce_op;
  amrex::ReduceData<amrex::Real,amrex::Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {

    const amrex::Box& bx = mfi.tilebox();

    amrex::Array4<amrex::Real const> const& mf_arr   = mf.const_array(mfi);
    amrex::Array4<amrex::Real const> const& epg_arr  = ep_g.const_array(mfi);
    amrex::Array4<amrex::Real const> const& vfrc_arr = volfrac.const_array(mfi);

    reduce_op.eval(bx, reduce_data, [comp,mf_arr,vfrc_arr,epg_arr]
        AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
    {
      return {mf_arr(i,j,k,comp) * vfrc_arr(i,j,k) * epg_arr(i,j,k),
              vfrc_arr(i,j,k) * epg_arr(i,j,k)};
    });
  }

  ReduceTuple host_tuple = reduce_data.value(reduce_op);

  amrex::Real sum = amrex::get<0>(host_tuple);
  amrex::Real vol = amrex::get<1>(host_tuple);

  amrex::ParallelDescriptor::ReduceRealSum(&sum, 1);
  amrex::ParallelDescriptor::ReduceRealSum(&vol, 1);

  return sum/vol;
}


static amrex::Real volWgtAvgBox (const amrex::MultiFab & mf,
                                 int comp,
                                 const amrex::MultiFab& volfrac,
                                 const amrex::Box a_bx)
{
  BL_PROFILE("mfix::volWgtSumBox()");

//  const amrex::MultiFab& volfrac =  ebfactory[lev]->getVolFrac();

  bool reduction_computed(false);

  amrex::ReduceOps<amrex::ReduceOpSum,amrex::ReduceOpSum> reduce_op;
  amrex::ReduceData<amrex::Real,amrex::Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {

    const amrex::Box bx = mfi.tilebox() & a_bx;

    if (bx.ok()) {
      reduction_computed = true;

      amrex::Array4<amrex::Real const> const& mf_arr   = mf.const_array(mfi);
      amrex::Array4<amrex::Real const> const& vfrc_arr = volfrac.const_array(mfi);

      reduce_op.eval(bx, reduce_data, [comp,mf_arr,vfrc_arr]
          AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
      {
        return {mf_arr(i,j,k,comp) * vfrc_arr(i,j,k), vfrc_arr(i,j,k)};
      });
    }
  }

  ReduceTuple host_tuple = reduce_data.value(reduce_op);

  amrex::Real sum(0.);
  amrex::Real vol(0.);
  
  if (reduction_computed) {
    sum = amrex::get<0>(host_tuple);
    vol = amrex::get<1>(host_tuple);
  }

  amrex::ParallelDescriptor::ReduceRealSum(&sum, 1);
  amrex::ParallelDescriptor::ReduceRealSum(&vol, 1);

  return sum/vol;
}

};


#endif
