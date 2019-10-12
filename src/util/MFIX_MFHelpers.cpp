#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <MFIX_MFHelpers.H>

using namespace std;
using namespace amrex;



//
// Creates an exact copy ("shape" + value) of "mold"
//
unique_ptr< MultiFab >
MFHelpers::createFrom (MultiFab& mold)
{
    AMREX_ASSERT(mold.ok());

    unique_ptr< MultiFab > mf;

    mf.reset(new MultiFab(mold.boxArray(), mold.DistributionMap(),
                          mold.nComp(), mold.nGrow(), MFInfo(),
                          mold.Factory() ) );

    MultiFab::Copy(*mf, mold, 0, 0, mold.nComp(), mold.nGrow());

    return mf;
}


//
// Creates a copy of "mold" and initializes it to "val".
//
unique_ptr< MultiFab >
MFHelpers::createFrom ( MultiFab& mold, Real val )
{
    AMREX_ASSERT(mold.ok());

    unique_ptr< MultiFab > mf;

    mf.reset(new MultiFab(mold.boxArray(), mold.DistributionMap(),
                          mold.nComp(), mold.nGrow(), MFInfo(),
                          mold.Factory() ) );
    mf -> setVal(val);

    return mf;
}

//
// Creates a copy of "mold" with nGrow ghosts and initializes it to "val".
//
unique_ptr< MultiFab >
MFHelpers::createFrom (MultiFab& mold, Real val, int nGrow )
{
    AMREX_ASSERT(mold.ok());

    unique_ptr< MultiFab > mf;

    mf.reset(new MultiFab(mold.boxArray(), mold.DistributionMap(),
                          mold.nComp(), nGrow, MFInfo(),
                          mold.Factory() ) );
    mf -> setVal(val);

    return mf;

}


//
// Copies src to dst including all components and ghost cells
//
void
MFHelpers::copy (MultiFab& dst, MultiFab& src)
{
    AMREX_ASSERT(dst.ok());
    AMREX_ASSERT(src.ok());
    AMREX_ASSERT(dst.nComp()==src.nComp());
    AMREX_ASSERT(dst.nGrow()==src.nGrow());

    MultiFab::Copy(dst, src, 0, 0,  dst.nComp(), dst.nGrow());
}

//
// Copies src to dst including all components and nGrow ghost cells
//
void
MFHelpers::copy (MultiFab& dst, MultiFab& src, int nGrow)
{
    AMREX_ASSERT(dst.ok());
    AMREX_ASSERT(src.ok());
    AMREX_ASSERT(dst.nComp()==src.nComp());
    AMREX_ASSERT(nGrow >= 0);
    AMREX_ASSERT(dst.nGrow()>=nGrow);
    AMREX_ASSERT(src.nGrow()>=nGrow);

    MultiFab::Copy(dst, src, 0, 0,  dst.nComp(), nGrow);
}
