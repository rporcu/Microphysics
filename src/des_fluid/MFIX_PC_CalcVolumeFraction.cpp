#include <AMReX.H>
#include "AMReX_Particles.H"
#include "AMReX_RealVect.H"
#include <iostream>
#include <MFIXParticleContainer.H>
#include <AMReX_LoadBalanceKD.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EB_F.H>

#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBSupport.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>

#include <math.h>

#include "mfix_F.H"
#include "mfix_des_F.H"
#include "mfix_eb_F.H"
#include "mfix_util_F.H"
#include "mfix_des_K.H"
#include "MFIX_DEM_Parms.H"

using namespace amrex;
using namespace std;

void MFIXParticleContainer::
CalcVolumeFraction(const Vector<std::unique_ptr<MultiFab>> & mf_to_be_filled,
                   const Vector<std::unique_ptr<EBFArrayBoxFactory>> & ebfactory,
                   const Vector<std::unique_ptr<IArrayBox>> & bc_ilo,
                   const Vector<std::unique_ptr<IArrayBox>> & bc_ihi,
                   const Vector<std::unique_ptr<IArrayBox>> & bc_jlo,
                   const Vector<std::unique_ptr<IArrayBox>> & bc_jhi,
                   const Vector<std::unique_ptr<IArrayBox>> & bc_klo,
                   const Vector<std::unique_ptr<IArrayBox>> & bc_khi,
                   const BcList & bc_list,
                   int nghost )
{
    // NOTE: ebfactory HAS to be the PARTICLE EB factory!
    int fortran_volume_comp = 5;
    PICDeposition(mf_to_be_filled, ebfactory,
                  bc_ilo, bc_ihi, bc_jlo, bc_jhi, bc_klo,bc_khi,
                  bc_list, fortran_volume_comp,nghost);

    for (int lev = 0; lev < nlev; lev++)
    {
        // Now define this mf = (1 - particle_vol)
        mf_to_be_filled[lev]->mult(-1.0,mf_to_be_filled[lev]->nGrow());
        mf_to_be_filled[lev]->plus( 1.0,mf_to_be_filled[lev]->nGrow());

        // We set ep_g to 1 rather than 0 in covered cells so that when we divide by ep_g
        //    following the projection we don't have to protect against divide by 0.
        EB_set_covered(*mf_to_be_filled[lev],1.0);

        // Impose a lower bound on volume fraction
        CapSolidsVolFrac(*mf_to_be_filled[lev]);
    }

    // HACK -- we really should average down (ep_g * volfrac) not ep_g.
    for (int lev = nlev - 1; lev > 0; lev --)
    {
        amrex::EB_average_down(* mf_to_be_filled[lev], * mf_to_be_filled[lev - 1],
                               0, 1, m_gdb->refRatio(lev - 1));
    }
}

void MFIXParticleContainer::CapSolidsVolFrac(amrex::MultiFab& mf_to_be_filled)
{
    for (MFIter mfi(mf_to_be_filled); mfi.isValid(); ++mfi) {
       //const Box& sbx = mf_to_be_filled[mfi].box(); // UNUSED_VARIABLE

//       const Real max_pack = 0.42;
       //const Real max_pack = 0.21; // UNUSED_VARIABLE

       //Array4<Real> const& ep_g = mf_to_be_filled.array(mfi); // UNUSED_VARIABLE

// These lines are commented because this code represents a functionality which
// may be added in future developments
//       AMREX_HOST_DEVICE_FOR_3D(sbx, i, j, k,
//       {
//         ep_g(i,j,k) = std::max(max_pack, ep_g(i,j,k));
//       });
    }
}


void MFIXParticleContainer::
PICDeposition(const amrex::Vector< std::unique_ptr<MultiFab> >& mf_to_be_filled,
              const amrex::Vector< std::unique_ptr<EBFArrayBoxFactory>  >& ebfactory,
              const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ilo,
              const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ihi,
              const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jlo,
              const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jhi,
              const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_klo,
              const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_khi,
              const BcList & bc_list,
              int fortran_particle_comp, int nghost )
{
    BL_PROFILE("MFIXParticleContainer::PICDeposition()");

    int ncomp = 1;

    MultiFab* mf_pointer[nlev];

    // Start the timers ...
    const Real strttime = ParallelDescriptor::second();

    if (nlev > 2)
      amrex::Abort("For right now MFIXParticleContainer::PICDeposition can only handle up to 2 levels");

    for (int lev = 0; lev < nlev; lev++)
    {
       if (lev == 0 && OnSameGrids(lev, *mf_to_be_filled[lev])) {
          // If we are already working with the internal mf defined on the
          // particle_box_array, then we just work with this.
          mf_pointer[lev] = mf_to_be_filled[lev].get();

       } else if (lev == 0 && !OnSameGrids(lev, *mf_to_be_filled[lev]))  {
          // If mf_to_be_filled is not defined on the particle_box_array, then we need
          // to make a temporary here and copy into mf_to_be_filled at the end.
          mf_pointer[lev] = new MultiFab(ParticleBoxArray(lev), ParticleDistributionMap(lev),
                                         ncomp, mf_to_be_filled[lev]->nGrow());

       } else {
          // If lev > 0 we make a temporary at the coarse resolution
          BoxArray ba_crse(amrex::coarsen(ParticleBoxArray(lev),this->m_gdb->refRatio(0)));
          mf_pointer[lev] = new MultiFab(ba_crse, ParticleDistributionMap(lev),
                                         ncomp, 1);
       }

       // We must have ghost cells for each FAB so that a particle in one grid can spread
       // its effect to an adjacent grid by first putting the value into ghost cells of its
       // own grid.  The mf->sumBoundary call then adds the value from one grid's ghost cell
       // to another grid's valid region.
       if (mf_pointer[lev]->nGrow() < 1)
          amrex::Error("Must have at least one ghost cell when in CalcVolumeFraction");

       mf_pointer[lev]->setVal(0.0,0,1,mf_pointer[lev]->nGrow());
    }

    // We always use the coarse dx
    const Geometry& gm          = Geom(0);
    const auto      plo         = gm.ProbLoArray();
    const auto      dx          = gm.CellSizeArray();
    const auto      dxi         = gm.InvCellSizeArray();

    using ParConstIter = ParConstIter<realData::count,intData::count,0,0>;

    const FabArray<EBCellFlagFab>* flags;

    for (int lev = 0; lev < nlev; lev++)
    {
        Vector<int> ngrow = {1,1,1};
        std::unique_ptr<EBFArrayBoxFactory> crse_factory;

        if (lev == 0) {
            flags   = &(ebfactory[lev]->getMultiEBCellFlagFab());
        } else {
            // std::unique_ptr<EBFArrayBoxFactory>
            crse_factory = makeEBFabFactory(
                  gm, mf_pointer[lev]->boxArray(), mf_pointer[lev]->DistributionMap(),
                  ngrow, EBSupport::volume);
            flags   = &(crse_factory->getMultiEBCellFlagFab());
        }

        const MultiFab* volfrac = (lev == 0) ? &(ebfactory[lev]->getVolFrac()) : &(crse_factory->getVolFrac());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       {
        FArrayBox local_vol;
        for (ParConstIter pti(*this, lev); pti.isValid(); ++pti) {
            const auto& particles = pti.GetArrayOfStructs();
            const ParticleType* pstruct = particles().dataPtr();

            const long nrp = pti.numParticles();
            FArrayBox& fab = (*mf_pointer[lev])[pti];
#ifdef _OPENMP
            Box tile_box = pti.tilebox();
            tile_box.grow(1);
            local_vol.resize(tile_box,ncomp);
            local_vol = 0.0;
#else
            //const Box& box = fab.box(); // UNUSED_VARIABLE
#endif

            const Box& bx  = pti.tilebox(); // I need a box without ghosts

            if ((*flags)[pti].getType(bx) != FabType::covered )
            {
                auto volarr = fab.array();
                auto flagsarr = (*flags)[pti].array();
                auto vfrac = (*volfrac)[pti].array();

                AMREX_FOR_1D ( nrp, ip,
                {
                    const ParticleType& p = pstruct[ip];

                    amrex::Real reg_cell_vol = dx[0]*dx[1]*dx[2];

                    amrex::Real x = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
                    amrex::Real y = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
                    amrex::Real z = (p.pos(2) - plo[2]) * dxi[2] + 0.5;

                    int i = std::floor(x);
                    int j = std::floor(y);
                    int k = std::floor(z);

                    amrex::Real wx_hi = x - i;
                    amrex::Real wy_hi = y - j;
                    amrex::Real wz_hi = z - k;

                    amrex::Real wx_lo = 1.0 - wx_hi;
                    amrex::Real wy_lo = 1.0 - wy_hi;
                    amrex::Real wz_lo = 1.0 - wz_hi;

                    amrex::Real weights[2][2][2];

                    weights[0][0][0] = wx_lo * wy_lo * wz_lo;
                    weights[0][0][1] = wx_lo * wy_lo * wz_hi;
                    weights[0][1][0] = wx_lo * wy_hi * wz_lo;
                    weights[0][1][1] = wx_lo * wy_hi * wz_hi;
                    weights[1][0][0] = wx_hi * wy_lo * wz_lo;
                    weights[1][0][1] = wx_hi * wy_lo * wz_hi;
                    weights[1][1][0] = wx_hi * wy_hi * wz_lo;
                    weights[1][1][1] = wx_hi * wy_hi * wz_hi;

                    amrex::Real total_weight = 0.0;
                    for (int ii = 0; ii <= 1; ++ii)
                        for (int jj = 0; jj <= 1; ++jj)
                            for (int kk = 0; kk <= 1; ++kk)
                                total_weight += weights[ii][jj][kk] * vfrac(i-1+ii,j-1+jj,k-1+kk);

                    for (int ii = 0; ii <= 1; ++ii)
                        for (int jj = 0; jj <= 1; ++jj)
                            for (int kk = 0; kk <= 1; ++kk)
                                weights[ii][jj][kk] /= total_weight;

                    amrex::Real pvol = p.rdata(realData::volume) / reg_cell_vol;

                    for (int ii = -1; ii <= 0; ++ii) {
                        for (int jj = -1; jj <= 0; ++jj) {
                            for (int kk = -1; kk <= 0; ++kk) {
                                if (flagsarr(i+ii,j+jj,k+kk).isCovered())
                                    continue;
                                amrex::Gpu::Atomic::Add(&volarr(i+ii,j+jj,k+kk),
                                                        weights[ii+1][jj+1][kk+1]*pvol);
                            }
                        }
                    }
                });

                Gpu::streamSynchronize();
            }

#ifdef _OPENMP
            fab.atomicAdd(local_vol, tile_box, tile_box, 0, 0, ncomp);
#endif
        }
       }

       // Move any field deposited outside the domain back into the domain
       // when BC is pressure inlet and mass inflow.
       Box domain(Geom(lev).Domain());

       const int minf = bc_list.get_minf();
       const int pinf = bc_list.get_pinf();

       for (MFIter mfi(*mf_pointer[lev]); mfi.isValid(); ++mfi) {

         const Box& sbx = (*mf_pointer[lev])[mfi].box();

         // Extract the lower and upper boundaries of Box sbx and Domain
         const IntVect sbx_lo(sbx.loVect()), sbx_hi(sbx.hiVect());
         const amrex::Dim3& dom_lo = amrex::lbound(domain);
         const amrex::Dim3& dom_hi = amrex::ubound(domain);

         // Create a 2D Box collapsing sbx on x-direction
         IntVect sbx_yz_hi(sbx.hiVect());
         sbx_yz_hi[0] = sbx_lo[0];
         const Box sbx_yz(sbx_lo, sbx_yz_hi);

         // Create a 2D Box collapsing sbx on y-direction
         IntVect sbx_xz_hi(sbx.hiVect());
         sbx_xz_hi[1] = sbx_lo[1];
         const Box sbx_xz(sbx_lo, sbx_xz_hi);

         // Create a 2D Box collapsing sbx on z-direction
         IntVect sbx_xy_hi(sbx.hiVect());
         sbx_xy_hi[2] = sbx_lo[2];
         const Box sbx_xy(sbx_lo, sbx_xy_hi);

         Array4<int> const& bc_ilo_type = bc_ilo[lev]->array();
         Array4<int> const& bc_ihi_type = bc_ihi[lev]->array();
         Array4<int> const& bc_jlo_type = bc_jlo[lev]->array();
         Array4<int> const& bc_jhi_type = bc_jhi[lev]->array();
         Array4<int> const& bc_klo_type = bc_klo[lev]->array();
         Array4<int> const& bc_khi_type = bc_khi[lev]->array();

         Array4<Real> const& vol = mf_pointer[lev]->array(mfi);

         if(sbx_lo[0] < dom_lo.x)
         {
           const int ilo = dom_lo.x;
           AMREX_HOST_DEVICE_FOR_3D(sbx_yz, i, j, k,
           {
             if(bc_ilo_type(dom_lo.x-1,j,k,0) == pinf or
                bc_ilo_type(dom_lo.x-1,j,k,0) == minf)
             {
               vol(ilo,j,k) += vol(ilo-1,j,k);
               vol(ilo-1,j,k) = 0;
             }
           });
         }

         if(sbx_hi[0] > dom_hi.x)
         {
           const int ihi = dom_hi.x;
           AMREX_HOST_DEVICE_FOR_3D(sbx_yz, i, j, k,
           {
             if(bc_ihi_type(dom_hi.x+1,j,k,0) == pinf or
                bc_ihi_type(dom_hi.x+1,j,k,0) == minf)
             {
               vol(ihi,j,k) += vol(ihi+1,j,k);
               vol(ihi+1,j,k) = 0;
             }
           });
         }

         if(sbx_lo[1] < dom_lo.y)
         {
           const int jlo = dom_lo.y;
           AMREX_HOST_DEVICE_FOR_3D(sbx_xz, i, j, k,
           {
             if(bc_jlo_type(i,dom_lo.y-1,k,0) == pinf or
                bc_jlo_type(i,dom_lo.y-1,k,0) == minf)
             {
               vol(i,jlo,k) += vol(i,jlo-1,k);
               vol(i,jlo-1,k) = 0;
             }
           });
         }

         if(sbx_hi[1] > dom_hi.y)
         {
           const int jhi = dom_hi.y;
           AMREX_HOST_DEVICE_FOR_3D(sbx_xz, i, j, k,
           {
             if(bc_jhi_type(i,dom_hi.y+1,k,0) == pinf or
                bc_jhi_type(i,dom_hi.y+1,k,0) == minf)
             {
               vol(i,jhi,k) += vol(i,jhi+1,k);
               vol(i,jhi+1,k) = 0;
             }
           });
         }

         if(sbx_lo[2] < dom_lo.z)
         {
           const int klo = dom_lo.z;
           AMREX_HOST_DEVICE_FOR_3D(sbx_xy, i, j, k,
           {
             if(bc_klo_type(i,j,dom_lo.z-1,0) == pinf or
                bc_klo_type(i,j,dom_lo.z-1,0) == minf)
             {
               vol(i,j,klo) += vol(i,j,klo-1);
               vol(i,j,klo-1) = 0;
             }
           });
         }

         if(sbx_hi[2] > dom_hi.z)
         {
           const int khi = dom_hi.z;
           AMREX_HOST_DEVICE_FOR_3D(sbx_xy, i, j, k,
           {
             if(bc_khi_type(i,j,dom_hi.z+1,0) == pinf or
                bc_khi_type(i,j,dom_hi.z+1,0) == minf)
             {
               vol(i,j,khi) += vol(i,j,khi+1);
               vol(i,j,khi+1) = 0;
             }
           });
         }

         Gpu::streamSynchronize();
       }
    }

    int  src_nghost = 1;
    int dest_nghost = 0;
    for (int lev = 1; lev < nlev; lev++)
        mf_pointer[0]->copy(*mf_pointer[lev],0,0,ncomp,src_nghost,dest_nghost,gm.periodicity(),FabArrayBase::ADD);

    mf_pointer[0]->SumBoundary(gm.periodicity());

    if (nlev > 1)
    {
        IntVect ref_ratio(this->m_gdb->refRatio(0));

        // Now interpolate from the coarse grid to define the fine grid ep-g
        Interpolater* mapper = &cell_cons_interp;
        int lo_bc[] = {BCType::foextrap, BCType::foextrap, BCType::foextrap};
        int hi_bc[] = {BCType::foextrap, BCType::foextrap, BCType::foextrap};
        Vector<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

        BndryFuncArray bfunc(phifill);

        Real time = 0.0;
        for (int lev = 1; lev < nlev; lev++)
        {
            PhysBCFunct<BndryFuncArray> cphysbc(Geom(lev-1), bcs, bfunc);
            PhysBCFunct<BndryFuncArray> fphysbc(Geom(lev  ), bcs, bfunc);
            mf_to_be_filled[lev]->setVal(0.0);
            amrex::InterpFromCoarseLevel(*mf_to_be_filled[lev], time, *mf_pointer[lev-1],
                                         0, 0, 1, Geom(lev-1), Geom(lev),
                                         cphysbc, 0, fphysbc, 0,
                                         ref_ratio, mapper,
                                         bcs, 0);
        }
    }

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled. I believe that we don't
    // need any information in ghost cells so we don't copy those.

    if (mf_pointer[0] != mf_to_be_filled[0].get())
       mf_to_be_filled[0]->copy(*mf_pointer[0],0,0,ncomp);

    for (int lev = 0; lev < nlev; lev++)
       if (mf_pointer[lev] != mf_to_be_filled[lev].get())
          delete mf_pointer[lev];

    if (m_verbose > 1) {
      Real stoptime = ParallelDescriptor::second() - strttime;

      ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

      amrex::Print() << "MFIXParticleContainer::PICDeposition time: " << stoptime << '\n';
    }
}
