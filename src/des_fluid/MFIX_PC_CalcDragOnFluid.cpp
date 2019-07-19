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
CalcDragOnFluid(const Vector<std::unique_ptr<MultiFab>> & drag_mf,
                const Vector<std::unique_ptr<EBFArrayBoxFactory>> & ebfactory,
                const Vector<std::unique_ptr<IArrayBox>> & bc_ilo,
                const Vector<std::unique_ptr<IArrayBox>> & bc_ihi,
                const Vector<std::unique_ptr<IArrayBox>> & bc_jlo,
                const Vector<std::unique_ptr<IArrayBox>> & bc_jhi,
                const Vector<std::unique_ptr<IArrayBox>> & bc_klo,
                const Vector<std::unique_ptr<IArrayBox>> & bc_khi,
                int nghost )
{
    int fortran_beta_comp = 15;
    int fortran_vel_comp  =  9;
    PICMultiDeposition(drag_mf, ebfactory,
                       bc_ilo, bc_ihi, bc_jlo, bc_jhi, bc_klo, bc_khi,
                       fortran_beta_comp, fortran_vel_comp, nghost);
}

void MFIXParticleContainer::
PICMultiDeposition(const amrex::Vector< std::unique_ptr<MultiFab> >& drag_mf,
                   const amrex::Vector< std::unique_ptr<EBFArrayBoxFactory>  >& ebfactory,
                   const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ilo,
                   const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ihi,
                   const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jlo,
                   const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jhi,
                   const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_klo,
                   const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_khi,
                   int fortran_beta_comp, int fortran_vel_comp, int nghost )
{
    BL_PROFILE("MFIXParticleContainer::PICMultiDeposition()");

    const Real strttime = ParallelDescriptor::second();

    if (nlev > 2)
        amrex::Abort("For right now MFIXParticleContainer::PICMultiDeposition can only handle up to 2 levels");

    MultiFab*  drag_ptr[nlev];

    for (int lev = 0; lev < nlev; lev++)
    {
       if (lev == 0 && OnSameGrids(lev, *drag_mf[lev])) {
          // If we are already working with the internal mf defined on the
          // particle_box_array, then we just work with this.
          drag_ptr[lev] = drag_mf[lev].get();

       } else if (lev == 0 && !OnSameGrids(lev, *drag_mf[lev]))  {
          // If beta_mf is not defined on the particle_box_array, then we need
          // to make a temporary here and copy into beta_mf at the end.
          drag_ptr[lev] = new MultiFab(ParticleBoxArray(lev), ParticleDistributionMap(lev),
                                       drag_mf[lev]->nComp(), drag_mf[lev]->nGrow());

       } else {
          // If lev > 0 we make a temporary at the coarse resolution
          BoxArray ba_crse(amrex::coarsen(ParticleBoxArray(lev),this->m_gdb->refRatio(0)));
          drag_ptr[lev] = new MultiFab(ba_crse, ParticleDistributionMap(lev),drag_mf[lev]->nComp(),1);
       }

       // We must have ghost cells for each FAB so that a particle in one grid can spread
       // its effect to an adjacent grid by first putting the value into ghost cells of its
       // own grid.  The mf->sumBoundary call then adds the value from one grid's ghost cell
       // to another grid's valid region.
       if (drag_ptr[lev]->nGrow() < 1)
          amrex::Error("Must have at least one ghost cell when in CalcVolumeFraction");

       drag_ptr[lev]->setVal(0.0,0,4,drag_ptr[lev]->nGrow());
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
                  gm, drag_ptr[lev]->boxArray(), drag_ptr[lev]->DistributionMap(),
                  ngrow, EBSupport::volume);
           flags   = &(crse_factory->getMultiEBCellFlagFab());
        }

        const MultiFab* volfrac = (lev == 0) ? &(ebfactory[lev]->getVolFrac()) : &(crse_factory->getVolFrac());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {

        //const int* lo; // SET_BUT_NOT_USED
        //const int* hi; // SET_BUT_NOT_USED

        FArrayBox local_vol;
         for (ParConstIter pti(*this, lev); pti.isValid(); ++pti) {

            const auto& particles = pti.GetArrayOfStructs();
            const ParticleType* pstruct = particles().dataPtr();
            const long nrp = pti.numParticles();

            FArrayBox& drag_fab = (*drag_ptr[lev])[pti];

#ifdef _OPENMP
            // Note that we actually grow the tilebox rather than calling growntilebox
            //     because we need the overlap even in the interior.
            Box grown_tilebox = pti.tilebox();
            grown_tilebox.grow(1);

            int ncomp = 4;
            local_vol.resize(grown_tilebox,ncomp);

            local_vol = 0.0;

            //lo = grown_tilebox.loVect(); // SET_BUT_NOT_USED
            //hi = grown_tilebox.hiVect(); // SET_BUT_NOT_USED
#else

            //const Box& bx  = drag_fab.box(); // UNUSED_VARIABLE

            //lo = bx.loVect(); // SET_BUT_NOT_USED
            //hi = bx.hiVect(); // SET_BUT_NOT_USED
#endif
            const Box& box = pti.tilebox(); // I need a box without ghosts

            if ((*flags)[pti].getType(box) != FabType::covered )
            {
                auto drag_arr = drag_fab.array();
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

                    amrex::Real pbeta = p.rdata(realData::dragx) / reg_cell_vol;
                    amrex::Real pvx   = p.rdata(realData::velx) * pbeta;
                    amrex::Real pvy   = p.rdata(realData::vely) * pbeta;
                    amrex::Real pvz   = p.rdata(realData::velz) * pbeta;

                    for (int ii = -1; ii <= 0; ++ii) {
                        for (int jj = -1; jj <= 0; ++jj) {
                            for (int kk = -1; kk <= 0; ++kk) {
                                if (flagsarr(i+ii,j+jj,k+kk).isCovered())
                                    continue;

                                amrex::Real weight_vol = weights[ii+1][jj+1][kk+1];

                                amrex::Gpu::Atomic::Add(&drag_arr(i+ii,j+jj,k+kk,0),weight_vol*pvx);
                                amrex::Gpu::Atomic::Add(&drag_arr(i+ii,j+jj,k+kk,1),weight_vol*pvy);
                                amrex::Gpu::Atomic::Add(&drag_arr(i+ii,j+jj,k+kk,2),weight_vol*pvz);
                                amrex::Gpu::Atomic::Add(&drag_arr(i+ii,j+jj,k+kk,3),weight_vol*pbeta);
                            }
                        }
                    }
                });
            }

#ifdef _OPENMP
            drag_fab.atomicAdd(local_vol, grown_tilebox, grown_tilebox, 0, 0, drag_fab.nComp());
#endif

         }
        }
    }

    int  src_nghost = 1;
    int dest_nghost = 0;
    for (int lev = 1; lev < nlev; lev++)
    {
        drag_ptr[0]->copy(*drag_ptr[lev],0,0,drag_ptr[0]->nComp(),src_nghost,dest_nghost,gm.periodicity(),FabArrayBase::ADD);
    }

    drag_ptr[0]->SumBoundary(gm.periodicity());

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
            drag_mf[lev]->setVal(0.0);
            amrex::InterpFromCoarseLevel(*drag_mf[lev], time, *drag_ptr[lev-1],
                                         0, 0, 1, Geom(lev-1), Geom(lev),
                                         cphysbc, 0, fphysbc, 0,
                                         ref_ratio, mapper,
                                         bcs, 0);
        }
    }

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled. I believe that we don't
    // need any information in ghost cells so we don't copy those.

    if (drag_ptr[0] != drag_mf[0].get())
    {
       drag_mf[0]->copy(*drag_ptr[0],0,0,drag_mf[0]->nComp());
    }

    for (int lev = 0; lev < nlev; lev++)
    {
       if (drag_ptr[lev] != drag_mf[lev].get())
          delete drag_ptr[lev];
    }

    if (m_verbose > 1) {
      Real stoptime = ParallelDescriptor::second() - strttime;

      ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

      amrex::Print() << "MFIXParticleContainer::PICMultiDeposition time: " << stoptime << '\n';
    }
}


                    
