#ifndef _EB_LEVELSET_H_
#define _EB_LEVELSET_H_

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_BoxArray.H>
#include <AMReX_EBFabFactory.H>
#include <MFIXParticleContainer.H>

using namespace amrex;

class LSFactory {
    private:
        std::unique_ptr<MultiFab> ls_phi;
        std::unique_ptr<iMultiFab> ls_valid;
        int amr_lev;
        const MFIXParticleContainer * mfix_pc;
        const EBFArrayBoxFactory * eb_factory;
        int ls_grid_refinement;

    public:
        LSFactory(int lev, int ref, const MFIXParticleContainer * pc, const EBFArrayBoxFactory * ebfactory);
        ~LSFactory();

        void update(MultiFab * dummy);
        const MultiFab * get_data(){return ls_phi.get();};
        const iMultiFab * get_valid(){return ls_valid.get();};
        const int get_refinement(){return ls_grid_refinement;};
};



#endif
