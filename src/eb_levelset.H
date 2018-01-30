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

        // Hack: figure out how to get level eet info out of GeometryShop
        // For now thoughm we manually hack a cylinder level-set function
        //double xy_centre[2] = {0.0016, 0.0016};
        //double R2 = 2.56e-6;
        
        //double value_wall(const double pos[3]);
        //double value_bottom(const double pos[3]);

    public:
        LSFactory(int lev, const MFIXParticleContainer * pc, const EBFArrayBoxFactory * ebfactory);
        ~LSFactory();

        void update(MultiFab * dummy);
        const MultiFab * get_data(){return ls_phi.get();};

};



#endif
