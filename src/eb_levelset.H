#ifndef _EB_LEVELSET_H_
#define _EB_LEVELSET_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_RealVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_BoxArray.H>
#include <MFIXParticleContainer.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_PolynomialIF.H>
#include <AMReX_EBFabFactory.H>


using namespace amrex;

class LSFactory {
    private:
        std::unique_ptr<MultiFab> ls_phi;
        std::unique_ptr<iMultiFab> ls_valid;
        std::unique_ptr<MultiFab> dummy;
        
        BoxArray phi_ba_refined, particle_ba_refined;

        int amr_lev;
        
        const MFIXParticleContainer * mfix_pc;
        
        const RealVect dx_vect;
        int ls_grid_refinement;

        std::unique_ptr<FArrayBox> eb_facets(const EBFArrayBoxFactory * eb_factory);

    public:
        LSFactory(int lev, int ref, const MFIXParticleContainer * pc, const Real * dx);
        ~LSFactory();

        void update_ebf(const EBFArrayBoxFactory * eb_factory);
        void update_ebis(const EBIndexSpace * eb_is);

        const MultiFab * get_data(){return ls_phi.get();};
        const iMultiFab * get_valid(){return ls_valid.get();};
        const int get_refinement(){return ls_grid_refinement;};

        const BoxArray * get_cc_ba(){return & particle_ba_refined;};
        const BoxArray * get_nd_ba(){return & phi_ba_refined;};
};


class PolynomialDF : public PolynomialIF {
    private:
        int order;

    public:
        PolynomialDF(const Vector<PolyTerm> & a_polynomial,
                     const bool             & a_inside);


        virtual Real value(const RealVect         & a_point,
                           const Vector<PolyTerm> & a_polynomial) const;


        virtual Real value(const RealVect & a_point) const;


        virtual BaseIF * newImplicitFunction() const;

};

#endif
