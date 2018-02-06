#ifndef _EB_LEVELSET_H_
#define _EB_LEVELSET_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_RealVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_BoxArray.H>
//#include <AMReX_EBFabFactory.H>
#include <MFIXParticleContainer.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_PolynomialIF.H>
#include <AMReX_IntersectionIF.H>

using namespace amrex;

class LSFactory {
    private:
        std::unique_ptr<MultiFab> ls_phi;
        std::unique_ptr<iMultiFab> ls_valid;
        int amr_lev;
        const MFIXParticleContainer * mfix_pc;
        // const EBFArrayBoxFactory * eb_factory;
        const EBIndexSpace * eb_is;
        const RealVect dx_vect;
        int ls_grid_refinement;

    public:
        LSFactory(int lev, int ref, const MFIXParticleContainer * pc, 
                  //const EBFArrayBoxFactory * ebfactory
                  const EBIndexSpace * ebis, const Real * dx
                  );
        ~LSFactory();

        void update(MultiFab * dummy);
        const MultiFab * get_data(){return ls_phi.get();};
        const iMultiFab * get_valid(){return ls_valid.get();};
        const int get_refinement(){return ls_grid_refinement;};
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
//
//class  IntersectionDF : public IntersectionIF {
//    public:
//        IntersectionDF(const Vector<BaseIF*> & a_impFuncs);
//
//    virtual Real value(const RealVect& a_point) const;
//
//    virtual BaseIF* newImplicitFunction() const;
//};

#endif
