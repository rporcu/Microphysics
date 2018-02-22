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
        // Data structures holding level set information
        std::unique_ptr<MultiFab>  ls_phi;
        std::unique_ptr<iMultiFab> ls_valid;

        // Dummy structures encloding (padded) EB index space
        std::unique_ptr<MultiFab> dummy;
        std::unique_ptr<MultiFab> eb_grid;
        BoxArray eb_ba_refined;
        const int eb_grid_ref, eb_grid_pad;

        // Refined box array (over which the level-set MultiFab is constructed)
        BoxArray phi_ba_refined, particle_ba_refined;
        const RealVect dx_vect, dx_eb_vect;
        const int ls_grid_ref, ls_grid_pad;

        // AMR information
        const int amr_lev;
        const MFIXParticleContainer * mfix_pc;


        // Construct eb_facets list from eb_factory
        //  -> eb_facets are currently stored in a 1-D FArrayBox
        //     (this is a bit of a hack, but we need a fortran-compatible list structure)
        std::unique_ptr<Vector<Real>> eb_facets(const EBFArrayBoxFactory * eb_factory);

        // Construct nodal MultiFab containing EBIS Implicit function velues
        std::unique_ptr<MultiFab> ebis_impfunc(const EBIndexSpace * eb_is);

        // Unpdate internal `ls_phi` and `ls_valid` based on the min of `ls_phi` and `ls_in`
        void update(const MultiFab * ls_in);

        void init_box();

    public:
        LSFactory(int lev, int ls_ref, int eb_ref, int ls_pad, int eb_pad, const MFIXParticleContainer * pc);
        ~LSFactory();

        void update_ebf(const EBFArrayBoxFactory * eb_factory, const EBIndexSpace * eb_is);
        void update_ebis(const EBIndexSpace * eb_is);

        const MultiFab * get_data(){return ls_phi.get();};
        const iMultiFab * get_valid(){return ls_valid.get();};
        const int get_ls_ref(){return ls_grid_ref;};
        const int get_eb_ref(){return eb_grid_ref;};

        const BoxArray * get_cc_ba(){return & particle_ba_refined;};
        const BoxArray * get_nd_ba(){return & phi_ba_refined;};
        
        const BoxArray * get_eb_ba(){return & eb_ba_refined;}
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
