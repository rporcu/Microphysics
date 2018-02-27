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
#include <AMReX_EBCellFlag.H>


using namespace amrex;



class LSFactory {
    private:
        // Data structures holding level set information
        std::unique_ptr<MultiFab>  ls_grid;
        std::unique_ptr<iMultiFab> ls_valid;

        // Dummy structures encloding (padded) EB index space
        std::unique_ptr<MultiFab> dummy, eb_grid;
        BoxArray eb_ba;
        const int eb_grid_ref, eb_grid_pad;

        // Refined box array (over which the level-set MultiFab is constructed)
        BoxArray ls_ba, cc_ba;
        const RealVect dx_vect, dx_eb_vect;
        const int ls_grid_ref, ls_grid_pad;

        // AMR information
        const int amr_lev;
        const MFIXParticleContainer * mfix_pc;


        // Construct eb_facets list from eb_factory
        //  -> eb_facets are currently stored in a 1-D FArrayBox
        //     (this is a bit of a hack, but we need a fortran-compatible list structure)
        std::unique_ptr<Vector<Real>> eb_facets(const EBFArrayBoxFactory & eb_factory);

        // Construct nodal MultiFab containing EBIS Implicit function velues
        std::unique_ptr<MultiFab> ebis_impfunc(const EBIndexSpace & eb_is);

        // Unpdate internal `ls_phi` and `ls_valid` based on the min of `ls_phi` and `ls_in`
        void update(const MultiFab & ls_in);

        void init_box();

    public:
        LSFactory(int lev, int ls_ref, int eb_ref, int ls_pad, int eb_pad, const MFIXParticleContainer * pc);
        ~LSFactory();

        void update_ebf(const EBFArrayBoxFactory & eb_factory, const EBIndexSpace & eb_is);
        void update_ebis(const EBIndexSpace & eb_is);

        const MultiFab * get_data() const {return ls_grid.get();};
        const iMultiFab * get_valid() const {return ls_valid.get();};
        const int get_ls_ref() const {return ls_grid_ref;};
        const int get_ls_pad() const {return ls_grid_pad;};
        const int get_eb_ref() const {return eb_grid_ref;};
        const int get_eb_pad() const {return eb_grid_pad;};
        const int get_amr_level() const {return amr_lev;};

        const BoxArray * get_ls_ba() const {return & ls_ba;};
        const BoxArray * get_cc_ba() const {return & cc_ba;};
        const BoxArray * get_eb_ba() const {return & eb_ba;};

        const MFIXParticleContainer & get_pc() const {return * mfix_pc;};
};



class LSUtility {
    public:
        static Geometry make_ls_geometry(const LSFactory & level_set){
            const MFIXParticleContainer & pc = level_set.get_pc();
            const int amr_lev = level_set.get_amr_level();

            Box dom_ls = pc.Geom(amr_lev).Domain();// geom[lev].Domain();
            dom_ls.refine(level_set.get_ls_ref());
            Geometry geom_ls(dom_ls);
            Box dom_ls_g = geom_ls.Domain();
            dom_ls_g.grow(level_set.get_ls_pad());
            geom_ls.Domain(dom_ls_g);
            amrex::Print() << "Refined geom_ls: " << geom_ls << std::endl;

            return geom_ls;
        };


        static Geometry make_eb_geometry(const LSFactory & level_set){
            const MFIXParticleContainer & pc = level_set.get_pc();
            const int amr_lev = level_set.get_amr_level();

            Box dom_eb = pc.Geom(amr_lev).Domain();// geom[lev].Domain();
            dom_eb.refine(level_set.get_eb_ref());
            Geometry geom_eb(dom_eb);
            Box dom_eb_g = geom_eb.Domain();
            dom_eb_g.grow(level_set.get_eb_pad());
            geom_eb.Domain(dom_eb_g);
            amrex::Print() << "EB geom_eb: " << geom_eb << std::endl;

            return geom_eb;
        };


        //static EBFArrayBoxFactory make_eb_factory(const LSFactory & level_set){
        //    Geometry geom_eb_g = LSUtility::make_eb_geometry(level_set);
        //    //Box dom_eb_g = geom_eb_g.Domain();
        //    //dom_eb_g.grow(4);
        //    //geom_eb_g.Domain(dom_eb_g);
        //    //amrex::Print() << "Grown geom_ls_g: " << geom_eb_g << std::endl;
        //    //BoxArray ba_g = * level_set.get_eb_ba();
        //    //ba_g.grow(4);
        //
        //
        //}


        static void PrintFlagType(const Box & tile_box, const EBCellFlagFab & flag){
            if(flag.getType(tile_box) == FabType::covered)
                amrex::Print() << "flag covered" << std::endl;
            else if(flag.getType(tile_box) == FabType::regular)
                amrex::Print() << "flag regular" << std::endl;
            else if(flag.getType(tile_box) == FabType::singlevalued)
                amrex::Print() << "flag singlevalued" << std::endl;
            else if(flag.getType(tile_box) == FabType::multivalued)
                amrex::Print() << "flag multivalued" << std::endl;
            else if(flag.getType(tile_box) == FabType::undefined)
                amrex::Print() << "flag undefined" << std::endl;
            else
                amrex::Print() << "flag ill-formed" << std::endl;
        }
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
