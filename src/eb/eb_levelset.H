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

        // MultiFab pointers storing level-set data
        //    -> ls_phi:   nodal MultiFab storing signed distance function to the nearest wall
        //    -> ls_valid: cell-centered iMultiFab storing integer flags assuming the following values:
        //         -1 : not all nodes around cell have been initialized
        //          0 : none of the cell's neighbours contain negative vlaues of ls_phi on its nodes
        //          1 : the cell is in the neighbourhood of phi < 0

        std::unique_ptr<MultiFab>  ls_grid;
        std::unique_ptr<iMultiFab> ls_valid;

        // AMR information
        const int amr_lev;

        // Refined box array (over which the level-set MultiFab is constructed)
        BoxArray ls_ba, cc_ba, eb_ba;
        const int ls_grid_ref, eb_grid_ref;
        // Data structures encloding (padded) EB index space
        std::unique_ptr<MultiFab> eb_grid;

        // Geometry objects reflecting refined (and padded) BoxArrays. This is
        // required to prevent PeriodicFill operations from "spilling" refined
        // indices over periodic boundaries
        Geometry geom_ls, geom_eb;

        // Padding of the leve-set and EB grids:
        const int ls_grid_pad, eb_grid_pad;

        // Particle container (used to syncronize DistributionMapping)
        const MFIXParticleContainer * mfix_pc;

        // Cell size of the refined grids
        const RealVect dx_vect, dx_eb_vect;

        // Construct eb_facets list from eb_factory
        std::unique_ptr<Vector<Real>> eb_facets(const EBFArrayBoxFactory & eb_factory);

        // Construct nodal MultiFab containing EBIS Implicit function velues
        std::unique_ptr<MultiFab> ebis_impfunc(const EBIndexSpace & eb_is);

        void update_ba();
        void init_geom();

    public:
        LSFactory(int lev, int ls_ref, int eb_ref, int ls_pad, int eb_pad, const MFIXParticleContainer * pc);
        LSFactory(const LSFactory & other);
        ~LSFactory();

        void regrid();

        void invert();

        void intersection_ebf(const EBFArrayBoxFactory & eb_factory, const EBIndexSpace & eb_is);
        void union_ebf(const EBFArrayBoxFactory & eb_factory, const EBIndexSpace & eb_is);
        void intersection_ebis(const EBIndexSpace & eb_is);
        void union_ebis(const EBIndexSpace & eb_is);

        // Unpdate internal `ls_phi` and `ls_valid` based on the min of
        // `ls_phi` and `ls_in`
        void update_intersection(const MultiFab & ls_in);
        void update_union(const MultiFab & ls_in);

        // Return level-set data
        const MultiFab *  get_data() const {return ls_grid.get();};
        const iMultiFab * get_valid() const {return ls_valid.get();};

        // Copy level-set data
        std::unique_ptr<MultiFab> copy_data() const;
        std::unique_ptr<iMultiFab> copy_valid() const;

        // Return grid parameters
        const int get_ls_ref() const {return ls_grid_ref;};
        const int get_ls_pad() const {return ls_grid_pad;};
        const int get_eb_ref() const {return eb_grid_ref;};
        const int get_eb_pad() const {return eb_grid_pad;};

        // Return AMR level (wrt the MFIXParticleContainer pointer (pc)
        const int get_amr_level() const {return amr_lev;};

        // Return BoxArrays:
        //    -> ls: levelset (nodal)
        //    -> cc: cell-centered (dual to ls) used by valid
        //    -> eb: BoxArray used by EBIS (cell-centered)
        const BoxArray & get_ls_ba() const {return ls_ba;};
        const BoxArray & get_cc_ba() const {return cc_ba;};
        const BoxArray & get_eb_ba() const {return eb_ba;};

        const MFIXParticleContainer & get_pc() const {return * mfix_pc;};
};



class LSUtility {
    public:
        static Geometry make_ls_geometry(const LSFactory & level_set){
            const MFIXParticleContainer & pc = level_set.get_pc();
            const int amr_lev = level_set.get_amr_level();

            // Set up refined geometry
            Box dom_ls = pc.Geom(amr_lev).Domain();
            dom_ls.refine(level_set.get_ls_ref());
            Geometry geom_ls(dom_ls);

            // This is a bit of a hack: we want to grow the geometry's domain
            // box, without influencing the cell size.  I.e. we want
            // Geometry::Domain() to extend beyond probLo and probHi.
            Box dom_ls_g = geom_ls.Domain();
            dom_ls_g.grow(level_set.get_ls_pad());
            geom_ls.Domain(dom_ls_g);
            //amrex::Print() << "Refined geom_ls: " << geom_ls << std::endl;

            return geom_ls;
        }


        static Geometry make_eb_geometry(const LSFactory & level_set){
            const MFIXParticleContainer & pc = level_set.get_pc();
            const int amr_lev = level_set.get_amr_level();

            // Set up refined geometry
            Box dom_eb = pc.Geom(amr_lev).Domain();
            dom_eb.refine(level_set.get_eb_ref());
            Geometry geom_eb(dom_eb);

            // This is a bit of a hack: we want to grow the geometry's domain
            // box, without influencing the cell size.  I.e. we want
            // Geometry::Domain() to extend beyond probLo and probHi.
            Box dom_eb_g = geom_eb.Domain();
            dom_eb_g.grow(level_set.get_eb_pad());
            geom_eb.Domain(dom_eb_g);
            //amrex::Print() << "EB geom_eb: " << geom_eb << std::endl;

            return geom_eb;
        }


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
