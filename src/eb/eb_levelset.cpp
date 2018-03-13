#include "eb_levelset.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_RealVect.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_MultiCutFab.H>
#include "AMReX_BoxIterator.H"
#include <AMReX_EBCellFlag.H>
#include <mfix_F.H>


LSFactory::LSFactory(int lev, int ls_ref, int eb_ref, int ls_pad, int eb_pad, const MFIXParticleContainer * pc)
    : amr_lev(lev), ls_grid_ref(ls_ref), eb_grid_ref(eb_ref), ls_grid_pad(ls_pad), eb_grid_pad(eb_pad), mfix_pc(pc),
    dx_vect(AMREX_D_DECL(pc->Geom(lev).CellSize()[0]/ls_ref,
                         pc->Geom(lev).CellSize()[1]/ls_ref,
                         pc->Geom(lev).CellSize()[2]/ls_ref)),
    dx_eb_vect(AMREX_D_DECL(pc->Geom(lev).CellSize()[0]/eb_ref,
                            pc->Geom(lev).CellSize()[1]/eb_ref,
                            pc->Geom(lev).CellSize()[2]/eb_ref))
{
    // DistributionMapping inherited from MFIXParticleContainer
    const DistributionMapping & dm = mfix_pc -> ParticleDistributionMap(amr_lev);
    // Init geometry over which the level set and EB are defined, note that BoxArrays and DistributionMapping are
    // inherited from the MFIXParticleContainer => the DistributionMapping above applies to refined BoxArrays. Also
    // note that the Geometry objects are updated (for periodic fill operations)
    init_geom();

    // Initialize MultiFab pointers storing level-set data
    //    -> ls_phi:   nodal MultiFab storing signed distance function to the nearest wall
    //    -> ls_valid: cell-centered iMultiFab storing integer flags assuming the following values:
    //         -1 : not all nodes around cell have been initialized
    //          0 : none of the cell's neighbours contain negative vlaues of ls_phi on its nodes
    //          1 : the cell is in the neighbourhood of phi < 0

    ls_grid  = std::unique_ptr<MultiFab>(new MultiFab);
    ls_valid = std::unique_ptr<iMultiFab>(new iMultiFab);
    // Temporary MultiFab used for generating EB factories.
    eb_grid = std::unique_ptr<MultiFab>(new MultiFab);

    // Define ls_grid and ls_valid, growing them by ls_pad
    // Note: box arrays (such as ls_ba) are initialized in init_box()
    ls_grid->define(ls_ba, dm, 1, ls_pad);
    ls_valid->define(ls_ba, dm, 1, ls_pad);
    ls_valid->setVal(-1);

    // Define eb_grid, growing it by eb_pad
    eb_grid->define(eb_ba, dm, 1, eb_pad);

    // Initialize by setting all ls_phi = huge(c_real)
    for(MFIter mfi( * ls_grid, true); mfi.isValid(); ++mfi){
        Box tile_box   = mfi.growntilebox();
        auto & ls_tile = (* ls_grid)[mfi];

        // Initialize in fortran land
        init_levelset(tile_box.loVect(), tile_box.hiVect(),
                      ls_tile.dataPtr(), ls_tile.loVect(),  ls_tile.hiVect());
    }
}


LSFactory::~LSFactory() {
    ls_grid.reset();
    ls_valid.reset();
    eb_grid.reset();
}


void LSFactory::update_ba() {
    // Refined versions of both the cell-centered (particle) and nodal (phi) BoxArrays
    // Note: BoxArrays are inherited from MFIXParticleContainer => the DistributionMapping inherited from
    // MFIXParticleContainer still applies to the refined BoxArrays
    const BoxArray & particle_ba = mfix_pc->ParticleBoxArray(amr_lev);
    const BoxArray & phi_ba      = amrex::convert(particle_ba, IntVect{1,1,1});

    ls_ba = phi_ba;
    ls_ba.refine(ls_grid_ref);

    cc_ba = particle_ba;
    cc_ba.refine(ls_grid_ref);

    eb_ba = particle_ba;
    eb_ba.refine(eb_grid_ref);
}



void LSFactory::init_geom() {
    // Initialize Geometry objects for the level set and the EB, note that the Geometry objects reflect the refined (and
    // padded) box arrays, preventing periodic fill operations from "spilling over" from refined/padded indices.
    update_ba();

    geom_ls = LSUtility::make_ls_geometry(*this);
    geom_eb = LSUtility::make_eb_geometry(*this);
}



std::unique_ptr<Vector<Real>> LSFactory::eb_facets(const EBFArrayBoxFactory & eb_factory) {
    // 1-D list of eb-facet data. Format:
    // { px_1, py_1, pz_1, nx_1, ny_1, nz_1, px_2, py_2, ... , nz_N }
    //   ^                 ^
    //   |                 +---- {nx, ny, nz} is the normal vector pointing _towards_ the facet
    //   +-----------------------{px, py, pz} is the position vector of the facet centre
    std::unique_ptr<Vector<Real>> facet_list;


    const DistributionMapping & dm = mfix_pc->ParticleDistributionMap(amr_lev);


    /***************************************************************************
     *                                                                         *
     * Access EB Cut-Cell data:                                                *
     *                                                                         *
     ***************************************************************************/

    MultiFab dummy(eb_ba, dm, 1, eb_grid_pad, MFInfo(), eb_factory);
    // Area fraction data
    std::array<const MultiCutFab*, AMREX_SPACEDIM> areafrac = eb_factory.getAreaFrac();
    // EB boundary-centre data
    const MultiCutFab * bndrycent = &(eb_factory.getBndryCent());


    /***************************************************************************
     *                                                                         *
     * Compute normals data (which are stored on MultiFab over the eb_ba Grid) *
     *                                                                         *
     ***************************************************************************/

    MultiFab normal(eb_ba, dm, 3, eb_grid_pad);

    // while computing normals, count EB-facets
    int n_facets = 0;

    for(MFIter mfi(normal, true); mfi.isValid(); ++mfi) {
        Box tile_box = mfi.growntilebox();
        const int * lo = tile_box.loVect();
        const int * hi = tile_box.hiVect();

        const auto & sfab = dynamic_cast <EBFArrayBox const&>(dummy[mfi]);
        const auto & flag = sfab.getEBCellFlagFab();

        // Need to count number of eb-facets (in order to allocate FArrayBox)
        count_eb_facets(lo, hi, flag.dataPtr(), flag.loVect(), flag.hiVect(), & n_facets);

        // Target for compute_normals(...)
        auto & norm_tile = normal[mfi];
        // Area fractions in x, y, and z directions
        const auto & af_x_tile = (* areafrac[0])[mfi];
        const auto & af_y_tile = (* areafrac[1])[mfi];
        const auto & af_z_tile = (* areafrac[2])[mfi];

        //if (flag.getType(amrex::grow(tile_box,1)) == FabType::singlevalued) {
        if (flag.getType(tile_box) == FabType::singlevalued) {
           BL_PROFILE_VAR("compute_normals()", compute_normals);
           compute_normals(lo,                  hi,
                           flag.dataPtr(),      flag.loVect(),      flag.hiVect(),
                           norm_tile.dataPtr(), norm_tile.loVect(), norm_tile.hiVect(),
                           af_x_tile.dataPtr(), af_x_tile.loVect(), af_x_tile.hiVect(),
                           af_y_tile.dataPtr(), af_y_tile.loVect(), af_y_tile.hiVect(),
                           af_z_tile.dataPtr(), af_z_tile.loVect(), af_z_tile.hiVect());
           BL_PROFILE_VAR_STOP(compute_normals);
        }
    }

    normal.FillBoundary(geom_eb.periodicity());


    /***************************************************************************
     *                                                                         *
     * Compute EB-facet centres data (which are stored in a 1-D array)         *
     *                                                                         *
     ***************************************************************************/


    facet_list = std::unique_ptr<Vector<Real>>(new Vector<Real>(6 * n_facets));

    int c_facets = 0;
    for(MFIter mfi(normal, true); mfi.isValid(); ++mfi) {
        Box tile_box = mfi.growntilebox();

        const auto & sfab = dynamic_cast <EBFArrayBox const&>(dummy[mfi]);
        const auto & flag = sfab.getEBCellFlagFab();

        const auto & norm_tile =normal[mfi];
        const auto & bcent_tile = (* bndrycent)[mfi];

        int facet_list_size = facet_list->size();

        eb_as_list(tile_box.loVect(),     tile_box.hiVect(),    & c_facets,
                   flag.dataPtr(),        flag.loVect(),        flag.hiVect(),
                   norm_tile.dataPtr(),   norm_tile.loVect(),   norm_tile.hiVect(),
                   bcent_tile.dataPtr(),  bcent_tile.loVect(),  bcent_tile.hiVect(),
                   facet_list->dataPtr(), & facet_list_size,
                   dx_eb_vect.dataPtr());
    }
    return facet_list;
}


std::unique_ptr<MultiFab> LSFactory::ebis_impfunc(const EBIndexSpace & eb_is) {
    std::unique_ptr<MultiFab> mf_impfunc = std::unique_ptr<MultiFab>(new MultiFab);
    const DistributionMapping & dm = mfix_pc->ParticleDistributionMap(amr_lev);
    mf_impfunc->define(ls_ba, dm, 1, ls_grid_pad);

    for(MFIter mfi(* mf_impfunc, true); mfi.isValid(); ++ mfi)
        eb_is.fillNodeFarrayBoxFromImplicitFunction((* mf_impfunc)[mfi], dx_vect);

    mf_impfunc->FillBoundary(geom_ls.periodicity());
    return mf_impfunc;
}


void LSFactory::update(const MultiFab & ls_in) {
    for(MFIter mfi( * ls_grid, true); mfi.isValid(); ++mfi){
        Box tile_box = mfi.growntilebox();

        const auto & ls_in_tile = ls_in[mfi];
        auto & v_tile = (* ls_valid)[mfi];
        auto & ls_tile = (* ls_grid)[mfi];

        update_levelset(tile_box.loVect(),    tile_box.hiVect(),
                        ls_in_tile.dataPtr(), ls_in_tile.loVect(), ls_in_tile.hiVect(),
                        v_tile.dataPtr(),     v_tile.loVect(),     v_tile.hiVect(),
                        ls_tile.dataPtr(),    ls_tile.loVect(),    ls_tile.hiVect(),
                        dx_vect.dataPtr(),    & ls_grid_pad);
    }

    ls_grid->FillBoundary(geom_ls.periodicity());
    ls_valid->FillBoundary(geom_ls.periodicity());
}


void LSFactory::regrid(){
    // Regrids the level-set data whenever the MFIXParticleContainer's DistributionMapping has changed:
    //      -> Loads the updated DistributionMapping from MFIXParticleContainer
    //      -> Rebuilds the nodal levelset (ls_ba), cell-centered valid (cc_ba), and eb (eb_ba) BoxArrays
    //          -> ls_ba, cc_ba, and eb_ba are all inherited from MFIXParticleContainer::ParticleBoxArray
    //  =>  make sure that the AMReX level of the MFIXParticleContainer has been regridded before calling this method
    const DistributionMapping & dm = mfix_pc -> ParticleDistributionMap(amr_lev);
    update_ba();

    int ng = ls_grid_pad;
    std::unique_ptr<MultiFab> ls_grid_new(new MultiFab(ls_ba, dm, 1, ng));

    ls_grid_new->copy(* ls_grid, 0, 0, 1, ng, ng);
    ls_grid_new->FillBoundary(geom_ls.periodicity());
    ls_grid = std::move(ls_grid_new);

    std::unique_ptr<iMultiFab> ls_valid_new(new iMultiFab(ls_ba, dm, 1, ng));

    ls_valid_new->copy(* ls_valid, 0, 0, 1, ng, ng);
    ls_valid_new->FillBoundary(geom_ls.periodicity());
    ls_valid = std::move(ls_valid_new);
}


void LSFactory::update_ebf(const EBFArrayBoxFactory & eb_factory, const EBIndexSpace & eb_is) {

    // Generate facets (TODO: in future these can also be provided by user)
    std::unique_ptr<Vector<Real>> facets = eb_facets(eb_factory);
    int len_facets = facets->size();
    // Generate implicit function (used to determine the interior of EB)
    std::unique_ptr<MultiFab> impfunct = ebis_impfunc(eb_is);
    impfunct->FillBoundary(geom_ls.periodicity());

    // Local MultiFab storing level-set data for this eb_factory
    MultiFab eb_ls;
    iMultiFab eb_valid;

    const DistributionMapping & dm = mfix_pc -> ParticleDistributionMap(amr_lev);
    eb_ls.define(ls_ba, dm, 1, ls_grid_pad);
    eb_valid.define(ls_ba, dm, 1, ls_grid_pad);
    eb_valid.setVal(0);

    // Fill local MultiFab with eb_factory's level-set data. Note the role of eb_valid:
    //  -> eb_valid = 1 if the corresponding eb_ls location could be projected onto the eb-facets
    //  -> eb_valid = 0 if eb_ls is the fall-back (euclidian) distance to the nearest eb-facet
    for(MFIter mfi(eb_ls, true); mfi.isValid(); ++mfi){
        Box tile_box = mfi.growntilebox();
        const int * lo = tile_box.loVect();
        const int * hi = tile_box.hiVect();

        auto & v_tile = eb_valid[mfi];
        auto & ls_tile = eb_ls[mfi];
        const auto & if_tile = (* impfunct)[mfi];

        fill_levelset_eb(lo,                hi,
                         facets->dataPtr(), & len_facets,
                         v_tile.dataPtr(),  v_tile.loVect(),  v_tile.hiVect(),
                         ls_tile.dataPtr(), ls_tile.loVect(), ls_tile.hiVect(),
                         dx_vect.dataPtr(), dx_eb_vect.dataPtr());

        validate_levelset(lo,                hi,               & ls_grid_ref,
                          if_tile.dataPtr(), if_tile.loVect(), if_tile.hiVect(),
                          v_tile.dataPtr(),  v_tile.loVect(),  v_tile.hiVect(),
                          ls_tile.dataPtr(), ls_tile.loVect(), ls_tile.hiVect());

    }

    // Update LSFactory using local eb level-set
    update(eb_ls);
}


void LSFactory::update_ebis(const EBIndexSpace & eb_is) {
    std::unique_ptr<MultiFab> mf_impfunc = ebis_impfunc(eb_is);

    // GeometryService convetion:
    //      -- implicit_function(r) < 0 : r in fluid (outside of EB)
    //      -- implicit_function(r) > 0 : r not in fluid (inside EB)
    //   => If implicit_function is a signed-distance function, we need to invert sign
    for(MFIter mfi( * mf_impfunc, true); mfi.isValid(); ++ mfi){
        FArrayBox & a_fab = (* mf_impfunc)[mfi];
        for(BoxIterator bit(mfi.growntilebox()); bit.ok(); ++bit)
            a_fab(bit(), 0) = - a_fab(bit(), 0);
    }

    update(* mf_impfunc);
}


PolynomialDF::PolynomialDF(const Vector<PolyTerm> & a_polynomial, const bool & a_inside)
             :PolynomialIF(a_polynomial, a_inside)
{
    int size = a_polynomial.size();
    order = 0;
    for(int iterm = 0; iterm < size; iterm++){
        int cur_order = 0;
        for(int idir = 0; idir < SpaceDim; idir++){
            cur_order += a_polynomial[iterm].powers[idir];
        }
        order = cur_order > order ? cur_order : order;
    }
}


Real PolynomialDF::value(const RealVect & a_point, const Vector<PolyTerm> & a_polynomial) const {
    Real retval = 0;

    int size = a_polynomial.size();
    Real terms[order + 1];
    for(int i = 0; i <= order; i++)
        terms[i] = 0;

    // Collect like powers as terms
    for(int iterm = 0; iterm < size; iterm++){
        PolyTerm pterm = a_polynomial[iterm];
        Real coeff     = pterm.coef;
        Real cur       = coeff;
        int cur_order  = 0;
        for(int idir = 0; idir < SpaceDim; idir++){
            cur *= pow(a_point[idir], pterm.powers[idir]);
            cur_order += pterm.powers[idir];
        }
        terms[cur_order] += cur;
    }

    // Evaluate distance function term-by-term:
    Real sg_t0 = terms[0] < 0 ? -1. : 1.;
    retval = sg_t0 * sqrt(sg_t0 * terms[0]); // compatibility for standard PolynomialIF:
                                             // spheres, cylinders have r^2 as 0-order term
                                             // -> hence take sqrt on terms[0] and itterate starting from term 1
    for(int i = 1; i <= order; i++){
        retval += pow(terms[i], 1./(double) i);
    }

    // Change the sign to change inside to outside
    if (!m_inside)
      retval = -retval;

    return retval;
};


Real PolynomialDF::value(const RealVect & a_point) const {
    return value(a_point,m_polynomial);
}


BaseIF * PolynomialDF::newImplicitFunction() const {
    PolynomialIF * polynomialPtr = new PolynomialDF(m_polynomial,
                                                    m_inside);

    return static_cast<BaseIF*>(polynomialPtr);
}
