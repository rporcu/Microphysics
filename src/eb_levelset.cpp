#include "eb_levelset.H"

#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_MultiCutFab.H>

#include "AMReX_BoxIterator.H"

#include <mfix_F.H>


LSFactory::LSFactory(int lev, int ref, const MFIXParticleContainer * pc, const Real * dx)
    : amr_lev(lev), ls_grid_refinement(ref), mfix_pc(pc),
    dx_vect(AMREX_D_DECL(mfix_pc->Geom(amr_lev).CellSize()[0]/ref,
                         mfix_pc->Geom(amr_lev).CellSize()[1]/ref,
                         mfix_pc->Geom(amr_lev).CellSize()[2]/ref))
{
    // initialize MultiFab pointers storing level-set data
    //    -> ls_phi:   nodal MultiFab storing signed distance function to the nearest wall
    //    -> ls_valid: cell-centered iMultiFab storing integer flags assuming the following values: 
    //         -1 : not all nodes around cell have been initialized
    //          0 : none of the cell's neighbours contain negative vlaues of ls_phi on its nodes
    //          1 : the cell is in the neighbourhood of phi < 0
    //
    ls_phi   = std::unique_ptr<MultiFab>(new MultiFab);
    ls_valid = std::unique_ptr<iMultiFab>(new iMultiFab);


    // Keep refined versions of both the cell-centered (particle) and nodal (phi) BoxArrays
    const BoxArray & particle_ba   = mfix_pc -> ParticleBoxArray(lev);
    const BoxArray & phi_ba        = amrex::convert(particle_ba, IntVect{1,1,1});
    const DistributionMapping & dm = mfix_pc -> ParticleDistributionMap(lev);

    phi_ba_refined = phi_ba;
    phi_ba_refined.refine(ref);
    particle_ba_refined = particle_ba;
    particle_ba_refined.refine(ref);

    // Define ls_phi and ls_valid, growing them by twice the refinement ratio
    ls_phi -> define(phi_ba_refined, dm, 1, ref); // grow by 1 x ref, as it is nodal
    ls_valid -> define(particle_ba_refined, dm, 1, 2 * ref);

    // Initialize by setting all ls_valid = -1, and ls_phi = huge(c_real)
    for(MFIter mfi( * ls_phi, true); mfi.isValid(); ++mfi){
        Box tile_box   = mfi.tilebox();
        auto & v_tile  = (* ls_valid)[mfi];
        auto & ls_tile = (* ls_phi)[mfi];

        // Initialize in fortran land
        init_levelset(tile_box.loVect(), tile_box.hiVect(), & ls_grid_refinement,
                      v_tile.dataPtr(),  v_tile.loVect(),   v_tile.hiVect(),
                      ls_tile.dataPtr(), ls_tile.loVect(),  ls_tile.hiVect());
    }

    // temporary dummy variable used for storing eb-factory flag bits
    dummy = std::unique_ptr<MultiFab>(new MultiFab);
}


LSFactory::~LSFactory() {
    ls_phi.reset();
    ls_valid.reset();
}


std::unique_ptr<FArrayBox> LSFactory::eb_facets(const EBFArrayBoxFactory * eb_factory) {
    // Container for normal data
    std::unique_ptr<MultiFab> normal = std::unique_ptr<MultiFab>(new MultiFab);
    std::unique_ptr<FArrayBox> facet_list;

    // Only call the routine for wall collisions if the box has a wall
    if (eb_factory != NULL) {
        dummy->define(particle_ba_refined, mfix_pc->ParticleDistributionMap(amr_lev),
                1, 0, MFInfo(), * eb_factory);
        std::array<const MultiCutFab*, AMREX_SPACEDIM> areafrac = eb_factory->getAreaFrac();
        const MultiCutFab * bndrycent = &(eb_factory->getBndryCent());

        int n_facets = 0;

        // We pre-compute the normals
        normal->define(particle_ba_refined, mfix_pc->ParticleDistributionMap(amr_lev), 3, 2);

        for(MFIter mfi(* normal, true); mfi.isValid(); ++mfi) {
            Box tile_box = mfi.tilebox();
            const int * lo = tile_box.loVect();
            const int * hi = tile_box.hiVect();

            const auto & sfab = dynamic_cast <EBFArrayBox const&>((*dummy)[mfi]);
            const auto & flag = sfab.getEBCellFlagFab();

            // Need to count number of eb-facets (in order to allocate FArrayBox)
            count_eb_facets(lo, hi, flag.dataPtr(), flag.loVect(), flag.hiVect(), & n_facets);

            // Target for compute_normals(...)
            auto & norm_tile = (* normal)[mfi];
            // Area fractions in x, y, and z directions
            const auto & af_x_tile = (* areafrac[0])[mfi];
            const auto & af_y_tile = (* areafrac[1])[mfi];
            const auto & af_z_tile = (* areafrac[2])[mfi];

            if (flag.getType(amrex::grow(tile_box,1)) == FabType::singlevalued) {
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
        normal->FillBoundary(mfix_pc->Geom(0).periodicity());

        Box dom_list(IntVect{0,0,0}, IntVect{n_facets,0,0});
        facet_list = std::unique_ptr<FArrayBox>(new FArrayBox(dom_list, 6));

        int c_facets = 0;
        for(MFIter mfi( * normal, true); mfi.isValid(); ++mfi) {
            Box tile_box = mfi.tilebox();

            const auto & sfab = dynamic_cast <EBFArrayBox const&>((*dummy)[mfi]);
            const auto & flag = sfab.getEBCellFlagFab();

            const auto & norm_tile = (* normal)[mfi];
            const auto & bcent_tile = (* bndrycent)[mfi];

            // Target: facet_list
            eb_as_list(tile_box.loVect(),     tile_box.hiVect(),    & c_facets,
                       flag.dataPtr(),        flag.loVect(),        flag.hiVect(),
                       norm_tile.dataPtr(),   norm_tile.loVect(),   norm_tile.hiVect(),
                       bcent_tile.dataPtr(),  bcent_tile.loVect(),  bcent_tile.hiVect(),
                       facet_list->dataPtr(), facet_list->loVect(), facet_list->hiVect(),
                       mfix_pc->Geom(amr_lev).CellSize(), & ls_grid_refinement);
        }
    }
    return facet_list;
}


std::unique_ptr<MultiFab> LSFactory::ebis_impfunc(const EBIndexSpace * eb_is) {
    std::unique_ptr<MultiFab> mf_impfunc = std::unique_ptr<MultiFab>(new MultiFab);
    const DistributionMapping & dm = mfix_pc->ParticleDistributionMap(amr_lev);
    mf_impfunc->define(phi_ba_refined, dm, 1, ls_grid_refinement);

    for(MFIter mfi(* mf_impfunc, true); mfi.isValid(); ++ mfi)
        eb_is->fillNodeFarrayBoxFromImplicitFunction((* mf_impfunc)[mfi], dx_vect);

    mf_impfunc->FillBoundary(mfix_pc -> Geom(0).periodicity());
    return mf_impfunc;
}


void LSFactory::update(const MultiFab * ls_in) {
    for(MFIter mfi( * ls_phi, true); mfi.isValid(); ++mfi){
        Box tile_box = mfi.tilebox();

        const auto & ls_in_tile = (* ls_in)[mfi];
        auto & v_tile = (* ls_valid)[mfi];
        auto & ls_tile = (* ls_phi)[mfi];

        update_levelset(tile_box.loVect(),    tile_box.hiVect(),
                        ls_in_tile.dataPtr(), ls_in_tile.loVect(), ls_in_tile.hiVect(),
                        v_tile.dataPtr(),     v_tile.loVect(),     v_tile.hiVect(),
                        ls_tile.dataPtr(),    ls_tile.loVect(),    ls_tile.hiVect(),
                        mfix_pc->Geom(amr_lev).CellSize(),         & ls_grid_refinement);
    }

    ls_phi -> FillBoundary(mfix_pc -> Geom(0).periodicity());
    ls_valid -> FillBoundary(mfix_pc -> Geom(0).periodicity());
}


void LSFactory::update_ebf(const EBFArrayBoxFactory * eb_factory, const EBIndexSpace * eb_is) {
    if(eb_factory == NULL) return;

    // Generate facets (TODO: in future these can also be provided by user)
    std::unique_ptr<FArrayBox> facets = eb_facets(eb_factory);
    // Generate implicit function (used to determine the interior of EB)
    std::unique_ptr<MultiFab> impfunct = ebis_impfunc(eb_is);
    impfunct->FillBoundary(mfix_pc->Geom(0).periodicity());

    // Local MultiFab storing level-set data for this eb_factory
    MultiFab eb_ls;
    iMultiFab eb_valid;

    const DistributionMapping & dm = mfix_pc -> ParticleDistributionMap(amr_lev);
    eb_ls.define(phi_ba_refined, dm, 1, ls_grid_refinement);
    eb_valid.define(phi_ba_refined, dm, 1, ls_grid_refinement);
    eb_valid.setVal(0);

    // Fill local MultiFab with eb_factory's level-set data. Note the role of eb_valid:
    //  -> eb_valid = 1 if the corresponding eb_ls location could be projected onto the eb-facets
    //  -> eb_valid = 0 if eb_ls is the fall-back (euclidian) distance to the nearest eb-facet
    for(MFIter mfi(eb_ls, true); mfi.isValid(); ++mfi){
        Box tile_box = mfi.tilebox();
        const int * lo = tile_box.loVect();
        const int * hi = tile_box.hiVect();

        const auto & sfab = dynamic_cast<EBFArrayBox const &>((* dummy)[mfi]);
        const auto & flag = sfab.getEBCellFlagFab();

        if(flag.getType(amrex::grow(tile_box, 1)) == FabType::singlevalued){
            auto & v_tile = eb_valid[mfi];
            auto & ls_tile = eb_ls[mfi];
            const auto & if_tile = (* impfunct)[mfi];

            fill_levelset_eb(lo,                hi,
                             facets->dataPtr(), facets->loVect(), facets->hiVect(),
                             v_tile.dataPtr(),  v_tile.loVect(),  v_tile.hiVect(),
                             ls_tile.dataPtr(), ls_tile.loVect(), ls_tile.hiVect(),
                             mfix_pc->Geom(amr_lev).CellSize(),   & ls_grid_refinement);

            validate_levelset(lo,                hi,               & ls_grid_refinement,
                              if_tile.dataPtr(), if_tile.loVect(), if_tile.hiVect(),
                              v_tile.dataPtr(),  v_tile.loVect(),  v_tile.hiVect(),
                              ls_tile.dataPtr(), ls_tile.loVect(), ls_tile.hiVect());

        }
    }

    // Update LSFactory using local eb level-set
    update(& eb_ls);
}


void LSFactory::update_ebis(const EBIndexSpace * eb_is) {
    std::unique_ptr<MultiFab> mf_impfunc = ebis_impfunc(eb_is);

    // GeometryService convetion:
    //      -- implicit_function(r) < 0 : r in fluid (outside of EB)
    //      -- implicit_function(r) > 0 : r not in fluid (inside EB)
    //   => If implicit_function is a signed-distance function, we need to invert sign
    for(MFIter mfi( * mf_impfunc, true); mfi.isValid(); ++ mfi){
        FArrayBox & a_fab = (* mf_impfunc)[mfi];
        for(BoxIterator bit(a_fab.box()); bit.ok(); ++bit)
                a_fab(bit(), 0) = - a_fab(bit(), 0);
    }

    update(mf_impfunc.get());
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
