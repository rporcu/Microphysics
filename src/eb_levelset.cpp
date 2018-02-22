#include "eb_levelset.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_RealVect.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_MultiCutFab.H>

#include "AMReX_BoxIterator.H"

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
    // Initialize MultiFab pointers storing level-set data
    //    -> ls_phi:   nodal MultiFab storing signed distance function to the nearest wall
    //    -> ls_valid: cell-centered iMultiFab storing integer flags assuming the following values: 
    //         -1 : not all nodes around cell have been initialized
    //          0 : none of the cell's neighbours contain negative vlaues of ls_phi on its nodes
    //          1 : the cell is in the neighbourhood of phi < 0
    //
    ls_phi   = std::unique_ptr<MultiFab>(new MultiFab);
    ls_valid = std::unique_ptr<iMultiFab>(new iMultiFab);
    

    // Refined versions of both the cell-centered (particle) and nodal (phi) BoxArrays
    const BoxArray & particle_ba   = mfix_pc -> ParticleBoxArray(lev);
    const BoxArray & phi_ba        = amrex::convert(particle_ba, IntVect{1,1,1});
    const DistributionMapping & dm = mfix_pc -> ParticleDistributionMap(lev);

    ls_ba = phi_ba;
    ls_ba.refine(ls_ref);
    //ls_ba.grow(ls_pad);
    
    cc_ba = particle_ba;
    cc_ba.refine(ls_ref);
    //cc_ba.grow(ls_pad);

    // Define ls_phi and ls_valid, growing them by twice the refinement ratio
    ls_phi->define(ls_ba, dm, 1, ls_pad);
    ls_valid->define(ls_ba, dm, 1, ls_pad + 1);
    ls_valid->setVal(-1);

    // Initialize by setting all ls_valid = -1, and ls_phi = huge(c_real)
    for(MFIter mfi( * ls_phi, true); mfi.isValid(); ++mfi){
        Box tile_box   = mfi.growntilebox();
        auto & v_tile  = (* ls_valid)[mfi];
        auto & ls_tile = (* ls_phi)[mfi];

        // Initialize in fortran land
        init_levelset(tile_box.loVect(), tile_box.hiVect(),
                      ls_tile.dataPtr(), ls_tile.loVect(),  ls_tile.hiVect());
    }

    // Temporary dummy variable used for storing eb-factory flag bits
    dummy = std::unique_ptr<MultiFab>(new MultiFab);
    
    // Temporary MultiFab used for generating EB factories.
    eb_grid = std::unique_ptr<MultiFab>(new MultiFab);
    eb_ba = particle_ba;
    eb_ba.refine(eb_grid_ref);
    //eb_ba.grow(eb_pad);

    eb_grid->define(eb_ba, mfix_pc->ParticleDistributionMap(amr_lev), 1, eb_grid_pad + 1);
}


LSFactory::~LSFactory() {
    ls_phi.reset();
    ls_valid.reset();
    eb_grid.reset();
}


void LSFactory::init_box(){

}


std::unique_ptr<Vector<Real>> LSFactory::eb_facets(const EBFArrayBoxFactory * eb_factory) {
    // Container for normal data
    std::unique_ptr<MultiFab> normal = std::unique_ptr<MultiFab>(new MultiFab);
    std::unique_ptr<Vector<Real>> facet_list;

    // Only call the routine for wall collisions if the box has a wall
    if (eb_factory != NULL) {
        dummy->define(eb_ba_refined, mfix_pc->ParticleDistributionMap(amr_lev), 1, eb_grid_pad, MFInfo(), * eb_factory);

        std::array<const MultiCutFab*, AMREX_SPACEDIM> areafrac = eb_factory->getAreaFrac();
        const MultiCutFab * bndrycent = &(eb_factory->getBndryCent());

        int n_facets = 0;

        // We pre-compute the normals
        normal->define(eb_ba_refined, mfix_pc->ParticleDistributionMap(amr_lev), 3, eb_grid_pad);

        for(MFIter mfi(* normal, true); mfi.isValid(); ++mfi) {
            Box tile_box = mfi.growntilebox();
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
        normal->FillBoundary(mfix_pc->Geom(0).periodicity());

        facet_list = std::unique_ptr<Vector<Real>>(new Vector<Real>(6 * n_facets));
        
        int c_facets = 0;
        for(MFIter mfi( * normal, true); mfi.isValid(); ++mfi) {
            Box tile_box = mfi.growntilebox();

            const auto & sfab = dynamic_cast <EBFArrayBox const&>((*dummy)[mfi]);
            const auto & flag = sfab.getEBCellFlagFab();

            const auto & norm_tile = (* normal)[mfi];
            const auto & bcent_tile = (* bndrycent)[mfi];

            // Target: facet_list
            int facet_list_size = facet_list->size();

            eb_as_list(tile_box.loVect(),     tile_box.hiVect(),    & c_facets,
                       flag.dataPtr(),        flag.loVect(),        flag.hiVect(),
                       norm_tile.dataPtr(),   norm_tile.loVect(),   norm_tile.hiVect(),
                       bcent_tile.dataPtr(),  bcent_tile.loVect(),  bcent_tile.hiVect(),
                       facet_list->dataPtr(), & facet_list_size,
                       dx_eb_vect.dataPtr());
        }
    }
    return facet_list;
}


std::unique_ptr<MultiFab> LSFactory::ebis_impfunc(const EBIndexSpace * eb_is) {
    std::unique_ptr<MultiFab> mf_impfunc = std::unique_ptr<MultiFab>(new MultiFab);
    const DistributionMapping & dm = mfix_pc->ParticleDistributionMap(amr_lev);
    mf_impfunc->define(phi_ba_refined, dm, 1, ls_grid_pad);

    for(MFIter mfi(* mf_impfunc, true); mfi.isValid(); ++ mfi)
        eb_is->fillNodeFarrayBoxFromImplicitFunction((* mf_impfunc)[mfi], dx_vect);

    mf_impfunc->FillBoundary(mfix_pc -> Geom(0).periodicity());
    return mf_impfunc;
}


void LSFactory::update(const MultiFab * ls_in) {
    for(MFIter mfi( * ls_phi, true); mfi.isValid(); ++mfi){
        Box tile_box = mfi.growntilebox();

        const auto & ls_in_tile = (* ls_in)[mfi];
        auto & v_tile = (* ls_valid)[mfi];
        auto & ls_tile = (* ls_phi)[mfi];

        update_levelset(tile_box.loVect(),    tile_box.hiVect(),
                        ls_in_tile.dataPtr(), ls_in_tile.loVect(), ls_in_tile.hiVect(),
                        v_tile.dataPtr(),     v_tile.loVect(),     v_tile.hiVect(),
                        ls_tile.dataPtr(),    ls_tile.loVect(),    ls_tile.hiVect(),
                        dx_vect.dataPtr(),    & ls_grid_pad);
    }

    ls_phi -> FillBoundary(mfix_pc -> Geom(0).periodicity());
    ls_valid -> FillBoundary(mfix_pc -> Geom(0).periodicity());
}


void LSFactory::update_ebf(const EBFArrayBoxFactory * eb_factory, const EBIndexSpace * eb_is) {
    if(eb_factory == NULL) return;

    // Generate facets (TODO: in future these can also be provided by user)
    amrex::Print() << "generating facets" << std::endl;
    std::unique_ptr<Vector<Real>> facets = eb_facets(eb_factory);
    amrex::Print() << "done generating facets" << std::endl;
    int len_facets = facets->size();
    // Generate implicit function (used to determine the interior of EB)
    std::unique_ptr<MultiFab> impfunct = ebis_impfunc(eb_is);
    impfunct->FillBoundary(mfix_pc->Geom(0).periodicity());

    // Local MultiFab storing level-set data for this eb_factory
    MultiFab eb_ls;
    iMultiFab eb_valid;

    const DistributionMapping & dm = mfix_pc -> ParticleDistributionMap(amr_lev);
    eb_ls.define(phi_ba_refined, dm, 1, ls_grid_pad);
    eb_valid.define(phi_ba_refined, dm, 1, ls_grid_pad);
    eb_valid.setVal(0);

    amrex::Print() << "filling eb-ls" << std::endl;
    // Fill local MultiFab with eb_factory's level-set data. Note the role of eb_valid:
    //  -> eb_valid = 1 if the corresponding eb_ls location could be projected onto the eb-facets
    //  -> eb_valid = 0 if eb_ls is the fall-back (euclidian) distance to the nearest eb-facet
    for(MFIter mfi(eb_ls, true); mfi.isValid(); ++mfi){
        Box tile_box = mfi.growntilebox();
        const int * lo = tile_box.loVect();
        const int * hi = tile_box.hiVect();

        amrex::Print() << "flag" << std::endl;
        const auto & sfab = dynamic_cast<EBFArrayBox const &>((* dummy)[mfi]);
        const auto & flag = sfab.getEBCellFlagFab();
        amrex::Print() << "done flag" << std::endl;

        // TODO: figure out why this test returns false ...
        //if(flag.getType(tile_box) == FabType::singlevalued){
            amrex::Print() << "v_tile" << std::endl;
            auto & v_tile = eb_valid[mfi];
            amrex::Print() << "ls_tile" << std::endl;
            auto & ls_tile = eb_ls[mfi];
            amrex::Print() << "if_tile" << std::endl;
            const auto & if_tile = (* impfunct)[mfi];

            amrex::Print() << "fill" << std::endl;
            fill_levelset_eb(lo,                hi,
                             facets->dataPtr(), & len_facets,
                             v_tile.dataPtr(),  v_tile.loVect(),  v_tile.hiVect(),
                             ls_tile.dataPtr(), ls_tile.loVect(), ls_tile.hiVect(),
                             dx_vect.dataPtr(), dx_eb_vect.dataPtr());

            amrex::Print() << "validate" << std::endl;
            validate_levelset(lo,                hi,               & ls_grid_ref,
                              if_tile.dataPtr(), if_tile.loVect(), if_tile.hiVect(),
                              v_tile.dataPtr(),  v_tile.loVect(),  v_tile.hiVect(),
                              ls_tile.dataPtr(), ls_tile.loVect(), ls_tile.hiVect());

            amrex::Print() << "mfi iter" << std::endl;
        //}
        amrex::Print() << "done" << std::endl;
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
