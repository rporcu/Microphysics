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

    ls_phi   = std::unique_ptr<MultiFab>(new MultiFab);
    ls_valid = std::unique_ptr<iMultiFab>(new iMultiFab);

    const BoxArray & particle_ba   = mfix_pc -> ParticleBoxArray(lev);
    const BoxArray & phi_ba        = amrex::convert(particle_ba, IntVect{1,1,1});
    const DistributionMapping & dm = mfix_pc -> ParticleDistributionMap(lev);

    phi_ba_refined = phi_ba;
    phi_ba_refined.refine(ref);
    particle_ba_refined = particle_ba;
    particle_ba_refined.refine(ref);

    ls_phi -> define(phi_ba_refined, dm, 1, ref);
    ls_valid -> define(particle_ba_refined, dm, 1, 2 * ref);
    //ls_valid -> setVal(0);

    for(MFIter mfi( * ls_phi, true); mfi.isValid(); ++mfi){
        Box tile_box = mfi.tilebox();
        const int * lo = tile_box.loVect();
        const int * hi = tile_box.hiVect();

        init_levelset(lo, hi, & ls_grid_refinement,
                      (* ls_valid)[mfi].dataPtr(), (* ls_valid)[mfi].loVect(), (* ls_valid)[mfi].hiVect(),
                      (* ls_phi)[mfi].dataPtr(), (* ls_phi)[mfi].loVect(), (* ls_phi)[mfi].hiVect());
    }

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
            const int* lo = tile_box.loVect();
            const int* hi = tile_box.hiVect();
            
            const auto& sfab = dynamic_cast <EBFArrayBox const&>((*dummy)[mfi]);
            const auto& flag = sfab.getEBCellFlagFab();

            count_eb_facets(lo, hi, flag.dataPtr(), flag.loVect(), flag.hiVect(), & n_facets);
            
            if (flag.getType(amrex::grow(tile_box,1)) == FabType::singlevalued) {
               BL_PROFILE_VAR("compute_normals()", compute_normals);
               compute_normals(lo, hi, flag.dataPtr(), flag.loVect(), flag.hiVect(),
                               (*normal)[mfi].dataPtr(),
                               (*normal)[mfi].loVect(), (*normal)[mfi].hiVect(),
                               (*areafrac[0])[mfi].dataPtr(),
                               (*areafrac[0])[mfi].loVect(), (*areafrac[0])[mfi].hiVect(),
                               (*areafrac[1])[mfi].dataPtr(),
                               (*areafrac[1])[mfi].loVect(), (*areafrac[1])[mfi].hiVect(),
                               (*areafrac[2])[mfi].dataPtr(),
                               (*areafrac[2])[mfi].loVect(), (*areafrac[2])[mfi].hiVect());
               BL_PROFILE_VAR_STOP(compute_normals);
            }

        }
        normal->FillBoundary(mfix_pc->Geom(0).periodicity());

        amrex::Print() << "n_facets=" << n_facets << std::endl;

        Box dom_list(IntVect{0,0,0}, IntVect{n_facets,0,0});
        facet_list = std::unique_ptr<FArrayBox>(new FArrayBox(dom_list, 6));

        int c_facets = 0;
        for(MFIter mfi( * normal, true); mfi.isValid(); ++mfi) {
            Box tile_box = mfi.tilebox();
            const int* lo = tile_box.loVect();
            const int* hi = tile_box.hiVect();

            const auto& sfab = dynamic_cast <EBFArrayBox const&>((*dummy)[mfi]);
            const auto& flag = sfab.getEBCellFlagFab();


            eb_as_list(lo, hi, & ls_grid_refinement,
                       flag.dataPtr(), flag.loVect(), flag.hiVect(),
                       (*normal)[mfi].dataPtr(),
                       (*normal)[mfi].loVect(), (*normal)[mfi].hiVect(),
                       (* bndrycent)[mfi].dataPtr(), 
                       (* bndrycent)[mfi].loVect(), (* bndrycent)[mfi].hiVect(),
                       facet_list->dataPtr(),
                       facet_list->loVect(), facet_list->hiVect(),
                       mfix_pc->Geom(amr_lev).CellSize(), & c_facets);
        }
        amrex::Print() << "c_facets=" << c_facets << std::endl;

    }
    return facet_list;
}



void LSFactory::update_ebf(const EBFArrayBoxFactory * eb_factory) {
    if(eb_factory == NULL)
        return;
    
    //dummy->define(mfix_pc -> ParticleBoxArray(amr_lev), mfix_pc -> ParticleDistributionMap(amr_lev),
    //              1, 0, MFInfo(), * eb_factory);
    //std::array<const MultiCutFab*, AMREX_SPACEDIM> areafrac = eb_factory->getAreaFrac();

    std::unique_ptr<FArrayBox> facets = eb_facets(eb_factory);

    for(MFIter mfi( * ls_phi, true); mfi.isValid(); ++ mfi){
        Box tile_box = mfi.tilebox();
        const int * lo = tile_box.loVect();
        const int * hi = tile_box.hiVect();
    
        //const auto & sfab = dynamic_cast<EBFArrayBox const &>((* dummy)[mfi]);
        //const auto & flag = sfab.getEBCellFlagFab();
    
        //if(flag.getType(amrex::grow(tile_box,1)) == FabType::singlevalued){
            //const MultiCutFab * bndrycent = &(eb_factory->getBndryCent());

            fill_levelset(lo, hi, & ls_grid_refinement,
                          facets->dataPtr(), facets->loVect(), facets->hiVect(),
                          (* ls_valid)[mfi].dataPtr(), (* ls_valid)[mfi].loVect(), (* ls_valid)[mfi].hiVect(),
                          (* ls_phi)[mfi].dataPtr(), (* ls_phi)[mfi].loVect(), (* ls_phi)[mfi].hiVect(),
                          mfix_pc->Geom(amr_lev).CellSize());
        //}
    }
    
    ls_phi -> FillBoundary(mfix_pc -> Geom(0).periodicity());
    ls_valid -> FillBoundary(mfix_pc -> Geom(0).periodicity());
}


void LSFactory::update_ebis(const EBIndexSpace * eb_is){
    for(MFIter mfi( * ls_phi, true); mfi.isValid(); ++ mfi){
        eb_is->fillNodeFarrayBoxFromImplicitFunction(( * ls_phi)[mfi], dx_vect);
    }

    for(MFIter mfi( * ls_phi, true); mfi.isValid(); ++ mfi){
        FArrayBox & a_fab = (* ls_phi)[mfi];
        for(BoxIterator bit(a_fab.box()); bit.ok(); ++bit)
                a_fab(bit(), 0) = - a_fab(bit(), 0);
    }

    ls_phi -> FillBoundary(mfix_pc -> Geom(0).periodicity());
    ls_valid -> FillBoundary(mfix_pc -> Geom(0).periodicity());
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
