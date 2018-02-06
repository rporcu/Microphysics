#include "eb_levelset.H"

#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_EBFArrayBox.H>
//#include <AMReX_EBFabFactory.H>
#include <AMReX_EBIndexSpace.H>


#include "AMReX_BoxIterator.H"

#include <mfix_F.H>


LSFactory::LSFactory(int lev, int ref, const MFIXParticleContainer * pc, 
        const EBIndexSpace * ebis, const Real * dx
        //const EBFArrayBoxFactory * ebfactory
        )
    : amr_lev(lev), ls_grid_refinement(ref), mfix_pc(pc), 
    eb_is(ebis), dx_vect(AMREX_D_DECL(mfix_pc->Geom(amr_lev).CellSize()[0]/ref, 
                                      mfix_pc->Geom(amr_lev).CellSize()[1]/ref, 
                                      mfix_pc->Geom(amr_lev).CellSize()[2]/ref)) 
    // eb_factory(ebfactory) 
{
    //amrex::Print() << "dx_vect" << dx_vect << "\n";

    ls_phi   = std::unique_ptr<MultiFab>(new MultiFab);
    ls_valid = std::unique_ptr<iMultiFab>(new iMultiFab);

    const BoxArray & particle_ba   = mfix_pc -> ParticleBoxArray(lev);
    const BoxArray & phi_ba        = amrex::convert(particle_ba, IntVect{1,1,1});
    const DistributionMapping & dm = mfix_pc -> ParticleDistributionMap(lev);

    BoxArray phi_ba2 = phi_ba;
    phi_ba2.refine(ref);
    BoxArray particle_ba2 = particle_ba;
    particle_ba2.refine(ref);
    
    ls_phi -> define(phi_ba2, dm, 1, ref);
    ls_valid -> define(particle_ba2, dm, 1, 2 * ref);
    ls_valid -> setVal(0);
}


LSFactory::~LSFactory() {
    ls_phi.reset();
    ls_valid.reset();
}

 
void LSFactory::update(MultiFab * dummy){
    //if(eb_factory == NULL)
    //    return;
    //
    //dummy->define(mfix_pc -> ParticleBoxArray(amr_lev), mfix_pc -> ParticleDistributionMap(amr_lev), 
    //              1, 0, MFInfo(), * eb_factory);
    //
    //std::array<const MultiCutFab*, AMREX_SPACEDIM> areafrac = eb_factory->getAreaFrac();
    //
    //for(MFIter mfi( * ls_phi, true); mfi.isValid(); ++ mfi){
    //    Box tile_box = mfi.tilebox();
    //    const int * lo = tile_box.loVect();
    //    const int * hi = tile_box.hiVect();
    //
    //    //const auto & sfab = dynamic_cast<EBFArrayBox const &>((* dummy)[mfi]);
    //    //const auto & flag = sfab.getEBCellFlagFab();
    //
    //    //if(flag.getType(amrex::grow(tile_box,1)) == FabType::singlevalued){
    //        fill_levelset(lo, hi, & ls_grid_refinement,
    //                      (* ls_valid)[mfi].dataPtr(), (* ls_valid)[mfi].loVect(), (* ls_valid)[mfi].hiVect(),
    //                      (* ls_phi)[mfi].dataPtr(), (* ls_phi)[mfi].loVect(), (* ls_phi)[mfi].hiVect(),
    //                      mfix_pc->Geom(amr_lev).CellSize());
    //    //}
    //}

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
             :PolynomialIF(a_polynomial, a_inside){
    
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


  Real PolynomialDF::value(const RealVect& a_point) const
  {
    return value(a_point,m_polynomial);
  }
  ///
  BaseIF* PolynomialDF::newImplicitFunction() const
  {
    PolynomialIF* polynomialPtr = new PolynomialDF(m_polynomial,
                                                   m_inside);

    return static_cast<BaseIF*>(polynomialPtr);
  }


//Real IntersectionDF::value(const RealVect& a_point) const
//{
//  // Maximum of the implicit functions values
//  Real retval;
//
//  retval = 1.0; 
//  std::cout << "pt=" << a_point << std::endl;
//  // Find the maximum value and return it
//  if (m_numFuncs > 0)
//  {
//    retval = m_impFuncs[0]->value(a_point);
//    std::cout << retval << std::endl;
//
//    for (int ifunc = 1; ifunc < m_numFuncs; ifunc++)
//    {
//      Real cur;
//
//      cur = m_impFuncs[ifunc]->value(a_point);
//      std::cout << cur << std::endl;
//      if (cur < retval)
//        retval = cur;
//      
//    }
//  }
//  std::cout << " =>" << retval << std::endl;
//  return retval;
//}
//
//
//BaseIF* IntersectionDF::newImplicitFunction() const
//{
//  IntersectionIF* intersectionPtr = new IntersectionDF(m_impFuncs);
//
//  return static_cast<BaseIF*>(intersectionPtr);
//}
//
//IntersectionDF::IntersectionDF(const Vector<BaseIF *>& a_impFuncs)
//               : IntersectionIF(a_impFuncs)
//{
//  // Number of implicit function in intersection
//  m_numFuncs = a_impFuncs.size();
//
//  // Vector of implicit function pointers
//  m_impFuncs.resize(m_numFuncs);
//
//  // Make copies of the implicit functions
//
//  for (int ifunc = 0; ifunc < m_numFuncs; ifunc++)
//  {
//    if (a_impFuncs[ifunc] == NULL)
//    {
//      m_impFuncs[ifunc] = NULL;
//    }
//    else
//    {
//      m_impFuncs[ifunc] = a_impFuncs[ifunc]->newImplicitFunction();
//    }
//  }
//}
