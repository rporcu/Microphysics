#include "eb_levelset.H"
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>

#include <mfix_F.H>

double LSFactory::value_wall(const double pos[3]){
    return R2 - ( (pos[0] - xy_centre[0]) * (pos[0] - xy_centre[0]) +
                  (pos[1] - xy_centre[1]) * (pos[1] - xy_centre[1]) );
}

double LSFactory::value_bottom(const double pos[3]){
    return pos[2];
}


LSFactory::LSFactory(int lev, const MFIXParticleContainer * pc, const EBFArrayBoxFactory * ebfactory)
    : amr_lev(lev), mfix_pc(pc), eb_factory(ebfactory) 
{
    ls_phi   = std::unique_ptr<MultiFab>(new MultiFab);
    ls_valid = std::unique_ptr<iMultiFab>(new iMultiFab);

    const BoxArray & particle_ba   = mfix_pc -> ParticleBoxArray(lev);
    const BoxArray & phi_ba        = amrex::convert(particle_ba, IntVect{1,1,1});
    const DistributionMapping & dm = mfix_pc -> ParticleDistributionMap(lev);

    ls_phi -> define(phi_ba, dm, 1, 1);
    ls_valid -> define(particle_ba, dm, 1, 1);
    ls_valid -> setVal(0);
}

LSFactory::~LSFactory() {
    ls_phi.reset();
    ls_valid.reset();
}

void LSFactory::update(MultiFab * dummy){
    if(eb_factory == NULL)
        return;

    dummy->define(mfix_pc -> ParticleBoxArray(amr_lev), mfix_pc -> ParticleDistributionMap(amr_lev), 
                  1, 0, MFInfo(), * eb_factory);
   
    std::array<const MultiCutFab*, AMREX_SPACEDIM> areafrac = eb_factory->getAreaFrac();

    for(MFIter mfi( * ls_phi, true); mfi.isValid(); ++ mfi){
        Box tile_box = mfi.tilebox();
        const int * lo = tile_box.loVect();
        const int * hi = tile_box.hiVect();

        const auto & sfab = dynamic_cast<EBFArrayBox const &>((* dummy)[mfi]);
        const auto & flag = sfab.getEBCellFlagFab();

        if(flag.getType(amrex::grow(tile_box,1)) == FabType::singlevalued){
            compute_levelset(lo, hi,
                             flag.dataPtr(), flag.loVect(), flag.hiVect(),
                             (* ls_valid)[mfi].dataPtr(), (* ls_valid)[mfi].loVect(), (* ls_valid)[mfi].hiVect(),
                             (* ls_phi)[mfi].dataPtr(), (* ls_phi)[mfi].loVect(), (* ls_phi)[mfi].hiVect() );
        }
    }

    ls_phi -> FillBoundary(mfix_pc -> Geom(0).periodicity());
    ls_valid -> FillBoundary(mfix_pc -> Geom(0).periodicity());
}
