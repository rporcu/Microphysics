#include <AMReX_EBFArrayBox.H>
#include <AMReX_Box.H>

#include <mfix.H>
#include <mfix_des_K.H>

//using namespace amrex;

void MFIXParticleContainer::MFIX_PC_ImposeWalls (int lev,
                                                 amrex::EBFArrayBoxFactory* ebfactory,
                                                 const int ls_refinement,
                                                 const amrex::MultiFab* ls_phi)
{
    BL_PROFILE_VAR("MFIXParticleContainer::MFIX_PC_Imposewalls()",mfix_pc_impose_walls);

    /****************************************************************************
     * Get particle EB geometric info
     ***************************************************************************/
    const FabArray<EBCellFlagFab>* flags = &(ebfactory->getMultiEBCellFlagFab());

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const Box& bx = pti.tilebox();
        PairIndex index(pti.index(), pti.LocalTileIndex());

        if ((ebfactory != NULL) and
           ((*flags)[pti].getType(amrex::grow(bx,1)) == FabType::singlevalued))  {


            auto& plev = GetParticles(lev);
            auto& ptile = plev[index];
            auto& aos   = ptile.GetArrayOfStructs();
            ParticleType* pstruct = aos().dataPtr();

            const int nrp = GetParticles(lev)[index].numRealParticles();

            auto& geom = this->Geom(lev);
            const auto dxi = geom.InvCellSizeArray();
            const auto plo = geom.ProbLoArray();
            const auto phiarr = ls_phi->array(pti);

            const amrex::Real En = 0.85;
            const amrex::Real Et = 1.00;

            amrex::ParallelFor(nrp,
                [pstruct,ls_refinement,phiarr,plo,dxi,En,Et]
                AMREX_GPU_DEVICE (int pid) noexcept
                {
                    ParticleType& particle = pstruct[pid];
                    amrex::Real radius = particle.rdata(realData::radius) *
                        std::cbrt(particle.rdata(realData::statwt));

                    amrex::Real ls_value = interp_level_set(particle, ls_refinement, phiarr, plo, dxi);

                    const amrex::Real overlap = radius - ls_value;

                    // The particle is actually touching the wall. Reflect it.
                    if (overlap > 0.)
                    {
                        amrex::RealVect normal(0.);
                        level_set_normal(particle, ls_refinement, normal, phiarr, plo, dxi);

                        // Reflect the particle.
                        particle.pos(0) += overlap*normal[0];
                        particle.pos(1) += overlap*normal[1];
                        particle.pos(2) += overlap*normal[2];

                        // Plane ref point
                        const amrex::Real Nw_Vp = normal[0]*particle.rdata(realData::velx)
                                                + normal[1]*particle.rdata(realData::vely)
                                                + normal[2]*particle.rdata(realData::velz);

                        // Parcel normal velocity
                        const amrex::RealVect Vpn = {Nw_Vp*normal[0],
                                                     Nw_Vp*normal[1],
                                                     Nw_Vp*normal[2]};

                        // Parcel tangential velocity
                        const amrex::RealVect Vpt = {particle.rdata(realData::velx) - Vpn[0],
                                                     particle.rdata(realData::vely) - Vpn[1],
                                                     particle.rdata(realData::velz) - Vpn[2]};

                        // Rebound parcel if moving towards wall.
                        if(Nw_Vp < 0.) {
                            particle.rdata(realData::velx) = -En*Vpn[0] + Et*Vpt[0];
                            particle.rdata(realData::vely) = -En*Vpn[1] + Et*Vpt[1];
                            particle.rdata(realData::velz) = -En*Vpn[2] + Et*Vpt[2];

                        } else {
                            particle.rdata(realData::velx) = Vpn[0] + Et*Vpt[0];
                            particle.rdata(realData::vely) = Vpn[1] + Et*Vpt[1];
                            particle.rdata(realData::velz) = Vpn[2] + Et*Vpt[2];
                        }

                    }
            });

            }
    }
    BL_PROFILE_VAR_STOP(mfix_pc_impose_walls);
}
