#include <AMReX_EBFArrayBox.H>
#include <AMReX_Box.H>

#include <mfix.H>
#include <mfix_des_K.H>

using namespace amrex;

void MFIXParticleContainer::MFIX_PC_ImposeWalls (int lev,
                                                 EBFArrayBoxFactory* ebfactory,
                                                 const int ls_refinement,
                                                 const MultiFab* ls_phi,
                                                 MultiFab* cost,
                                                 std::string& knapsack_weight_type)
{
    BL_PROFILE_VAR("MFIXParticleContainer::MFIX_PC_Imposewalls()",mfix_pc_impose_walls);

    /****************************************************************************
     * Get particle EB geometric info
     ***************************************************************************/
    const FabArray<EBCellFlagFab>* flags = &(ebfactory->getMultiEBCellFlagFab());

    const Real* dx = Geom(lev).CellSize();

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        // Timer used for load-balancing
        Real wt = ParallelDescriptor::second();

        const Box& bx = pti.tilebox();
        PairIndex index(pti.index(), pti.LocalTileIndex());

        // Determine if this particle tile actually has any walls
        bool has_wall = false;

        if ((ebfactory != NULL) &&
           ((*flags)[pti].getType(amrex::grow(bx,1)) == FabType::singlevalued))  {

          has_wall = true;

        } else {
          // We need this test for the case of an inflow boundary:
          // inflow does not appear in the EBFactory but
          // the particles see it as a wall

          // Create the nodal refined box based on the current particle tile
          Box refined_box(amrex::convert(amrex::refine(bx,ls_refinement), IntVect{1,1,1}));

          // Set tol to 1/2 dx
          Real tol = amrex::min(dx[0], amrex::min(dx[1], dx[2])) / 2;

          Real ls_min_over_box = ((*ls_phi)[pti]).min<RunOn::Gpu>(refined_box,0);

          if (ls_min_over_box < tol) has_wall = true;

        }

        if( has_wall ) {

            auto& plev = GetParticles(lev);
            auto& ptile = plev[index];
            auto& aos   = ptile.GetArrayOfStructs();
            ParticleType* pstruct = aos().dataPtr();

            auto& soa = ptile.GetStructOfArrays();
            auto p_realarray = soa.realarray();

            const int nrp = GetParticles(lev)[index].numRealParticles();

            auto& geom = this->Geom(lev);
            const auto dxi = geom.InvCellSizeArray();
            const auto plo = geom.ProbLoArray();
            const auto& phiarr = ls_phi->array(pti);

            const Real En = 0.85;
            const Real Et = 1.00;

            amrex::ParallelFor(nrp,
                [pstruct,p_realarray,ls_refinement,phiarr,plo,dxi,En,Et]
                AMREX_GPU_DEVICE (int pid) noexcept
                {
                    ParticleType& particle = pstruct[pid];
                    Real radius = p_realarray[SoArealData::radius][pid] *
                        std::cbrt(p_realarray[SoArealData::statwt][pid]);

                    Real ls_value = interp_level_set(particle.pos(), ls_refinement, phiarr, plo, dxi);

                    const Real overlap = radius - ls_value;

                    // The particle is actually touching the wall. Reflect it.
                    if (overlap > 0.)
                    {
                        RealVect normal(0.);
                        level_set_normal(particle.pos(), ls_refinement, normal, phiarr, plo, dxi);

                        // Reflect the particle.
                        particle.pos(0) += overlap*normal[0];
                        particle.pos(1) += overlap*normal[1];
                        particle.pos(2) += overlap*normal[2];

                        // Plane ref point
                        const Real Nw_Vp = normal[0]*p_realarray[SoArealData::velx][pid]
                                                + normal[1]*p_realarray[SoArealData::vely][pid]
                                                + normal[2]*p_realarray[SoArealData::velz][pid];

                        // Parcel normal velocity
                        const RealVect Vpn = {Nw_Vp*normal[0],
                                                     Nw_Vp*normal[1],
                                                     Nw_Vp*normal[2]};

                        // Parcel tangential velocity
                        const RealVect Vpt = {p_realarray[SoArealData::velx][pid] - Vpn[0],
                                                     p_realarray[SoArealData::vely][pid] - Vpn[1],
                                                     p_realarray[SoArealData::velz][pid] - Vpn[2]};

                        // Rebound parcel if moving towards wall.
                        if(Nw_Vp < 0.) {
                            p_realarray[SoArealData::velx][pid] = -En*Vpn[0] + Et*Vpt[0];
                            p_realarray[SoArealData::vely][pid] = -En*Vpn[1] + Et*Vpt[1];
                            p_realarray[SoArealData::velz][pid] = -En*Vpn[2] + Et*Vpt[2];

                        } else {
                            p_realarray[SoArealData::velx][pid] = Vpn[0] + Et*Vpt[0];
                            p_realarray[SoArealData::vely][pid] = Vpn[1] + Et*Vpt[1];
                            p_realarray[SoArealData::velz][pid] = Vpn[2] + Et*Vpt[2];
                        }

                    }
            });


          /********************************************************************
           * Update runtime cost (used in load-balancing)                     *
           *******************************************************************/
          if (cost)
          {
            // Runtime cost is either (weighted by tile box size):
            //   * time spent
            //   * number of particles
            const Box& tbx = pti.tilebox();
            if (knapsack_weight_type == "RunTimeCosts")
            {
              wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
            }
            else if (knapsack_weight_type == "NumParticles")
            {
              wt = nrp / tbx.d_numPts();
            }
            (*cost)[pti].plus<RunOn::Device>(wt, tbx);
          }

        }

    }
    BL_PROFILE_VAR_STOP(mfix_pc_impose_walls);
}
