#include <AMReX_EBFArrayBox.H>
#include <AMReX_Box.H>

#include <mfix.H>
#include <mfix_des_K.H>

//using namespace amrex;

void MFIXParticleContainer::MFIX_PC_ImposeWalls (int lev,
                                                 amrex::EBFArrayBoxFactory* ebfactory,
                                                 const int ls_refinement,
                                                 const amrex::MultiFab* ls_phi,
                                                 amrex::MultiFab* cost,
                                                 std::string& knapsack_weight_type)
{
    BL_PROFILE("MFIXParticleContainer::MFIX_PC_Imposewalls()");

    /****************************************************************************
     * Get particle EB geometric info
     ***************************************************************************/
    const FabArray<EBCellFlagFab>* flags = &(ebfactory->getMultiEBCellFlagFab());

    auto& geom = this->Geom(lev);
    const Real* dx = Geom(lev).CellSize();
    const auto dxi = geom.InvCellSizeArray();
    const auto p_lo = geom.ProbLoArray();
    const auto p_hi = geom.ProbHiArray();

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        // Timer used for load-balancing
        amrex::Real wt = ParallelDescriptor::second();

        const Box& bx = pti.tilebox();
        PairIndex index(pti.index(), pti.LocalTileIndex());

        // Determine if this particle tile actually has any walls
        bool has_wall = false;

        if ((ebfactory != NULL) and
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


        auto& plev = GetParticles(lev);
        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        ParticleType* pstruct = aos().dataPtr();

        const int nrp = GetParticles(lev)[index].numRealParticles();

        if( has_wall ) {

            const auto& phiarr = ls_phi->array(pti);

            const amrex::Real En = 0.85;
            const amrex::Real Et = 1.00;

            amrex::ParallelFor(nrp,
              [pstruct,ls_refinement,phiarr,p_lo,dxi,En,Et]
              AMREX_GPU_DEVICE (int pid) noexcept
              {
                  ParticleType& particle = pstruct[pid];
                  amrex::Real radius = particle.rdata(realData::radius) *
                      std::cbrt(particle.rdata(realData::statwt));

                  amrex::Real ls_value = interp_level_set(particle, ls_refinement, phiarr, p_lo, dxi);

                  const amrex::Real overlap = radius - ls_value;

                  // The particle is actually touching the wall. Reflect it.
                  if (overlap > 0.)
                  {
                      amrex::RealVect normal(0.);
                      level_set_normal(particle, ls_refinement, normal, phiarr, p_lo, dxi);

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
            amrex::Gpu::synchronize();

        }

        amrex::Gpu::synchronize();

        if (ebfactory != NULL) {

          const auto& flagsarr = (*flags)[pti].array();

          amrex::ParallelFor(nrp,
            [pstruct,flagsarr,p_lo,p_hi,dxi]
            AMREX_GPU_DEVICE (int pid) noexcept
            {
                ParticleType& p = pstruct[pid];


                if(p_lo[0] < p.pos(0) and p.pos(0) < p_hi[0] and
                   p_lo[1] < p.pos(1) and p.pos(1) < p_hi[1] and
                   p_lo[2] < p.pos(2) and p.pos(2) < p_hi[2]) {

                  amrex::Real x = (p.pos(0) - p_lo[0]) * dxi[0];
                  amrex::Real y = (p.pos(1) - p_lo[1]) * dxi[1];
                  amrex::Real z = (p.pos(2) - p_lo[2]) * dxi[2];

                  int i = static_cast<int>(amrex::Math::floor(x));
                  int j = static_cast<int>(amrex::Math::floor(y));
                  int k = static_cast<int>(amrex::Math::floor(z));

                  if(flagsarr(i,j,k).isCovered())
                    p.id() = -1;

                } else {
                  p.id() = -1;

                }


            });
          amrex::Gpu::synchronize();
        }

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

        amrex::Gpu::synchronize();
    }

}
