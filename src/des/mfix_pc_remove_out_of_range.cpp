#include <AMReX.H>

#include <mfix.H>
#include <mfix_des_K.H>
#include <mfix_pc.H>
#include <mfix_dem.H>
#include <mfix_calc_cell.H>
#include <mfix_utils.H>

using namespace amrex;

void MFIXParticleContainer::RemoveOutOfRange (int lev,
                                              const EBFArrayBoxFactory * ebfactory,
                                              const MultiFab * ls_phi,
                                              int ls_refinement)
{
    // Only call the routine for wall collisions if we actually have walls
    if (ebfactory != NULL) {

        const Real* cell_size  = Geom(lev).CellSize();

        const RealVect dx(cell_size[0], cell_size[1], cell_size[2]);
        const GpuArray<Real,3> plo = Geom(lev).ProbLoArray();

        const Geometry& gm  = Geom(0);
        const auto      p_lo = gm.ProbLoArray();
        const auto      dxi = gm.InvCellSizeArray();

        // This holds the mesh spacing of the level set, which may be finer than
        // the local mesh spacing
        const GpuArray<Real,3> dx_ls{dx[0]/ls_refinement, dx[1]/ls_refinement, dx[2]/ls_refinement};

        const FabArray<EBCellFlagFab>* flags = &(ebfactory->getMultiEBCellFlagFab());

        const Real inv_ep_cp = (m_pic.solve()) ? 1.0/m_pic.ep_cp() : 1.0;

        for (MFIXParIter pti(* this, lev); pti.isValid(); ++pti)
        {
            // Real particles

            const Box & bx = pti.tilebox();

            // Remove particles outside of or touching the walls
            if ((*flags)[pti].getType(bx) != FabType::regular)
            {
                auto& aos = pti.GetArrayOfStructs();
                ParticleType* pstruct = aos().dataPtr();

                auto& soa = pti.GetStructOfArrays();
                auto p_realarray = soa.realarray();

                const int np = pti.numParticles();

                if ((*flags)[pti].getType(bx) == FabType::covered)
                {
                    amrex::ParallelFor(np, [pstruct] AMREX_GPU_DEVICE (int ip) noexcept
                    {
                      ParticleType& p = pstruct[ip];
                      p.id() = -1;
                    });
                }
                else
                {
                    const auto& flag_fab =  flags->array(pti);
                    const auto&  phi_arr = ls_phi->array(pti);

                    amrex::ParallelFor(np, [pstruct,p_realarray,plo,dx,flag_fab,inv_ep_cp,
                    dx_ls,phi_arr, ls_refinement, p_lo, dxi, cg_dem=m_dem.cg_dem()]
                    AMREX_GPU_DEVICE (int ip) noexcept
                    {
                        ParticleType& p = pstruct[ip];

                        const Real* plo_ptr = plo.data();
                        int icell = static_cast<int>(amrex::Math::floor( ( p.pos(0) - plo_ptr[0] ) / dx[0] ));
                        int jcell = static_cast<int>(amrex::Math::floor( ( p.pos(1) - plo_ptr[1] ) / dx[1] ));
                        int kcell = static_cast<int>(amrex::Math::floor( ( p.pos(2) - plo_ptr[2] ) / dx[2] ));

                        if (flag_fab(icell,jcell,kcell).isCovered())
                        {
                            p.id() = -1;
                        }
                        else
                        { // Interpolates level-set from nodal phi to position pos


                          RealVect pos(p.pos());
                          Real ls_value = interp_level_set(pos, ls_refinement, phi_arr, p_lo, dxi);

                          Real radius = p_realarray[SoArealData::radius][ip] *
                            std::cbrt(p_realarray[SoArealData::statwt][ip] * inv_ep_cp);

                          if (cg_dem) {
                            radius = radius/std::cbrt(p_realarray[SoArealData::statwt][ip]);
                          }

                          const Real overlap = radius - ls_value;

                          if (overlap > 0.0) {
                            p.id() = -1;
                          }
                        }
                    });
                }
            }
        }

        Redistribute();

        long fin_np = 0;
        for (MFIXParIter pti(* this, lev); pti.isValid(); ++pti) {
            long np = pti.numParticles();
            fin_np += np;
        }

        ParallelDescriptor::ReduceLongSum(fin_np);

        if (mfix::m_run_type == RunType::PIC2DEM)
          amrex::Print() << "Final number of particles on level "
                         << lev << ": " << fin_np << std::endl;

        m_total_numparticle = fin_np;
    }

    if (mfix::m_run_type != RunType::PIC2DEM)
      ReportParticleGenerationStats(lev);

}



void MFIXParticleContainer::
ReportParticleGenerationStats (int lev)
{

  // Store particle count totals by IC region
  std::vector<long> total_np(m_initial_conditions.ic().size(), 0);

  constexpr Real tolerance = std::numeric_limits<Real>::epsilon();

  for (int icv(0); icv < m_initial_conditions.ic().size(); icv++) {

    if (Math::abs(m_initial_conditions.ic(icv).fluid.volfrac-1.) > tolerance) {

      const RealBox* ic_region = m_initial_conditions.ic(icv).region;

      // Reduce sum operation for np, Tp
      ReduceOps<ReduceOpSum> reduce_op;
      ReduceData<long> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

        Box bx = pti.tilebox();
        RealBox tile_region(bx, Geom(lev).CellSize(), Geom(lev).ProbLo());

        if (tile_region.intersects( *ic_region )) {

          const int np         = NumberOfParticles(pti);
          const AoS &particles = pti.GetArrayOfStructs();
          const ParticleType* pstruct = particles().dataPtr();

          const RealBox ic_rbx(ic_region->lo(), ic_region->hi());

          reduce_op.eval(np, reduce_data, [pstruct,ic_rbx]
            AMREX_GPU_DEVICE (int p_id) -> ReduceTuple
          {
            constexpr long l_0 = static_cast<long>(0);
            constexpr long l_1 = static_cast<long>(1);

            return {ic_rbx.contains(pstruct[p_id].pos()) ? l_1 : l_0 };
          });

        } // tile_region intersects ic_region
      } // MFIXParIter loop

      ReduceTuple host_tuple = reduce_data.value();
      total_np[icv] += amrex::get<0>(host_tuple);

    } // ep_g < 1
  } // loop over ICs

  ParallelDescriptor::ReduceLongSum(total_np.data(), m_initial_conditions.ic().size());

  amrex::ParmParse pp("ic");

  std::vector<std::string> input_regions;
  pp.queryarr("regions", input_regions);


  amrex::Print() << "\n";
  amrex::Print() << "  IC Region                Generated           Removed         Remaining\n";
  amrex::Print() << "  ****************  ****************  ****************  ****************\n";

  long sum_np_gen(0);
  long sum_np_dff(0);
  long sum_np_tot(0);

  for (int icv(0); icv < m_initial_conditions.ic().size(); icv++) {

    const long np_gen = m_initial_conditions.get_particle_count(icv);
    const long np_dff = np_gen - total_np[icv];

    sum_np_gen += np_gen;
    sum_np_dff += np_dff;
    sum_np_tot += total_np[icv];

    amrex::Print() << "  " << std::left  << std::setw(16) << input_regions[icv]
                   << "  " << std::right << std::setw(16) << MfixIO::FormatWithCommas(np_gen)
                   << "  " << std::right << std::setw(16) << MfixIO::FormatWithCommas(np_dff)
                   << "  " << std::right << std::setw(16) << MfixIO::FormatWithCommas(total_np[icv])
                   << "\n";
  }
  amrex::Print() << "  ----------------  ----------------  ----------------  ----------------\n";

  amrex::Print() << "  " << std::left  << std::setw(16) << "Total"
                 << "  " << std::right << std::setw(16) << MfixIO::FormatWithCommas(sum_np_gen)
                 << "  " << std::right << std::setw(16) << MfixIO::FormatWithCommas(sum_np_dff)
                 << "  " << std::right << std::setw(16) << MfixIO::FormatWithCommas(sum_np_tot)
                 << "\n";

  amrex::Print() << "  ****************  ****************  ****************  ****************\n";

}
