#include <mfix_des_K.H>

#include <mfix_solids_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_reactions_parms.H>
#include <mfix_bc_list.H>
#include <mfix_bc_parms.H>
#include <mfix_solvers.H>
#include <mfix_calc_cell.H>

using namespace amrex;

void MFIXParticleContainer::mfix_pc_inflow (int lev,
                                            Real dt,
                                            Real time,
                                            EBFArrayBoxFactory* factory)
{
  BL_PROFILE("mfix_dem::mfix_pc_inflow()");


  ParallelDescriptor::Barrier();


  const Real*  dx = Geom(lev).CellSize();
  const Real* idx = Geom(lev).InvCellSize();
  const Real  dx2 = dx[0]*dx[0];

  int total_np = 0;

  constexpr Real tolerance = std::numeric_limits<Real>::epsilon();
  const auto& flags = factory->getMultiEBCellFlagFab();

  const int myProc = ParallelDescriptor::MyProc();

  // Loop over BCs
  for (size_t bcv(0); bcv < BC::bc.size(); ++bcv) {

    // BCs with ep_g < 1 and are EB
    if (BC::bc[bcv].type == BCList::eb &&
        (1.0 - BC::bc[bcv].fluid.volfrac) > tolerance ) {

      const int  has_normal = BC::bc[bcv].eb.has_normal;
      amrex::GpuArray<amrex::Real,3> normal{0.};
      if (has_normal) {
         normal[0] = BC::bc[bcv].eb.normal[0];
         normal[1] = BC::bc[bcv].eb.normal[1];
         normal[2] = BC::bc[bcv].eb.normal[2];
      }

      constexpr Real pad = std::numeric_limits<float>::epsilon();
      const Real normal_tol = BC::bc[bcv].eb.normal_tol;

      const Real norm_tol_lo = Real(-1.) - (normal_tol + pad);
      const Real norm_tol_hi = Real(-1.) + (normal_tol + pad);

      const Box *ic_bx = calc_ic_box(Geom(lev), BC::bc[bcv].region);

      for (MFIter mfi = MakeMFIter(lev,true); mfi.isValid(); ++mfi) {

        //const Box& bx = pti.tilebox();
        const Box& bx = mfi.tilebox();
        FabType t = flags[mfi].getType(bx);

        if (t == FabType::singlevalued && ic_bx->intersects(bx)) {

          Array4<EBCellFlag const> const& flagsfab = flags.const_array(mfi);

          // EB normal and area
          Array4<Real const> const& bnorm = factory->getBndryNormal()[mfi].const_array();
          Array4<Real const> const& barea = factory->getBndryArea().const_array(mfi);
          Array4<Real const> const& bcent = factory->getBndryCent()[mfi].const_array();

          const Box bx_int = bx&(*ic_bx);

          ReduceOps<ReduceOpSum,
                    ReduceOpMin, ReduceOpMin, ReduceOpMin,
                    ReduceOpMax, ReduceOpMax, ReduceOpMax,
                    ReduceOpMin, ReduceOpMin, ReduceOpMin,
                    ReduceOpMax, ReduceOpMax, ReduceOpMax> reduce_op;

          ReduceData<Real, Real, Real, Real, Real, Real, Real,
                    int, int, int, int, int, int> reduce_data(reduce_op);

          using ReduceTuple = typename decltype(reduce_data)::Type;

          reduce_op.eval(bx_int, reduce_data, [flagsfab, barea, bnorm, bcent,
            has_normal, normal, norm_tol_lo, norm_tol_hi]
          AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
          {
            Real area(0.0);
            if(flagsfab(i,j,k).isSingleValued()) {
              area = barea(i,j,k);
              if(has_normal) {
                   const Real dotprod = bnorm(i,j,k,0)*normal[0]
                                      + bnorm(i,j,k,1)*normal[1]
                                      + bnorm(i,j,k,2)*normal[2];
                   area *= (norm_tol_lo <= dotprod) ? Real(1.0) : Real(0.0);
                   area *= (dotprod <= norm_tol_hi) ? Real(1.0) : Real(0.0);
              }
            }
            if (area > 0.99) {
              Real x = static_cast<Real>(i) + 0.5 + bcent(i,j,k,0);
              Real y = static_cast<Real>(j) + 0.5 + bcent(i,j,k,1);
              Real z = static_cast<Real>(k) + 0.5 + bcent(i,j,k,2);
              return {area, x, y, z, x, y, z, i, j, k, i, j, k};
            } else {
              constexpr Real dneR = 1.e32;
              constexpr int  dneI = 16384;

              return {area, dneR, dneR, dneR, -dneR, -dneR, -dneR,
                            dneI, dneI, dneI, -dneI, -dneI, -dneI};
            }
          });

          ReduceTuple htuple = reduce_data.value();

          const Real fab_area = dx2*amrex::get<0>(htuple);

          if (fab_area < tolerance) {
            amrex::AllPrint() << "This fab does not intersect!" << time << "\n";
            continue;
          }

          const GpuArray<Real, 3> plo = Geom(lev).ProbLoArray();

          const RealVect pnt_lo(plo[0] + dx[0]*amrex::get<1>(htuple),
                                plo[1] + dx[1]*amrex::get<2>(htuple),
                                plo[2] + dx[2]*amrex::get<3>(htuple));

          const RealVect pnt_hi(plo[0] + dx[0]*amrex::get<4>(htuple),
                                plo[1] + dx[1]*amrex::get<5>(htuple),
                                plo[2] + dx[2]*amrex::get<6>(htuple));

          const IntVect lo(amrex::get<7>(htuple),
                           amrex::get<8>(htuple),
                           amrex::get<9>(htuple));

          const IntVect hi(amrex::get<10>(htuple),
                           amrex::get<11>(htuple),
                           amrex::get<12>(htuple));

          // Currently we only support bringing in one solid type -- find it
          int type0(0);
          for(; type0 <= BC::bc[bcv].solids.size(); type0++)
            if(BC::bc[bcv].solids[type0].volfrac > tolerance)
              break;

          const SolidsPhase::SOLIDS_t solid = BC::bc[bcv].solids[type0];

          const Real dp   = solid.diameter.mean;
          const Real rhop = solid.density.mean;

          //const Real max_dp = (solid.diameter.distribution == "constant") ? dp : solid.diameter.max;
          //const Real min_dp = (solid.diameter.distribution == "constant") ? dp : solid.diameter.min;
          //const Real std_dp = (solid.diameter.distribution == "constant") ? 0. : solid.diameter.std;

          Real velmag = std::sqrt(solid.velocity[0]*solid.velocity[0]
                                + solid.velocity[1]*solid.velocity[1]
                                + solid.velocity[2]*solid.velocity[2]);

          const Real volflow = (velmag < tolerance) ? solid.volflow :
              solid.volfrac*fab_area*velmag;

          const Real pvol_per_step = volflow*dt + solid.vol_remainder;

          const Real pvol = (M_PI/6.0) * (dp*dp*dp);

          const int pcount = static_cast<int>(amrex::Math::floor(pvol_per_step/pvol));

          BC::bc[bcv].solids[type0].vol_remainder = pvol_per_step - pvol*static_cast<Real>(pcount);

          const int grid_id = mfi.index();
          const int tile_id = mfi.LocalTileIndex();

          // Create a particle container grid and add the particles to it
          auto& ptile = DefineAndReturnParticleTile(lev,mfi);

          const int old_np = ptile.numParticles();
          const int new_np = old_np + pcount;

          total_np += old_np;

          const RealVect eblen(pnt_hi - pnt_lo);

          const int rand_x = (eblen[0] > tolerance) ? 1 : 0;
          const int rand_y = (eblen[1] > tolerance) ? 1 : 0;
          const int rand_z = (eblen[2] > tolerance) ? 1 : 0;

#if 0
          amrex::AllPrint() << "\nvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n"
                            << "We intersect this box! " << time << "\n"
                            << "           has normal: " << has_normal << "\n"
                            << "                   lo: " << lo << "\n"
                            << "                   hi: " << hi << "\n"
                            << "               pnt_lo: " << pnt_lo[0] << "  " << pnt_lo[1] << "  " << pnt_lo[2] << "\n"
                            << "               pnt_hi: " << pnt_hi[0] << "  " << pnt_hi[1] << "  " << pnt_hi[2] << "\n"
                            << "    FAB   bc has area: " << fab_area << "\n"
                            << "    Total bc has area: " << BC::bc[bcv].eb.area << "\n"
                            << "  -- Adding solids  -- " << "\n"
                            << "             diameter: " << dp   << "\n"
                            << "              density  " << rhop << "\n"
                            << " ... particle stats .. " << "\n"
                            << "        initial count: " << old_np << "\n"
                            << "             creating: " << pcount << "\n"
                            << "                final: " << new_np << "\n"
                            << "     remainder volume: " << solid.vol_remainder << "\n"
                            << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
#endif

#if 1
          ptile.resize(new_np);

          // Add runtime-added components
          const int start = AoSrealData::count + SoArealData::count;
          for (int comp(0); comp < m_runtimeRealData.count; ++comp)
            ptile.push_back_real(start+comp, pcount, 0.);

          auto& particles = ptile.GetArrayOfStructs();
          ParticleType* pstruct = particles().dataPtr();

          auto& soa = ptile.GetStructOfArrays();
          auto p_realarray = soa.realarray();
          auto p_intarray = soa.intarray();

          auto ptile_data = ptile.getParticleTileData();

          const int nextID = ParticleType::NextID();
#endif

          // This loop is written for CPU only because it's only expected to
          // operate on a few particles.

          int added(0);
          int fails(0);
          for(; added < pcount && fails < 10000;) {

            const int ip = old_np + added;

            ParticleType& p = pstruct[ip];

            Real rad = 0.5*dp;
            Real mass = pvol * rhop;
            Real omoi = 2.5/(mass * rad*rad);

            Real xt = pnt_lo[0] + (rand_x ? eblen[0]*amrex::Random() : 0.);
            Real yt = pnt_lo[1] + (rand_y ? eblen[1]*amrex::Random() : 0.);
            Real zt = pnt_lo[2] + (rand_z ? eblen[2]*amrex::Random() : 0.);

            int i = static_cast<int>(amrex::Math::floor((xt - plo[0])*idx[0]));
            int j = static_cast<int>(amrex::Math::floor((yt - plo[1])*idx[1]));
            int k = static_cast<int>(amrex::Math::floor((zt - plo[2])*idx[2]));

            if(!flagsfab(i,j,k).isSingleValued()) {
              fails++;
              continue;
            }


            if (barea(i,j,k) < 0.95) {
              fails++;
              continue;
            }


            if(has_normal) {

              const Real dotprod = bnorm(i,j,k,0)*normal[0]
                                 + bnorm(i,j,k,1)*normal[1]
                                 + bnorm(i,j,k,2)*normal[2];

              if ((norm_tol_lo > dotprod) || (dotprod > norm_tol_hi)) {
                fails++;
                continue;
              }
            }

            // physical location of the eb centroid
            const Real bcx = (plo[0] + dx[0]*(static_cast<Real>(i) + 0.5 + bcent(i,j,k,0)));
            const Real bcy = (plo[1] + dx[1]*(static_cast<Real>(j) + 0.5 + bcent(i,j,k,1)));
            const Real bcz = (plo[2] + dx[2]*(static_cast<Real>(k) + 0.5 + bcent(i,j,k,2)));

            // distance of particle from eb face
            Real offset = bnorm(i,j,k,0)*(xt - bcx)
                        + bnorm(i,j,k,1)*(yt - bcy)
                        + bnorm(i,j,k,2)*(zt - bcz);

            Real shift = 1.21*(rad - offset);

            xt -= shift*bnorm(i,j,k,0);
            yt -= shift*bnorm(i,j,k,1);
            zt -= shift*bnorm(i,j,k,2);

            //amrex::AllPrint() << xt << " " << yt << " " << zt << " " << offset << "\n";

            added++;

#if 1
            p.pos(0) = xt;
            p.pos(1) = yt;
            p.pos(2) = zt;

            p_realarray[SoArealData::velx][ip] = -bnorm(i,j,k,0)*velmag;
            p_realarray[SoArealData::vely][ip] = -bnorm(i,j,k,1)*velmag;
            p_realarray[SoArealData::velz][ip] = -bnorm(i,j,k,2)*velmag;

            // Set other particle properties
            p_intarray[SoAintData::phase][ip] = type0;
            p_intarray[SoAintData::state][ip] = 1;
            p_realarray[SoArealData::volume][ip] = pvol;
            p_realarray[SoArealData::density][ip] = rhop;
            p_realarray[SoArealData::mass][ip] = mass;
            p_realarray[SoArealData::oneOverI][ip] = omoi;
            p_realarray[SoArealData::radius][ip] = rad;
            p_realarray[SoArealData::omegax][ip] = 0.;
            p_realarray[SoArealData::omegay][ip] = 0.;
            p_realarray[SoArealData::omegaz][ip] = 0.;
            p_realarray[SoArealData::statwt][ip] = 1.;
            p_realarray[SoArealData::dragcoeff][ip] = 0.;
            p_realarray[SoArealData::dragx][ip] = 0.;
            p_realarray[SoArealData::dragy][ip] = 0.;
            p_realarray[SoArealData::dragz][ip] = 0.;
            p_realarray[SoArealData::cp_s][ip] = 0.;
            p_realarray[SoArealData::temperature][ip] = 0.;
            p_realarray[SoArealData::convection][ip] = 0.;

            // Set id and cpu for this particle
            p.id()  = nextID + added;
            p.cpu() = myProc;
#if 0
            int start_idx = idx_X_sn;
            int end_idx   = idx_mass_sn_txfr;

            // Runtime added variables -- species mass fractions
            for (int idx(start_idx); idx < end_idx; ++idx) {
              // Copy data from particle to replicated one
              ptile_data.m_runtime_rdata[idx][ip] = ptile_data.m_runtime_rdata[idx][index];
            }

            start_idx = end_idx;
            end_idx   = idx_vel_s_txfr;

            // Runtime added variables -- species mass txfr rates
            for (int idx(start_idx); idx < end_idx; ++idx) {
              // Copy data from particle to replicated one
              ptile_data.m_runtime_rdata[idx][ip] = ptile_data.m_runtime_rdata[idx][index];
            }

            start_idx = end_idx;
            end_idx   = idx_h_s_txfr;

            // Runtime added variables -- species momentum txfr rate
            for (int idx(start_idx); idx < end_idx; ++idx) {
              // Copy data from particle to replicated one
              ptile_data.m_runtime_rdata[idx][ip] = ptile_data.m_runtime_rdata[idx][index];
            }

            start_idx = end_idx;
            end_idx   = idx_count;

            // Runtime added variables -- species energy txfr rate
            for (int idx(start_idx); idx < end_idx; ++idx) {
              // Copy data from particle to replicated one
              ptile_data.m_runtime_rdata[idx][ip] = ptile_data.m_runtime_rdata[idx][index];
            }
#endif
#endif
          }
          ParticleType::NextID(nextID + pcount);
          Gpu::synchronize();

          total_np += added;

        } // ic_bx intersects MFIter tile box

        else
          amrex::AllPrint() << "This box does not intersect!" << time << "\n";

      } // MFIter loop
    } // if eb and ep_g < 1.0
  } // loop over bcs



  Redistribute();

  ParallelDescriptor::ReduceIntSum(total_np,ParallelDescriptor::IOProcessorNumber());
  amrex::Print() << "Total number of particles: " <<
    time << " " << total_np << std::endl;


  ParallelDescriptor::Barrier();



}
