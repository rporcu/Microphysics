#include <mfix_pc.H>
#include <mfix_des_K.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_reactions.H>
#include <mfix_bc.H>
#include <mfix_solvers.H>
#include <mfix_calc_cell.H>

using namespace amrex;

namespace{

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
rand_prop (const amrex::Real var_mean,
           const amrex::Real var_min,
           const amrex::Real var_max,
           const amrex::Real var_std,
           const bool is_normal,
           RandomEngine const& engine) noexcept
{
  if (is_normal) {
    Real var = amrex::RandomNormal(var_mean, var_std, engine);
    while (var_min > var || var > var_max){
      var = amrex::RandomNormal(var_mean, var_std, engine);
    }
    return var;
  } else {
    return var_min + (var_max - var_min)*amrex::Random(engine);
  }
}
}//end namespace

void
MFIXParticleContainer::
mfix_pc_inflow (int const lev,
                int const is_dem,
                int const is_pic,
                Real dt,
                const int adv_enthalpy,
                EBFArrayBoxFactory* factory)
{
  BL_PROFILE("mfix_pc_inflow()");

  const int solve_species = solids.solve_species();
  const int solve_reactions = reactions.solve();

  const int nspecies_s = solids.nspecies();
  const int solid_is_a_mixture = solids.isMixture();

  const auto& solids_parms = solids.parameters();

  const auto  dx_array = Geom(lev).CellSizeArray();
  const auto idx_array = Geom(lev).InvCellSizeArray();
  const auto plo_array = Geom(lev).ProbLoArray();
  const auto phi_array = Geom(lev).ProbHiArray();

  const RealVect  dx( dx_array[0],  dx_array[1],  dx_array[2]);
  const RealVect idx(idx_array[0], idx_array[1], idx_array[2]);
  const RealVect plo(plo_array[0], plo_array[1], plo_array[2]);
  const RealVect phi(phi_array[0], phi_array[1], phi_array[2]);

  int total_np = 0;

  constexpr Real tolerance = std::numeric_limits<Real>::epsilon();
  const auto& flags = factory->getMultiEBCellFlagFab();

  const int MyProc = ParallelDescriptor::MyProc();

  Real const ep_cp = m_pic.ep_cp();
  Real const inv_ep_cp = (is_dem ? 1. : 1.0/ep_cp);

  const int isize = nspecies_s+2; // Tp and cps
  Gpu::HostVector<Real> h_bc_inputs(isize);

  // Loop over BCs
  for (int bcv(0); bcv < m_boundary_conditions.bc().size(); ++bcv) {

    int const is_mi = m_boundary_conditions.bc(bcv).type == BCList::minf;
    int const is_eb = m_boundary_conditions.bc(bcv).type == BCList::eb;

    // Better safe than sorry!
    if ( !is_mi && !is_eb ) { continue; }

    // BCs with ep_g < 1 and have solids
    if ( ((1.0 - m_boundary_conditions.bc(bcv).fluid.volfrac) < tolerance) &&
         (m_boundary_conditions.bc(bcv).solids.size() == 0) ) { continue; }

    int const has_normal( is_mi ? 0 : m_boundary_conditions.bc(bcv).eb.has_normal );

    amrex::GpuArray<amrex::Real,3> normal{0.};

    if (has_normal) {
       normal[0] = m_boundary_conditions.bc(bcv).eb.normal[0];
       normal[1] = m_boundary_conditions.bc(bcv).eb.normal[1];
       normal[2] = m_boundary_conditions.bc(bcv).eb.normal[2];
    }

    constexpr Real pad = std::numeric_limits<float>::epsilon();

    const Real normal_tol = m_boundary_conditions.bc(bcv).eb.normal_tol;

    const Real norm_tol_lo = Real(-1.) - (normal_tol + pad);
    const Real norm_tol_hi = Real(-1.) + (normal_tol + pad);

    // 0:x-lo, 1:x-hi, 2:y-lo, 3:y-hi, 4:z-lo, 5:z-hi
    int const face = m_boundary_conditions.get_dir(bcv);

    int const dir_lohi(face%2); // lo=0, hi=1
    int const dir((face-dir_lohi)/((int)2)); // dir=0,1,2 (x,y,z)

    Real const bc_area = m_boundary_conditions.get_bc_area(bcv);

    if (bc_area < tolerance) {
      Print() << "\n"
              << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
              << "               Particle flow: Area for BC " << bcv << "is zero!\n"
              << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
              << "\n";
      continue;
    }

    Box bc_bx = is_mi ?
      calc_bc_box(Geom(lev), m_boundary_conditions.bc(bcv).region, face):
      calc_bc_box(Geom(lev), m_boundary_conditions.bc(bcv).region,   -1);

    for(int lcs(0); lcs < m_boundary_conditions.bc(bcv).solids.size(); lcs++) {

      SOLIDS_t& bc_solid = m_boundary_conditions.bc(bcv).solids[lcs];

      if(bc_solid.volfrac < tolerance) { continue; }

      const int phase = bc_solid.phase;
      const int state = (is_dem) ? 0 : 1;

      const int dp_is_constant = bc_solid.diameter.is_constant();
      const int dp_is_normal   = bc_solid.diameter.is_normal();

      const Real mean_dp = bc_solid.diameter.get_mean();
      const Real  max_dp = bc_solid.diameter.get_max();
      const Real  min_dp = bc_solid.diameter.get_min();
      const Real  std_dp = bc_solid.diameter.get_stddev();

      const int rhop_is_constant = bc_solid.density.is_constant();
      const int rhop_is_normal   = bc_solid.density.is_normal();

      const Real mean_rhop = bc_solid.density.get_mean();
      const Real  max_rhop = bc_solid.density.get_max();
      const Real  min_rhop = bc_solid.density.get_min();
      const Real  std_rhop = bc_solid.density.get_stddev();

      const Real statwt = (is_pic) ? bc_solid.statwt : 1.0;

      Real Tp = bc_solid.temperature;
      Real cp_s(0.0);

      if (adv_enthalpy) {
        if(!solid_is_a_mixture) {
          cp_s = solids_parms.calc_cp_s<RunOn::Host>(Tp);

        } else {
          for (int n(0); n<nspecies_s; n++) {
            const Real Xs_n = bc_solid.species[n].mass_fraction;
            cp_s += (Xs_n * solids_parms.calc_cp_sn<RunOn::Host>(Tp,n));
          }
        }
      }
      h_bc_inputs[0] = cp_s;
      h_bc_inputs[1] = Tp;

      for (int n(0); n < nspecies_s; n++) {
        h_bc_inputs[2+n] = solve_species ? bc_solid.species[n].mass_fraction : 0.0;
      }

      Gpu::AsyncArray<Real> d_bc_inputs(h_bc_inputs.data(), h_bc_inputs.size());

      Real* p_bc_inputs = (adv_enthalpy || solve_species) ?  d_bc_inputs.data() : nullptr;

      Real const flow_area = (bc_solid.volfrac*bc_area);

      Real const bc_volflow = (bc_solid.volflow > tolerance) ?
          bc_solid.volflow : flow_area*bc_solid.velmag;

      Real const bc_velmag  = (bc_solid.velmag > tolerance) ?
          bc_solid.velmag : (bc_volflow / flow_area);

      RealVect const bc_velvec(bc_solid.velocity[0],bc_solid.velocity[1],bc_solid.velocity[2]);

      Real vol_remainder(bc_solid.vol_remainder);

      // Make sure that the computed (or provided) flowrate is positive.
      if (bc_volflow < tolerance) { continue; }

      for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {

        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();

        std::pair<int,int> index(grid_id, tile_id);

        Real const fab_area = m_boundary_conditions.get_bc_fab_area(bcv, index);

        if (fab_area < tolerance) {
          //amrex::AllPrint() << "This fab does not intersect!" << time << "\n";
          continue;
        }

        const Box& tile_box = mfi.tilebox();
        FabType t = flags[mfi].getType(tile_box);

        if (t == FabType::covered || !bc_bx.intersects(tile_box) ) { continue; }

        // Scale volflow by the fraction of the area contained in this fab.
        const Real fab_volflow = bc_volflow*(fab_area/bc_area);

        const Real pvol_per_step = fab_volflow*dt + vol_remainder;

        const Real mean_pvol = statwt * (M_PI/6.0) * (mean_dp*mean_dp*mean_dp);

        const int pcount = static_cast<int>(amrex::Math::floor(pvol_per_step/mean_pvol));

        if (pcount == 0) {
          vol_remainder = pvol_per_step;

        } else {

          auto& ptile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
          const int old_np = ptile.numParticles();

          // Increase the particle container to hold the added particles.
          const int new_np = old_np + pcount;
          ptile.resize(new_np);

          // Add components for each of the runtime variables
          const int start = SoArealData::count;
          for (int comp(0); comp < m_runtimeRealData.count; ++comp)
            ptile.push_back_real(start+comp, pcount, 0.);

          auto& ptile_aos  = ptile.GetArrayOfStructs();
          auto& ptile_soa  = ptile.GetStructOfArrays();
          auto  ptile_data = ptile.getParticleTileData();

          ParticleType* pstruct = ptile_aos().dataPtr();

          auto p_real = ptile_soa.realarray();
          auto p_int = ptile_soa.intarray();

          const int NextID = ParticleType::NextID();

          total_np += new_np;

          const int idx_X_sn = m_runtimeRealData.X_sn;
          const int idx_mass_txfr = m_runtimeRealData.mass_txfr;
          const int idx_vel_txfr = m_runtimeRealData.vel_txfr;
          const int idx_h_txfr = m_runtimeRealData.h_txfr;

          // DEM particles go OUTSIDE the EB surface, PIC parcels
          // are place inside the domain.
          Real const shift_dir   = ( is_dem ? -1.01 : 1.01);
          Real const shift_lohi  = ( (dir_lohi==0) ? 1.0 : -1.0);

          bool const regular_fab = (t == FabType::regular);
          Real const eff_rad_scale = (is_dem ? 1.0 : std::cbrt(statwt * inv_ep_cp));

          const Box bbx = tile_box&(bc_bx);

          const RealBox rbbx (bbx,Geom(lev).CellSize(),Geom(lev).ProbLo());

          Array4<EBCellFlag const> const& flagsfab = flags.const_array(mfi);

          // EB normal and face centroid
          Array4<Real const> const& bnorm = (t != FabType::singlevalued) ? Array4<Real const>{} :
                                            factory->getBndryNormal()[mfi].const_array();

          Array4<Real const> const& bcent = (t != FabType::singlevalued) ? Array4<Real const>{} :
                                             factory->getBndryCent()[mfi].const_array();

          amrex::ParallelForRNG(pcount, [pstruct, ptile_data, NextID, MyProc, old_np, bbx,
              p_real, p_int, phase, state, is_dem, is_pic, statwt, plo, phi, dx, idx, rbbx,
              regular_fab, is_mi, dir_lohi, dir, shift_dir, shift_lohi, eff_rad_scale,
              bc_velmag, bc_velvec, has_normal, normal, norm_tol_lo, norm_tol_hi, flagsfab, bnorm, bcent,
              mean_dp,   max_dp,   min_dp,   std_dp,     dp_is_constant,   dp_is_normal,
              mean_rhop, max_rhop, min_rhop, std_rhop, rhop_is_constant, rhop_is_normal,
              adv_enthalpy, solve_species, p_bc_inputs, nspecies_s, idx_X_sn,
              idx_mass_txfr, idx_vel_txfr, idx_h_txfr, solve_reactions]
            AMREX_GPU_DEVICE (int pid, amrex::RandomEngine const& engine) noexcept
          {
            const int ip = pid + old_np;

            ParticleType& p = pstruct[ip];

            p.id()  = NextID + pid;
            p.cpu() = MyProc;

            const Real rad = (dp_is_constant) ?  0.5*mean_dp :
              0.5*rand_prop(mean_dp, min_dp, max_dp, std_dp, dp_is_normal, engine);

            // Effective radius of the parcel same as radius for DEM
            Real const eff_rad = rad * eff_rad_scale;

            bool valid = false;
            while ( !valid ) {

              // Grab a random point within the bounding box
              p.pos(0) = rbbx.lo(0) + rbbx.length(0)*amrex::Random(engine);
              p.pos(1) = rbbx.lo(1) + rbbx.length(1)*amrex::Random(engine);
              p.pos(2) = rbbx.lo(2) + rbbx.length(2)*amrex::Random(engine);

              if (is_mi) {
                  Real const face_offset = (dir_lohi==0) ? plo[dir] : phi[dir];
                  p.pos(dir) = face_offset + shift_lohi*eff_rad;
              }

              // compute the cell index containing the particle center
              int i = static_cast<int>(amrex::Math::floor((p.pos(0) - plo[0])*idx[0]));
              int j = static_cast<int>(amrex::Math::floor((p.pos(1) - plo[1])*idx[1]));
              int k = static_cast<int>(amrex::Math::floor((p.pos(2) - plo[2])*idx[2]));

              if ( bbx.contains(IntVect(i,j,k))) {
                if ( regular_fab ) {
                  valid = true;
                } else if(is_mi && flagsfab(i,j,k).isRegular()) {
                    valid = true;
                } else if(flagsfab(i,j,k).isSingleValued()) {

                  if(is_mi) {

                    // physical location of the eb centroid
                    const Real bcx = (plo[0] + dx[0]*(static_cast<Real>(i) + 0.5 + bcent(i,j,k,0)));
                    const Real bcy = (plo[1] + dx[1]*(static_cast<Real>(j) + 0.5 + bcent(i,j,k,1)));
                    const Real bcz = (plo[2] + dx[2]*(static_cast<Real>(k) + 0.5 + bcent(i,j,k,2)));

                    // distance of particle from eb face
                    Real offset = bnorm(i,j,k,0)*(p.pos(0) - bcx)
                                + bnorm(i,j,k,1)*(p.pos(1) - bcy)
                                + bnorm(i,j,k,2)*(p.pos(2) - bcz);

                    valid = (offset > eff_rad);

                  } else {
                    valid = true;
                    if (has_normal) {
                      const Real dotprod = bnorm(i,j,k,0)*normal[0]
                                         + bnorm(i,j,k,1)*normal[1]
                                         + bnorm(i,j,k,2)*normal[2];

                      valid = ((norm_tol_lo <= dotprod) && (dotprod <= norm_tol_hi));
                    }
                  } // is not mi (is eb)
                } // regular fab
              } // bbx contains (i,j,k)

              if (valid) {

                if (is_mi) { // mass inflow (pic only)

                  //IntVect ijk(i,j,k);

                    //static_cast<Real>(ijk[dir]+dir_lohi)*dx[dir];

                  p_real[SoArealData::velx][ip] = bc_velvec[0];
                  p_real[SoArealData::vely][ip] = bc_velvec[1];
                  p_real[SoArealData::velz][ip] = bc_velvec[2];

                } else { // eb-inflow (dem or pic)

                  // physical location of the eb centroid
                  const Real bcx = (plo[0] + dx[0]*(static_cast<Real>(i) + 0.5 + bcent(i,j,k,0)));
                  const Real bcy = (plo[1] + dx[1]*(static_cast<Real>(j) + 0.5 + bcent(i,j,k,1)));
                  const Real bcz = (plo[2] + dx[2]*(static_cast<Real>(k) + 0.5 + bcent(i,j,k,2)));

                  // distance of particle from eb face
                  Real offset = bnorm(i,j,k,0)*(p.pos(0) - bcx)
                              + bnorm(i,j,k,1)*(p.pos(1) - bcy)
                              + bnorm(i,j,k,2)*(p.pos(2) - bcz);

                  Real shift = shift_dir*(eff_rad - offset);

                  p.pos(0) -= shift*bnorm(i,j,k,0);
                  p.pos(1) -= shift*bnorm(i,j,k,1);
                  p.pos(2) -= shift*bnorm(i,j,k,2);

                  p_real[SoArealData::velx][ip] = -bnorm(i,j,k,0)*bc_velmag;
                  p_real[SoArealData::vely][ip] = -bnorm(i,j,k,1)*bc_velmag;
                  p_real[SoArealData::velz][ip] = -bnorm(i,j,k,2)*bc_velmag;
                }
              } // is valid
            } // while not valid

            p_int[SoAintData::phase][ip] = phase;
            p_int[SoAintData::state][ip] = state;

            const Real rhop = (rhop_is_constant) ? mean_rhop :
              rand_prop(mean_rhop, min_rhop, max_rhop, std_rhop,
                rhop_is_normal, engine);

            const Real pvol = (4.0/3.0) * M_PI * (rad*rad*rad);
            const Real mass = pvol * rhop;
            const Real omoi = (is_pic ? 0. : (2.5/(mass * rad*rad)));

            p_real[SoArealData::volume][ip] = pvol;
            p_real[SoArealData::density][ip] = rhop;
            p_real[SoArealData::mass][ip] = mass;
            p_real[SoArealData::oneOverI][ip] = omoi;
            p_real[SoArealData::radius][ip] = rad;
            p_real[SoArealData::omegax][ip] = 0.;
            p_real[SoArealData::omegay][ip] = 0.;
            p_real[SoArealData::omegaz][ip] = 0.;
            p_real[SoArealData::statwt][ip] = statwt;
            p_real[SoArealData::dragcoeff][ip] = 0.;
            p_real[SoArealData::dragx][ip] = 0.;
            p_real[SoArealData::dragy][ip] = 0.;
            p_real[SoArealData::dragz][ip] = 0.;
            p_real[SoArealData::cp_s][ip] = 0.;
            p_real[SoArealData::temperature][ip] = 0.;
            p_real[SoArealData::convection][ip] = 0.;

            if(adv_enthalpy) {
              p_real[SoArealData::cp_s][ip] = p_bc_inputs[0];
              p_real[SoArealData::temperature][ip] = p_bc_inputs[1];
            }

            if(solve_species) {
              for (int n(0); n < nspecies_s; n++) {
                ptile_data.m_runtime_rdata[idx_X_sn+n][ip] = p_bc_inputs[2+n];
              }
            }

            if (solve_reactions) {
              for (int n_s(0); n_s < nspecies_s; ++n_s)
                ptile_data.m_runtime_rdata[idx_mass_txfr+n_s][ip] = 0;

              ptile_data.m_runtime_rdata[idx_vel_txfr+0][ip] = 0;
              ptile_data.m_runtime_rdata[idx_vel_txfr+1][ip] = 0;
              ptile_data.m_runtime_rdata[idx_vel_txfr+2][ip] = 0;

              ptile_data.m_runtime_rdata[idx_h_txfr][ip] = 0;
            }

          });

          // Update the particles NextID
          ParticleType::NextID(NextID+pcount);

          Gpu::synchronize();

          Real total_new_vol(mean_pvol*static_cast<Real>(pcount));
          if (!dp_is_constant) {
            // Reduce sum over new particles to take into account
            // that the diameters may be different.
            ReduceOps<ReduceOpSum> reduce_op;
            ReduceData<Real> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;

            reduce_op.eval(pcount, reduce_data, [old_np, p_real]
              AMREX_GPU_DEVICE (int pid) -> ReduceTuple
            {
              Real pvol = p_real[SoArealData::volume][pid + old_np];
              return pvol;
            });

            ReduceTuple tuple = reduce_data.value();
            total_new_vol = amrex::get<0>(tuple);
          }

          vol_remainder = pvol_per_step - total_new_vol;

        } // pcount == 0
      } // MFIter loop

      // Store the particle volume underflow/overflow so we can account
      // for it on the next time step.
      bc_solid.vol_remainder = vol_remainder;

    } // loop over solids phases

  } // loop over bcs

  Gpu::synchronize();

  Redistribute(0, 0, 0, 1);

  //ParallelDescriptor::ReduceIntSum(total_np,ParallelDescriptor::IOProcessorNumber());
  //amrex::Print() << "Total number of particles: " <<
  //  time << " " << total_np << std::endl;
}
