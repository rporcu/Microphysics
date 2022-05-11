#include <mfix_des_K.H>

#include <mfix_solids_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_reactions_parms.H>
#include <mfix_bc_list.H>
#include <mfix_bc_parms.H>
#include <mfix_solvers.H>
#include <mfix_calc_cell.H>

using namespace amrex;

static constexpr int is_constant = 0;
static constexpr int is_uniform  = 1;
static constexpr int is_normal   = 2;

namespace{

AMREX_INLINE
int get_distribution (const std::string dist_str )
{
  if (dist_str.compare("normal") == 0) {
    return is_normal;
  } else if (dist_str.compare("uniform") == 0) {
    return is_uniform;
  } else {
    return is_constant;
  }
};


AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real
rand_prop (const amrex::Real var_mean,
           const amrex::Real var_min,
           const amrex::Real var_max,
           const amrex::Real var_std,
           const int dist_type,
           RandomEngine const& engine) noexcept
{
  if (dist_type == is_normal) {
    Real var(var_min - var_max);

    while( var < var_min || var > var_max) {
      Real x(0.0);
      Real y(0.0);
      Real w(1.1);
      while(w > 1. || amrex::Math::abs(w-1.) < std::numeric_limits<Real>::epsilon()) {
        x = 2.*amrex::Random(engine) - 1.;
        y = 2.*amrex::Random(engine) - 1.;
        w = x*x + y*y;
      }
      w = std::sqrt((-2.*std::log(w)) / w);
      var = var_mean + var_std*x*w;
    }
    return var;
  } else if (dist_type == is_uniform) {
    return var_min + (var_max - var_min)*amrex::Random(engine);
  } else {
    return var_mean;
  }
};
}//end namespace

void MFIXParticleContainer::mfix_pc_inflow (int lev,
                                            Real dt,
                                            Real time,
                                            const int adv_enthalpy,
                                            EBFArrayBoxFactory* factory)
{
  BL_PROFILE("mfix_dem::mfix_pc_inflow()");

  const int solve_species = solids.solve_species;

  const int nspecies_s = solids.nspecies;
  const int solid_is_a_mixture = solids.is_a_mixture;

  auto& solids_parms = *solids.parameters;

  const auto  dx_array = Geom(lev).CellSizeArray();
  const auto idx_array = Geom(lev).InvCellSizeArray();
  const auto plo_array = Geom(lev).ProbLoArray();

  const RealVect  dx( dx_array[0],  dx_array[1],  dx_array[2]);
  const RealVect idx(idx_array[0], idx_array[1], idx_array[2]);
  const RealVect plo(plo_array[0], plo_array[1], plo_array[2]);

  const Real  dx2 = dx[0]*dx[0];

  int total_np = 0;

  constexpr Real tolerance = std::numeric_limits<Real>::epsilon();
  const auto& flags = factory->getMultiEBCellFlagFab();

  const int MyProc = ParallelDescriptor::MyProc();

  // Loop over BCs
  for (size_t bcv(0); bcv < BC::bc.size(); ++bcv) {

    // BCs with ep_g < 1 and are EB
    if (BC::bc[bcv].type == BCList::eb &&
        (1.0 - BC::bc[bcv].fluid.volfrac) > tolerance ) {

      // This shouldn't be needed but "just in case"
      if(BC::bc[bcv].solids.size() == 0)
        continue;

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

        const Box& tile_box = mfi.tilebox();
        FabType t = flags[mfi].getType(tile_box);

        if (t == FabType::singlevalued && ic_bx->intersects(tile_box)) {

          Array4<EBCellFlag const> const& flagsfab = flags.const_array(mfi);

          // EB normal and area
          Array4<Real const> const& bnorm = factory->getBndryNormal()[mfi].const_array();
          Array4<Real const> const& barea = factory->getBndryArea().const_array(mfi);
          Array4<Real const> const& bcent = factory->getBndryCent()[mfi].const_array();

          const Box bx_int = tile_box&(*ic_bx);

          Real fab_area(0.);
          IntVect bblo, bbhi;

          { // Limiting the scope of the reduce operation
            ReduceOps<ReduceOpSum,
                      ReduceOpMin, ReduceOpMin, ReduceOpMin,
                      ReduceOpMax, ReduceOpMax, ReduceOpMax> reduce_op;

            ReduceData<Real, int, int, int, int, int, int> reduce_data(reduce_op);
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
              return {area, i, j, k, i, j, k};
            });

            ReduceTuple htuple = reduce_data.value();

            fab_area = dx2*amrex::get<0>(htuple);

            bblo[0] = amrex::get<1>(htuple);
            bblo[1] = amrex::get<2>(htuple);
            bblo[2] = amrex::get<3>(htuple);

            bbhi[0] = amrex::get<4>(htuple);
            bbhi[1] = amrex::get<5>(htuple);
            bbhi[2] = amrex::get<6>(htuple);

          } // end reduce scope bound

          if (fab_area < tolerance) {
            //amrex::AllPrint() << "This fab does not intersect!" << time << "\n";
            continue;
          }

          const amrex::Box bbx(bblo, bbhi);

          const RealBox rbbx (bbx,Geom(lev).CellSize(),Geom(lev).ProbLo());

          const int isize = nspecies_s+2; // Tp and cps
          Gpu::HostVector<Real> h_bc_inputs(isize);

          for(int phase(0); phase < BC::bc[bcv].solids.size(); phase++) {
            if(BC::bc[bcv].solids[phase].volfrac > tolerance) {

              const SOLIDS_t solid = BC::bc[bcv].solids[phase];

              Real Tp(0.0), cp_s(0.0);
              if (adv_enthalpy) {
                if(!solid_is_a_mixture) {
                  cp_s = solids_parms.calc_cp_s<RunOn::Host>(phase,Tp);

                } else {
                  for (int n(0); n<nspecies_s; n++) {
                    const Real Xs_n = solid.species[n].mass_fraction;
                    cp_s += (Xs_n * solids_parms.calc_cp_sn<RunOn::Host>(Tp,n));
                  }
                }
              }
              h_bc_inputs[0] = cp_s;
              h_bc_inputs[1] = Tp;

              for (int n(0); n < nspecies_s; n++) {
                h_bc_inputs[2+n] = solve_species ? solid.species[n].mass_fraction : 0.0;
              }

              Gpu::AsyncArray<Real> d_bc_inputs(h_bc_inputs.data(), h_bc_inputs.size());

              Real* p_bc_inputs = (adv_enthalpy || solve_species) ?  d_bc_inputs.data() : nullptr;


              const int dist_dp = get_distribution(solid.diameter.distribution);

              const Real mean_dp = solid.diameter.mean;
              const Real  max_dp = (dist_dp == is_constant) ? mean_dp : solid.diameter.max;
              const Real  min_dp = (dist_dp == is_constant) ? mean_dp : solid.diameter.min;
              const Real  std_dp = (dist_dp == is_constant) ? 0. : solid.diameter.std;

              const int dist_rhop = get_distribution(solid.density.distribution);

              const Real mean_rhop = solid.density.mean;
              const Real  max_rhop = (dist_rhop == is_constant) ? mean_rhop : solid.density.max;
              const Real  min_rhop = (dist_rhop == is_constant) ? mean_rhop : solid.density.min;
              const Real  std_rhop = (dist_rhop == is_constant) ? 0.        : solid.density.std;

              const Real act_area = solid.volfrac*fab_area;

              const Real volflow = (solid.velmag < tolerance) ? solid.volflow :
                  act_area*solid.velmag;

              const Real velmag = amrex::max(solid.velmag, volflow/act_area);

              const Real pvol_per_step = volflow*dt + solid.vol_remainder;

              const Real mean_pvol = (M_PI/6.0) * (mean_dp*mean_dp*mean_dp);

              const int pcount = static_cast<int>(amrex::Math::floor(pvol_per_step/mean_pvol));

              const int grid_id = mfi.index();
              const int tile_id = mfi.LocalTileIndex();

              // Flags to skip random numbers for directions with zero length
              const int rand_x = (rbbx.length(0) > tolerance) ? 1 : 0;
              const int rand_y = (rbbx.length(1) > tolerance) ? 1 : 0;
              const int rand_z = (rbbx.length(2) > tolerance) ? 1 : 0;

              auto& ptile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
              const int old_np = ptile.numParticles();

              // Increase the particle container to hold the added particles.
              const int new_np = old_np + pcount;
              ptile.resize(new_np);

              // Add components for each of the runtime variables
              const int start = AoSrealData::count + SoArealData::count;
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

              amrex::ParallelForRNG(pcount, [pstruct, ptile_data, NextID, MyProc, old_np, bbx,
                p_real,p_int, phase, plo, dx, idx, rbbx,
                has_normal, normal, norm_tol_lo, norm_tol_hi, flagsfab, bnorm, bcent,
                rand_x, rand_y, rand_z, velmag,
                dist_dp,   mean_dp,   max_dp,   min_dp,   std_dp,
                dist_rhop, mean_rhop, max_rhop, min_rhop, std_rhop,
                adv_enthalpy, solve_species, p_bc_inputs, nspecies_s, idx_X_sn ]
                AMREX_GPU_DEVICE (int pid, amrex::RandomEngine const& engine) noexcept
              {
                const int ip = pid + old_np;

                ParticleType& p = pstruct[ip];

                p.id()  = NextID + pid;
                p.cpu() = MyProc;

                const Real rad = (dist_dp == is_constant) ?  0.5*mean_dp :
                  0.5*rand_prop(mean_dp, min_dp, max_dp, std_dp, dist_dp, engine);

                bool valid = false;
                while ( !valid ) {

                  // Grab a random point within the bounding box
                  p.pos(0) = rbbx.lo(0) + (rand_x ? rbbx.length(0)*amrex::Random(engine) : 0.);
                  p.pos(1) = rbbx.lo(1) + (rand_y ? rbbx.length(1)*amrex::Random(engine) : 0.);
                  p.pos(2) = rbbx.lo(2) + (rand_z ? rbbx.length(2)*amrex::Random(engine) : 0.);

                  // compute the cell index containing the particle center
                  int i = static_cast<int>(amrex::Math::floor((p.pos(0) - plo[0])*idx[0]));
                  int j = static_cast<int>(amrex::Math::floor((p.pos(1) - plo[1])*idx[1]));
                  int k = static_cast<int>(amrex::Math::floor((p.pos(2) - plo[2])*idx[2]));

                  if ( bbx.contains(IntVect(i,j,k))) {
                    if(flagsfab(i,j,k).isSingleValued()) {
                      valid = true;
                      if(has_normal) {

                        const Real dotprod = bnorm(i,j,k,0)*normal[0]
                                           + bnorm(i,j,k,1)*normal[1]
                                           + bnorm(i,j,k,2)*normal[2];

                        valid = ((norm_tol_lo <= dotprod) && (dotprod <= norm_tol_hi));
                      }
                    }
                  }

                  if (valid) {

                    // physical location of the eb centroid
                    const Real bcx = (plo[0] + dx[0]*(static_cast<Real>(i) + 0.5 + bcent(i,j,k,0)));
                    const Real bcy = (plo[1] + dx[1]*(static_cast<Real>(j) + 0.5 + bcent(i,j,k,1)));
                    const Real bcz = (plo[2] + dx[2]*(static_cast<Real>(k) + 0.5 + bcent(i,j,k,2)));

                    // distance of particle from eb face
                    Real offset = bnorm(i,j,k,0)*(p.pos(0) - bcx)
                                + bnorm(i,j,k,1)*(p.pos(1) - bcy)
                                + bnorm(i,j,k,2)*(p.pos(2) - bcz);

                    Real shift = -1.01*(rad - offset);

                    p.pos(0) -= shift*bnorm(i,j,k,0);
                    p.pos(1) -= shift*bnorm(i,j,k,1);
                    p.pos(2) -= shift*bnorm(i,j,k,2);

                    p_real[SoArealData::velx][ip] = -bnorm(i,j,k,0)*velmag;
                    p_real[SoArealData::vely][ip] = -bnorm(i,j,k,1)*velmag;
                    p_real[SoArealData::velz][ip] = -bnorm(i,j,k,2)*velmag;
                  }

                }

                p_int[SoAintData::phase][ip] = phase+1;
                p_int[SoAintData::state][ip] = 0;

                const Real rhop = (dist_rhop == is_constant) ? mean_rhop :
                  rand_prop(mean_rhop, min_rhop, max_rhop, std_rhop, dist_rhop, engine);

                const Real pvol = (4.0/3.0) * M_PI * (rad*rad*rad);
                const Real mass = pvol * rhop;
                const Real omoi = 2.5/(mass * rad*rad);

                p_real[SoArealData::volume][ip] = pvol;
                p_real[SoArealData::density][ip] = rhop;
                p_real[SoArealData::mass][ip] = mass;
                p_real[SoArealData::oneOverI][ip] = omoi;
                p_real[SoArealData::radius][ip] = rad;
                p_real[SoArealData::omegax][ip] = 0.;
                p_real[SoArealData::omegay][ip] = 0.;
                p_real[SoArealData::omegaz][ip] = 0.;
                p_real[SoArealData::statwt][ip] = 1.;
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
              });

              // Update the particles NextID
              ParticleType::NextID(NextID+pcount);

              Gpu::synchronize();

              Real total_new_vol(mean_pvol*static_cast<Real>(pcount));
              if (dist_dp != is_constant) {
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

              // Store the particle volume underflow/overflow so we can account
              // for it on the next time step.
              BC::bc[bcv].solids[phase].vol_remainder = pvol_per_step - total_new_vol;

            } // solids phase has volfrac > 0
          } // loop over solids phases
        } // ic_bx intersects MFIter tile box
      } // MFIter loop
    } // if bc_type is eb and ep_g < 1.0
  } // loop over bcs

  Gpu::synchronize();

  Redistribute(0, 0, 0, 1);

  ParallelDescriptor::ReduceIntSum(total_np,ParallelDescriptor::IOProcessorNumber());
  amrex::Print() << "Total number of particles: " <<
    time << " " << total_np << std::endl;
}
