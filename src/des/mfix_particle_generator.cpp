#include <mfix_particle_generator.H>
#include <AMReX_AmrParGDB.H>
#include <mfix_calc_cell.H>

#include <limits>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <mfix_ic_parms.H>
#include <mfix_solids_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  SUBROUTINE: GENERATE_PARTICLE_CONFIG                                !
//                                                                      !
//  Purpose: Generate particle configuration based on maximum particle  !
//           radius and filling from top to bottom within specified     !
//           bounds                                                     !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

using namespace amrex;

ParticlesGenerator::ParticlesGenerator ()
  : m_rdata(0)
  , m_idata(0)
{}

ParticlesGenerator::~ParticlesGenerator ()
{}

void
ParticlesGenerator::generate (int& pc,
                              const IntVect& lo,
                              const IntVect& hi,
                              const Real dx,
                              const Real dy,
                              const Real dz,
                              const amrex::GpuArray<Real, 3>& plo)
{
  const int init_pc(pc);

  const Real tolerance = std::numeric_limits<Real>::epsilon();

  int np(0);
  int icv(0);
  int type(0);

  // Get the IC index
  int icv0(0);

  for(; icv0 < IC::ic.size(); icv0++)
  {
    if( amrex::Math::abs(IC::ic[icv0].fluid.volfrac -1.0) > tolerance)
    {
      // Get the solids type index
      int type0(0);

      for(; type0 <= IC::ic[icv0].solids.size(); type0++)
        if(IC::ic[icv0].solids[type0].volfrac > tolerance)
          break;

      std::string ic_pack_type_str = IC::ic[icv0].packing;

      int np0(0);

      int cube_base(-1);
      const size_t cube_pos = ic_pack_type_str.find("-cube");
      if(cube_pos != std::string::npos) {
        const size_t cube_len = ic_pack_type_str.length();
        std::string str_base = ic_pack_type_str.erase(cube_pos,cube_len-cube_pos);
        cube_base = std::stoi( str_base );
      }


      if(DEM::solve){
        if(ic_pack_type_str.compare("hcp") == 0)
          hex_close_pack(icv0, type0, lo, hi, np0, pc, dx, dy, dz, plo);
        else if(ic_pack_type_str.compare("random") == 0)
          random_fill_dem(icv0, type0, lo, hi, np0, pc, dx, dy, dz, plo, false);
        else if(ic_pack_type_str.compare("pseudo_random") == 0)
          random_fill_dem(icv0, type0, lo, hi, np0, pc, dx, dy, dz, plo, true);
        else if(ic_pack_type_str.compare("oneper") == 0)
          n_cube_per_fill(icv0, type0, 1, lo, hi, np0, pc, dx, dy, dz, plo);
        else if(ic_pack_type_str.compare("eightper") == 0)
          n_cube_per_fill(icv0, type0, 2, lo, hi, np0, pc, dx, dy, dz, plo);
        else if(cube_base > 0)
          n_cube_per_fill(icv0, type0, cube_base, lo, hi, np0, pc, dx, dy, dz, plo);
        else
        {
          amrex::Print() << "Unknown particle generator fill type" << std::endl;
          exit(1000);
        }
      }

      if(PIC::solve){
        if(ic_pack_type_str.compare("random") == 0)
          random_fill_pic(icv0, type0, lo, hi, np0, pc, dx, dy, dz, plo, false);
        else if(ic_pack_type_str.compare("pseudo_random") == 0)
          random_fill_pic(icv0, type0, lo, hi, np0, pc, dx, dy, dz, plo, true);
        else if(ic_pack_type_str.compare("oneper") == 0)
          n_cube_per_fill(icv0, type0, 1, lo, hi, np0, pc, dx, dy, dz, plo);
        else if(ic_pack_type_str.compare("eightper") == 0)
          n_cube_per_fill(icv0, type0, 2, lo, hi, np0, pc, dx, dy, dz, plo);
        else if(cube_base > 0)
          n_cube_per_fill(icv0, type0, cube_base, lo, hi, np0, pc, dx, dy, dz, plo);
        else
        {
          amrex::Print() << "Unknown particle generator fill type" << std::endl;
          exit(1000);
        }
      }

      // HACK -- the original code assumed that only one IC region would have
      // particles. This saves the IC region and and type to use later.
      if(np0 > 0)
      {
        type = type0;
        icv  = icv0;
        np  += np0;
        break;
      }
    }
  }

  // No more work.
  if(np == 0)
    return;

  amrex::Gpu::DeviceVector<Real> dp(np, 0);
  amrex::Gpu::DeviceVector<Real> ro_s(np, 0);

  Real* p_dp = dp.data();
  Real* p_ro_s = ro_s.data();

  SolidsPhase::SOLIDS_t solid;
  solid = IC::ic[icv].solids[type];

  // Setup particle diameters
  std::string ic_dp_dist_str = solid.diameter.distribution;

  if(ic_dp_dist_str.compare("normal") == 0)
  {
    nor_rno(dp, solid.diameter.mean, solid.diameter.std,
            solid.diameter.min, solid.diameter.max);
  }
  else if(ic_dp_dist_str.compare("uniform") == 0)
  {
    uni_rno(dp, solid.diameter.min,  solid.diameter.max);
  }
  else
  {
    const Real ic_dp_mean = solid.diameter.mean;

    amrex::ParallelFor(np, [p_dp,ic_dp_mean]
      AMREX_GPU_DEVICE (int p) noexcept { p_dp[p] = ic_dp_mean; });
  }

  std::string ic_ro_s_dist_str = solid.density.distribution;

  if(ic_ro_s_dist_str.compare("normal") == 0)
  {
    nor_rno(ro_s, solid.density.mean, solid.density.std,
            solid.density.min, solid.density.max);
  }
  else if(ic_ro_s_dist_str.compare("uniform") == 0)
  {
    uni_rno(ro_s, solid.density.min, solid.density.max);
  }
  else
  {
    const Real ic_ro_s_mean = solid.density.mean;

    amrex::ParallelFor(np, [p_ro_s,ic_ro_s_mean]
      AMREX_GPU_DEVICE (int p) noexcept { p_ro_s[p] = ic_ro_s_mean; });
  }

  pc = init_pc;

  const Real ic_u_s = solid.velocity[0];
  const Real ic_v_s = solid.velocity[1];
  const Real ic_w_s = solid.velocity[2];

  const Real statwt = solid.statwt;

  Real* p_rdata = m_rdata.data();
  int* p_idata = m_idata.data();

  const int local_nr = this->nr;
  const int local_ni = this->ni;

  amrex::ParallelFor(np,
    [p_rdata,p_idata,p_dp,p_ro_s,pc,ic_u_s,ic_v_s,ic_w_s,statwt,type,local_nr,local_ni,
     local_cg_dem=DEM::cg_dem]
    AMREX_GPU_DEVICE (int p) noexcept
  {
    // Real pvol = (M_PI/6.0) * (dp[p]*dp[p]*dp[p]); // UNUSED_VARIABLE
    // If coarse-grain DEM is activated
    if (local_cg_dem)
    {
       p_dp[p] = std::cbrt(statwt) * p_dp[p];
    }
    const int local_pc = pc + p;

    //< Radius................. 4
    p_rdata[local_pc*local_nr + 3] = 0.5 * p_dp[p];

    //< Density................ 5
    p_rdata[local_pc*local_nr + 4] = p_ro_s[p];

    //< Linear velocity........ 6,7,8
    p_rdata[local_pc*local_nr + 5] = ic_u_s;
    p_rdata[local_pc*local_nr + 6] = ic_v_s;
    p_rdata[local_pc*local_nr + 7] = ic_w_s;

    //< statistical weight ..... 9
    p_rdata[local_pc*local_nr + 8] = statwt;

    //< Type................... 1
    p_idata[local_pc*local_ni + 0] = type+1;
  });

  pc += np;

  Gpu::synchronize();

  return;
}

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: hex_close_pack                                          !
//                                                                      !
//  Purpose: Generate initial solids packing based on hexagonal close   !
//           packing of mono-sized spheres.                             !
//                                                                      !
//  TODO: * generalize fill direction to follow gravity.                !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void
ParticlesGenerator::hex_close_pack (const int icv,
                                    const int type,
                                    const IntVect& lo,
                                    const IntVect& hi,
                                    int& np,
                                    int& pc,
                                    const Real dx,
                                    const Real dy,
                                    const Real dz,
                                    const amrex::GpuArray<Real, 3>& plo)
{
  // indices
  int i_w, i_e, j_s, j_n, k_b, k_t;

  RealVect ic_dlo, ic_dhi;
  Real max_dp, max_rp;

  const Real sqrt3 = std::sqrt(3.0);
  const Real sqrt6o3x2 = 2.0*std::sqrt(6.0)/3.0;

  IntVect max_seed, seed_lo, seed_hi, delta_bx;

  calc_cell_ic(dx, dy, dz,
               IC::ic[icv].region->lo(),
               IC::ic[icv].region->hi(),
               plo.data(),
               i_w, i_e, j_s, j_n, k_b, k_t);

  const Real y_s(IC::ic[icv].region->lo(1));
  const Real z_b(IC::ic[icv].region->lo(2));
  const Real y_n(IC::ic[icv].region->hi(1));
  const Real z_t(IC::ic[icv].region->hi(2));

  // Start/end of IC domain bounds
  ic_dlo[0] = (amrex::max(lo[0], i_w)) * dx;
  ic_dlo[1] = (amrex::max(lo[1], j_s)) * dy;
  ic_dlo[2] = (amrex::max(lo[2], k_b)) * dz;

  ic_dhi[0] = (amrex::min(hi[0], i_e)+1) * dx;
  ic_dhi[1] = (amrex::min(hi[1], j_n)+1) * dy;
  ic_dhi[2] = (amrex::min(hi[2], k_t)+1) * dz;

  // physical volume of IC region
  const Real ic_vol = IC::ic[icv].region->volume();

  const Real mean = IC::ic[icv].solids[type].diameter.mean;

  // Spacing is based on maximum particle size
  if(IC::ic[icv].solids[type].diameter.max > 0.0)
    max_dp = IC::ic[icv].solids[type].diameter.max;
  else
    max_dp = mean;

  max_rp = 0.5 * max_dp;

  // Particle count is based on mean particle size
  const long seed =
    static_cast<long>(ic_vol * IC::ic[icv].solids[type].volfrac / ((M_PI/6.0)*mean*mean*mean));

  // Total to seed over the whole IC region
  max_seed[1] = static_cast<int>((y_n - y_s - max_dp) / max_dp);
  max_seed[2] = static_cast<int>((z_t - z_b - max_dp) / (sqrt3*max_rp));
  max_seed[0] = static_cast<int>(seed / (max_seed[1]*max_seed[2]));

  // local grid seed loop hi/lo
  seed_lo[0] = static_cast<int>(std::round((ic_dlo[0] - i_w*dx) / ((sqrt6o3x2) * max_rp)));
  seed_lo[1] = static_cast<int>(std::round((ic_dlo[1] - j_s*dy) / max_dp));
  seed_lo[2] = static_cast<int>(std::round((ic_dlo[2] - k_b*dz) / (sqrt3 * max_rp)));

  seed_hi[0] = static_cast<int>(std::round((ic_dhi[0] - i_w*dx) / ((sqrt6o3x2) * max_rp) - seed_lo[1]*max_dp));
  seed_hi[1] = static_cast<int>(std::round((ic_dhi[1] - j_s*dy) /  max_dp - seed_lo[1]*max_dp));
  seed_hi[2] = static_cast<int>(std::round((ic_dhi[2] - k_b*dz) / (sqrt3 * max_rp) - seed_lo[1]*max_dp));

  seed_hi[0] = amrex::min(max_seed[0], seed_hi[0]-1);
  seed_hi[1] = amrex::min(max_seed[1], seed_hi[1]-1);
  seed_hi[2] = amrex::min(max_seed[2], seed_hi[2]-1);

  const Box bx(seed_lo, seed_hi);

  delta_bx[0] = amrex::max(0, seed_hi[0] - seed_lo[0] + 1);
  delta_bx[1] = amrex::max(0, seed_hi[1] - seed_lo[1] + 1);
  delta_bx[2] = amrex::max(0, seed_hi[2] - seed_lo[2] + 1);

  np = delta_bx[0] * delta_bx[1] * delta_bx[2];
  grow_pdata(pc + np);

  Real* p_rdata = m_rdata.data();

  const int local_nr = this->nr;

  amrex::ParallelFor(bx,
    [p_rdata,seed_lo,delta_bx,max_rp,i_w,j_s,k_b,local_nr,pc,sqrt6o3x2,sqrt3,dx,dy,dz,plo]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const Real* plo_ptr = plo.data();

      const int local_i = i - seed_lo[0];
      const int local_j = j - seed_lo[1];
      const int local_k = k - seed_lo[2];

      const int local_pc =
        pc + (local_j + local_k*delta_bx[1] + local_i*delta_bx[1]*delta_bx[2]);

      p_rdata[local_pc*local_nr + 0] = plo_ptr[0] + i_w*dx + max_rp*(1. + i*sqrt6o3x2);
      p_rdata[local_pc*local_nr + 1] = plo_ptr[1] + j_s*dy + max_rp*(1. + 2.*j + ((i+k)%2));
      p_rdata[local_pc*local_nr + 2] = plo_ptr[2] + k_b*dz + max_rp*(1. + sqrt3*(k+((i%2)/3.)));
    });

  pc += np;

  return;
}

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: n_cube_per_fill                                         !
//                                                                      !
//  Purpose: Generate initial solids packing based on putting one       !
//           per cell                                                   !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void
ParticlesGenerator::n_cube_per_fill (const int icv,
                                     const int type,
                                     const int base,
                                     const IntVect& lo,
                                     const IntVect& hi,
                                     int& np,
                                     int& pc,
                                     const Real dx,
                                     const Real dy,
                                     const Real dz,
                                     const amrex::GpuArray<Real, 3>& plo)
{
  // This assertion should not be needed.
  AMREX_ALWAYS_ASSERT(dx==dy);
  AMREX_ALWAYS_ASSERT(dx==dz);

  int i_w, i_e, j_s, j_n, k_b, k_t;
  calc_cell_ic(dx, dy, dz,
               IC::ic[icv].region->lo(),
               IC::ic[icv].region->hi(),
               plo.data(),
               i_w, i_e, j_s, j_n, k_b, k_t);

  const IntVect ic_lo(AMREX_D_DECL(i_w, j_s, k_b));
  const IntVect ic_hi(AMREX_D_DECL(i_e, j_n, k_t));

  const Box mfi_bx(lo, hi);
  const Box ic_bx(ic_lo, ic_hi);

  if(!ic_bx.intersects(mfi_bx)){
    amrex::Print() << "Nothing to do. No intersection!\n";
    return;
  }

  // Intersection of mfi box and ic box
  const Box bx_int = mfi_bx&ic_bx;

  const Real mean = IC::ic[icv].solids[type].diameter.mean;
  Real max_dp(mean);

  // Spacing is based on maximum particle size
  if(IC::ic[icv].solids[type].diameter.max > 0.0)
    max_dp = IC::ic[icv].solids[type].diameter.max;

  // Scaled dx for particle spacing
  const Real dxn = dx/static_cast<Real>(base);

  if (max_dp > dxn) {
    std::vector<std::string> regions;
    amrex::ParmParse pp("ic");
    pp.queryarr("regions", regions);
    amrex::Print() << "\n\n";
    amrex::Print() << "**************************************************************\n";
    amrex::Print() << "  Invalid particle initialization. The spacing for " << base << "-cube\n";
    amrex::Print() << "  packing is less than the maximum particle diameter. \n";
    amrex::Print() << "  Initial Condition Name: " << regions[icv] << "\n";
    amrex::Print() << "  Fix the inputs file.\n";
    amrex::Print() << "**************************************************************\n";
    amrex::Print() << "\n\n";
    amrex::Abort("Fix the inputs file.");
  }

  // indices
  amrex::IntVect seed_lo, seed_hi, delta_bx;

  seed_lo[0] = bx_int.smallEnd(0)*base;
  seed_lo[1] = bx_int.smallEnd(1)*base;
  seed_lo[2] = bx_int.smallEnd(2)*base;

  seed_hi[0] = (bx_int.bigEnd(0) + 1)*base - 1;
  seed_hi[1] = (bx_int.bigEnd(1) + 1)*base - 1;
  seed_hi[2] = (bx_int.bigEnd(2) + 1)*base - 1;

  const Box bx(seed_lo, seed_hi);

  delta_bx[0] = amrex::max(0, seed_hi[0] - seed_lo[0] + 1);
  delta_bx[1] = amrex::max(0, seed_hi[1] - seed_lo[1] + 1);
  delta_bx[2] = amrex::max(0, seed_hi[2] - seed_lo[2] + 1);

  np = delta_bx[0] * delta_bx[1] * delta_bx[2];

  grow_pdata(pc + np);

  Real* p_rdata = m_rdata.data();

  const int local_nr = this->nr;

  amrex::ParallelFor(bx, [p_rdata,seed_lo,delta_bx,dxn,plo,pc,local_nr]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const Real* plo_ptr = plo.data();

      const int local_i = i - seed_lo[0];
      const int local_j = j - seed_lo[1];
      const int local_k = k - seed_lo[2];

      const int local_pc =
        pc + (local_i + local_k*delta_bx[0] + local_j*delta_bx[0]*delta_bx[2]);

      p_rdata[local_pc*local_nr + 0] = plo_ptr[0] + (i + 0.5) * dxn;
      p_rdata[local_pc*local_nr + 1] = plo_ptr[1] + (j + 0.5) * dxn;
      p_rdata[local_pc*local_nr + 2] = plo_ptr[2] + (k + 0.5) * dxn;
    });

  pc += np;

  return;
}




//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: random fill                                             !
//                                                                      !
//  Purpose: Generate initial solids packing based on randomly placing  !
//           particles in the ic region.                                !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void
ParticlesGenerator::random_fill_dem (const int icv,
                                     const int type,
                                     const IntVect& lo,
                                     const IntVect& hi,
                                     int& np,
                                     int& pc,
                                     const Real dx,
                                     const Real dy,
                                     const Real dz,
                                     const amrex::GpuArray<Real, 3>& plo,
                                     const bool /*fix_seed*/)
{
    // indices
  int i_w, i_e, j_s, j_n, k_b, k_t;

  const int maxfails = 1000;

  RealVect ic_dlo, ic_dhi, ic_len;
  Real max_dp, max_rp;

  calc_cell_ic(dx, dy, dz,
               IC::ic[icv].region->lo(),
               IC::ic[icv].region->hi(),
               plo.data(),
               i_w, i_e, j_s, j_n, k_b, k_t);

  // Start/end of IC domain bounds
  ic_dlo[0] = (amrex::max(lo[0], i_w)) * dx;
  ic_dlo[1] = (amrex::max(lo[1], j_s)) * dy;
  ic_dlo[2] = (amrex::max(lo[2], k_b)) * dz;

  ic_dhi[0] = (amrex::min(hi[0], i_e)+1) * dx;
  ic_dhi[1] = (amrex::min(hi[1], j_n)+1) * dy;
  ic_dhi[2] = (amrex::min(hi[2], k_t)+1) * dz;

  // physical volume of IC region intersecting this grid
  const Real ic_vol = (ic_dhi[0] - ic_dlo[0]) *
                             (ic_dhi[1] - ic_dlo[1]) *
                             (ic_dhi[2] - ic_dlo[2]);

  Real mean = IC::ic[icv].solids[type].diameter.mean;

  // Spacing is based on maximum particle size
  if(IC::ic[icv].solids[type].diameter.max > 0.)
    max_dp = IC::ic[icv].solids[type].diameter.max;
  else
    max_dp = IC::ic[icv].solids[type].diameter.mean;

  // If coarse-grain DEM is activated
  if (DEM::cg_dem)
  {
     Real statwt = IC::ic[icv].solids[type].statwt;
     mean = std::cbrt(statwt) * mean;
     max_dp = std::cbrt(statwt) * max_dp;
  }

  // Particle count is based on mean particle size
  const int seed =
    static_cast<int>(ic_vol * IC::ic[icv].solids[type].volfrac / ((M_PI/6.0)*mean*mean*mean));

  max_rp = 0.5 * max_dp;

  ic_len[0] = ic_dhi[0] - ic_dlo[0] - max_dp;
  ic_len[1] = ic_dhi[1] - ic_dlo[1] - max_dp;
  ic_len[2] = ic_dhi[2] - ic_dlo[2] - max_dp;

  ic_dlo[0] += max_rp;
  ic_dlo[1] += max_rp;
  ic_dlo[2] += max_rp;

  const Real mindist = (1.01*max_dp)*(1.01*max_dp);

  const RealVect Oodx = {1/dx, 1/dy, 1/dz};

  int nb(8);

  IntVect delta_bx;
  delta_bx[0] = hi[0] - lo[0] + 1;
  delta_bx[1] = hi[1] - lo[1] + 1;
  delta_bx[2] = hi[2] - lo[2] + 1;

  Vector<int> pinc(delta_bx[0]*delta_bx[1]*delta_bx[2], 0);
  Vector<int> pbin(delta_bx[0]*delta_bx[1]*delta_bx[2]*nb, 0);

  IntVect max_fitting_parts;
  max_fitting_parts[0] = static_cast<int>(amrex::Math::ceil(ic_len[0]/max_dp));
  max_fitting_parts[1] = static_cast<int>(amrex::Math::ceil(ic_len[1]/max_dp));
  max_fitting_parts[2] = static_cast<int>(amrex::Math::ceil(ic_len[2]/max_dp));

  const int max_np = max_fitting_parts[0]*max_fitting_parts[1]*max_fitting_parts[2];
  grow_pdata(max_np);

  amrex::ResetRandomSeed(ParallelDescriptor::MyProc()+1);

  amrex::Gpu::HostVector<Real> h_rdata(m_rdata.size());
  Gpu::copyAsync(Gpu::deviceToHost, m_rdata.begin(), m_rdata.end(), h_rdata.begin());
  Gpu::synchronize();
  Real* p_rdata = h_rdata.data();

  np = 0;
  int iterations(0);

  while(iterations < seed)
  {
    bool stop = false;
    int fails(0);

    while(fails < maxfails && (!stop)) {
      RealVect pos;
      pos[0] = ic_dlo[0] + ic_len[0]*amrex::Random();
      pos[1] = ic_dlo[1] + ic_len[1]*amrex::Random();
      pos[2] = ic_dlo[2] + ic_len[2]*amrex::Random();

      // Grid containing the new particle
      IntVect idx;
      idx[0] = static_cast<int>(amrex::Math::floor(pos[0]*Oodx[0]));
      idx[1] = static_cast<int>(amrex::Math::floor(pos[1]*Oodx[1]));
      idx[2] = static_cast<int>(amrex::Math::floor(pos[2]*Oodx[2]));

      // Local grid search for collisions.
      int overlaps = 0;

      IntVect seed_lo, seed_hi;
      seed_lo[0] = amrex::max(lo[0], idx[0]-1);
      seed_lo[1] = amrex::max(lo[1], idx[1]-1);
      seed_lo[2] = amrex::max(lo[2], idx[2]-1);

      seed_hi[0] = amrex::min(hi[0], idx[0]+1);
      seed_hi[1] = amrex::min(hi[1], idx[1]+1);
      seed_hi[2] = amrex::min(hi[2], idx[2]+1);

      for(int k = seed_lo[2]; k <= seed_hi[2]; ++k) {
        for(int j = seed_lo[1]; j <= seed_hi[1]; ++j) {
          for(int i = seed_lo[0]; i <= seed_hi[0]; ++i) {
            const int local_i = i - lo[0];
            const int local_j = j - lo[1];
            const int local_k = k - lo[2];

            const int upper_limit =
              pinc[local_i + local_j*delta_bx[0] + local_k*delta_bx[0]*delta_bx[1]];

            for(int p = 0; p < upper_limit; p++)
            {
              const Real* plo_ptr = plo.data();

              const int local_pc =
                pbin[local_i + local_j*delta_bx[0] +
                     local_k*delta_bx[0]*delta_bx[1] +
                     p*delta_bx[0]*delta_bx[1]*delta_bx[2]];

              const Real dist_x = p_rdata[local_pc*nr + 0] - pos[0] - plo_ptr[0];
              const Real dist_y = p_rdata[local_pc*nr + 1] - pos[1] - plo_ptr[1];
              const Real dist_z = p_rdata[local_pc*nr + 2] - pos[2] - plo_ptr[2];

              const Real dist = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;

              if(dist < mindist)
                overlaps++;
            }
          }
        }
      }

      if(overlaps == 0)
      {
        p_rdata[np*nr + 0] = pos[0] + plo[0];
        p_rdata[np*nr + 1] = pos[1] + plo[1];
        p_rdata[np*nr + 2] = pos[2] + plo[2];

        const int local_idx_x = idx[0] - lo[0];
        const int local_idx_y = idx[1] - lo[1];
        const int local_idx_z = idx[2] - lo[2];

        const int local_idx =
          local_idx_x + local_idx_y*delta_bx[0] +
          local_idx_z*delta_bx[0]*delta_bx[1];

        const int pinc_value = pinc[local_idx];
        pbin[local_idx + pinc_value*delta_bx[0]*delta_bx[1]*delta_bx[2]] = np;

        pinc[local_idx] = pinc_value + 1;

        if(pinc_value + 2 >= nb)
        {
          nb += 2;
          pbin.resize(delta_bx[0]*delta_bx[1]*delta_bx[2]*nb, 0);
        }

        np++;

        stop = true;
      }
      else
      {
        fails++;
      }
    }

    iterations++;
  }

  pc += np;

  Gpu::copyAsync(Gpu::hostToDevice, h_rdata.begin(), h_rdata.end(), m_rdata.begin());
  Gpu::synchronize();

  return;
}

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: random fill                                             !
//                                                                      !
//  Purpose: Generate initial solids packing based on randomly placing  !
//           particles in the ic region.                                !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void
ParticlesGenerator::random_fill_pic (const int icv,
                                     const int type,
                                     const IntVect& lo,
                                     const IntVect& hi,
                                     int& np,
                                     int& pc,
                                     const Real dx,
                                     const Real dy,
                                     const Real dz,
                                     const amrex::GpuArray<Real, 3>& plo,
                                     const bool /*fix_seed*/)
{

  const Real* ic_rlo = IC::ic[icv].region->lo();
  const Real* ic_rhi = IC::ic[icv].region->hi();

  amrex::IntVect seed_lo, seed_hi;

  seed_lo[0] = amrex::max<int>(lo[0], static_cast<int>(amrex::Math::ceil((ic_rlo[0]-plo[0])/dx - .5)));
  seed_lo[1] = amrex::max<int>(lo[1], static_cast<int>(amrex::Math::ceil((ic_rlo[1]-plo[1])/dy - .5)));
  seed_lo[2] = amrex::max<int>(lo[2], static_cast<int>(amrex::Math::ceil((ic_rlo[2]-plo[2])/dz - .5)));

  seed_hi[0] = amrex::min<int>(hi[0], static_cast<int>(amrex::Math::floor((ic_rhi[0]-plo[0])/dx - .5)));
  seed_hi[1] = amrex::min<int>(hi[1], static_cast<int>(amrex::Math::floor((ic_rhi[1]-plo[1])/dy - .5)));
  seed_hi[2] = amrex::min<int>(hi[2], static_cast<int>(amrex::Math::floor((ic_rhi[2]-plo[2])/dz - .5)));

  const Box bx(seed_lo, seed_hi);

  amrex::IntVect delta_bx;

  delta_bx[0] = amrex::max(0, seed_hi[0] - seed_lo[0] + 1);
  delta_bx[1] = amrex::max(0, seed_hi[1] - seed_lo[1] + 1);
  delta_bx[2] = amrex::max(0, seed_hi[2] - seed_lo[2] + 1);

  // Number of cells in IC region
  const Real ic_cells = delta_bx[0] * delta_bx[1] * delta_bx[2];

  // Particle properties: Mean diameter, stat weight, volume
  const Real mean = IC::ic[icv].solids[type].diameter.mean;
  const Real statwt(IC::ic[icv].solids[type].statwt);
  const Real parcel_volume = statwt*((M_PI/6.0)*mean*mean*mean);

  // Parcels per cell to satisfy the IC condition
  const Real parcels_per_cell =
      (dx * dy * dz * IC::ic[icv].solids[type].volfrac / parcel_volume);

  const int whole_parcels_per_cell = static_cast<int>(amrex::Math::floor(parcels_per_cell));

  // This is the total number of particles that will be generated.
  np = static_cast<int>(ic_cells) * whole_parcels_per_cell;
  grow_pdata(pc + np);

  Real* p_rdata = m_rdata.data();
  const int local_nr = this->nr;

  amrex::ResetRandomSeed(ParallelDescriptor::MyProc()+1);

  amrex::ParallelForRNG(bx,
    [p_rdata,seed_lo,delta_bx,local_nr,dx,dy,dz,plo,pc,whole_parcels_per_cell]
    AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
    {
      const Real* plo_ptr = plo.data();

      const int local_i = i - seed_lo[0];
      const int local_j = j - seed_lo[1];
      const int local_k = k - seed_lo[2];

      const int local_pc = pc + whole_parcels_per_cell *
        (local_i + local_k*delta_bx[0] + local_j*delta_bx[0]*delta_bx[2]);

      const Real xlo = plo_ptr[0] + (i * dx);
      const Real ylo = plo_ptr[1] + (j * dy);
      const Real zlo = plo_ptr[2] + (k * dz);

      for(int pseed(0); pseed < whole_parcels_per_cell; pseed++){

        p_rdata[(local_pc + pseed)*local_nr + 0] = xlo + dx*amrex::Random(engine);
        p_rdata[(local_pc + pseed)*local_nr + 1] = ylo + dy*amrex::Random(engine);
        p_rdata[(local_pc + pseed)*local_nr + 2] = zlo + dz*amrex::Random(engine);
      }
    }
  );
  pc += np;

  // TODO: At this point, we have seeded the bulk of the parcels. However there
  // are two oversites --
  // 1) This will get us close to the requested solids volume fraction, but
  //    not exact. We need to deal with the remainder.
  // 2) If the IC region is dilute (with less than one parcel per cell), this
  //    will not create any parcles. This could be addressed by adding a
  //    stride to the above logic so we skip past the stride cells.

  return;
}




//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  SUBROUTINE: GENERATE_PARTICLE_CONFIG                                !
//                                                                      !
//  Purpose: Generate particle configuration based on maximum particle  !
//           radius and filling from top to bottom within specified     !
//           bounds                                                     !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void
ParticlesGenerator::generate_prop (const int nrp, ParticleTileType& particles)
{
  auto& aos = particles.GetArrayOfStructs();
  ParticleType* pstruct = aos().dataPtr();

  auto& soa = particles.GetStructOfArrays();
  auto p_realarray = soa.realarray();
  auto p_intarray = soa.intarray();

  Real* p_rdata = m_rdata.data();
  int* p_idata = m_idata.data();

  const int local_nr = this->nr;
  const int local_ni = this->ni;

  const Real picmulti = (DEM::solve) ? 1.0 : 0.0;

  amrex::ParallelFor(nrp, [pstruct,p_realarray,p_intarray,p_rdata,p_idata,
      local_nr,local_ni,picmulti]
    AMREX_GPU_DEVICE (int p) noexcept
  {
    ParticleType& part = pstruct[p];

    part.pos(0) = p_rdata[p*local_nr + 0];
    part.pos(1) = p_rdata[p*local_nr + 1];
    part.pos(2) = p_rdata[p*local_nr + 2];

    Real rad = p_rdata[p*local_nr + 3];
    Real rho = p_rdata[p*local_nr + 4];

    Real vol  = (4.0/3.0)*M_PI*rad*rad*rad;
    Real mass = vol * rho;
    Real omoi = 2.5/(mass * rad*rad);

    p_intarray[SoAintData::phase][p] = p_idata[p*local_ni + 0];
    p_intarray[SoAintData::state][p] = 1;

    p_realarray[SoArealData::velx][p] = p_rdata[p*local_nr + 5];
    p_realarray[SoArealData::vely][p] = p_rdata[p*local_nr + 6];
    p_realarray[SoArealData::velz][p] = p_rdata[p*local_nr + 7];

    p_realarray[SoArealData::statwt][p] = p_rdata[p*local_nr + 8];

    p_realarray[SoArealData::radius][p] = rad;
    p_realarray[SoArealData::density][p] = rho;

    p_realarray[SoArealData::volume][p] = vol;
    p_realarray[SoArealData::mass][p] = mass;
    p_realarray[SoArealData::oneOverI][p] = omoi*picmulti;

    p_realarray[SoArealData::omegax][p] = 0.0;
    p_realarray[SoArealData::omegay][p] = 0.0;
    p_realarray[SoArealData::omegaz][p] = 0.0;

    p_realarray[SoArealData::dragcoeff][p] = 0.0;

    p_realarray[SoArealData::dragx][p] = 0.0;
    p_realarray[SoArealData::dragy][p] = 0.0;
    p_realarray[SoArealData::dragz][p] = 0.0;

    p_intarray[SoAintData::phase][p] = p_idata[p*local_ni + 0];
    p_intarray[SoAintData::state][p] = 1;
  });

  Gpu::synchronize();

  m_rdata.clear();
  m_idata.clear();

  return;
}

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
//                                                                     !
//                                                                     !
//                                                                     !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
void
ParticlesGenerator::nor_rno (amrex::Gpu::DeviceVector<Real>& dp,
                             const Real mean,
                             const Real sigma,
                             const Real dp_min,
                             const Real dp_max)
{
  // Local variables
  //-----------------------------------------------
  bool debug = false;
  //-----------------------------------------------

  const unsigned int nsize = dp.size();

  const int nsize_half = static_cast<int>(amrex::Math::ceil(nsize/2.));

  const Real tolerance =
    std::numeric_limits<Real>::epsilon();

  Real* p_dp = dp.data();

  const int maxfails = 10000;
  int fails(0);

#ifdef AMREX_USE_GPU
  Gpu::DeviceScalar<int> fails_gpu(fails);
  int* p_fails = fails_gpu.dataPtr();
#endif

  amrex::ParallelForRNG(nsize_half,
    [p_dp,dp_min,dp_max,tolerance,sigma,mean,nsize,
#ifdef AMREX_USE_GPU
     p_fails]
#else
     &fails]
#endif
    AMREX_GPU_DEVICE (int i, RandomEngine const& engine) noexcept
    {
      Real x(0);
      Real y(0);

      Real dp1 = dp_min - 1;
      Real dp2 = dp_min - 1;

      int iterations(0);

      while(!(dp1 >= dp_min && dp1 <= dp_max && dp2 >= dp_min && dp2 <= dp_max) &&
            iterations < maxfails)
      {
        Real w(1.1);
        while(w > 1 || amrex::Math::abs(w-1) < tolerance)
        {
          x = 2*amrex::Random(engine) - 1;
          y = 2*amrex::Random(engine) - 1;
          w = x*x + y*y;
        }

        w = std::sqrt((-2.0 * std::log(w)) / w);

        dp1 = x * w * sigma + mean;
        dp2 = y * w * sigma + mean;

        p_dp[2*i] = dp1;

        if(size_t(2*i+1) < nsize)
          p_dp[2*i+1] = dp2;

        iterations++;
      }

      if(!(iterations < maxfails))
      {
#ifdef AMREX_USE_GPU
        Gpu::Atomic::Add(p_fails, 1);
#else
        fails++;
#endif
      }
    });

  Gpu::synchronize();

#ifdef AMREX_USE_GPU
  fails = fails_gpu.dataValue();
#endif

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fails == 0,
      "In ParticlesGenerator: nor_rno routine failed.");

  if(debug)
  {
    Real lmean(0), lvariance(0), lsigma(0);
    std::accumulate(dp.begin(), dp.end(), lmean);

    lmean /= nsize;

    for(size_t i = 0; i < nsize; i++)
      lvariance += (dp[i]-lmean)*(dp[i]-lmean);

    lvariance /= nsize;
    lsigma = std::sqrt(lvariance);

    std::string divider("    +");
    divider += "-----------";
    divider += "+-----------------";
    divider += "+-----------------";
    divider += "+";

    std::string header("    |");
    header += "           ";
    header += "|    Specified    |    Computed     |";

    amrex::Print() << "   " << std::endl;
    amrex::Print() << divider;
    amrex::Print() << header;
    amrex::Print() << divider;

    amrex::Print() << "    |" << "  Mean     "
                   << std::scientific << std::setw(15) << std::setprecision(6)
                   << "| " << mean << " "
                   << "| " << lmean << " |" << std::endl;

    amrex::Print() << divider;

    amrex::Print() << "    |" << "  Sigma    "
                   << std::scientific << std::setw(15) << std::setprecision(6)
                   << "| " << sigma << " "
                   << "| " << lsigma << " |" << std::endl;

    amrex::Print() << divider;
  }

  return;
}

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
//                                                                     !
//                                                                     !
//                                                                     !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
void
ParticlesGenerator::uni_rno (amrex::Gpu::DeviceVector<Real>& dp,
                             const Real dp_min,
                             const Real dp_max)
{
  Real lscale = dp_max - dp_min;

  const unsigned int nsize = dp.size();

  Real* p_dp = dp.data();

  amrex::ParallelForRNG(nsize, [p_dp,dp_min,lscale]
    AMREX_GPU_DEVICE (int lc, amrex::RandomEngine const& engine) noexcept { p_dp[lc] = dp_min + lscale*amrex::Random(engine); });

  return;
}

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
//                                                                     !
//                                                                     !
//                                                                     !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
void
ParticlesGenerator::grow_pdata (const int gsize)
{
  const size_t r_gsize = gsize*nr;
  const size_t i_gsize = gsize*ni;

  // Increase real data
  if(m_rdata.size() == 0)
  {
    const int size = amrex::max<int>(r_gsize, 1024*nr);

    m_rdata.resize(size, 0);
  }
  else if(r_gsize >= m_rdata.size())
  {
    const int size = amrex::max<int>(2*m_rdata.size(), r_gsize);

    m_rdata.resize(size, 0);
  }

  // Increase integer data
  if(m_idata.size() == 0)
  {
    const int size = amrex::max<int>(i_gsize, 1024*ni);

    m_idata.resize(size, 0);
  }
  else if(i_gsize >= m_idata.size())
  {
    const int size = amrex::max<int>(2*m_idata.size(), i_gsize);

    m_idata.resize(size, 0);
  }

  return;
}

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
//                                                                     !
//                                                                     !
//                                                                     !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
void
ParticlesGenerator::write (const int nrp,
                           ParticleTileType& particles,
                           const int nstep) const
{
  const auto& aos = particles.GetArrayOfStructs();
  const ParticleType* pstruct = aos().dataPtr();

  auto& soa = particles.GetStructOfArrays();
  auto p_realarray = soa.realarray();

  std::ostringstream nstep_stream;
  nstep_stream << nstep;

  std::string nstep_string(nstep_stream.str());
  std::string filename = "test-" + nstep_string + ".vtp";

  std::ofstream output_file;
  output_file.open(filename.c_str());

  // Write the necessary header information for a PolyData file type
  output_file << "<?xml version=\"1.0\"?>" << std::endl;
  output_file << "<VTKFile type=\"PolyData\"" << " "
              << "version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
  output_file << "   " << "<PolyData>" << std::endl;

  // Write Piece tag and identify the number of parts in the system.
  output_file << "      " << "<Piece NumberOfPoints=\""
              << std::fixed << std::setw(10) << std::setprecision(10)
              << std::setfill('0') << nrp << "\" NumberOfVerts=\"0\" "
              << "NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">" << std::endl;

  output_file << "         " << "<PointData>" << std::endl;

  output_file << "            "
    << "<DataArray type=\"Float32\" Name=\"radius\" NumberOfComponents=\"1\" format=\"ascii\">"
    << std::endl;

  for(int lc1 = 0; lc1 < nrp; ++lc1)
  {
    output_file << "               "
                << std::scientific << std::setw(13) << std::setprecision(6)
                << std::setfill(' ') << Real(p_realarray[SoArealData::radius][lc1]) << std::endl;
  }

  output_file << "            " << "</DataArray>" << std::endl;

  output_file << "            "
    << "<DataArray type=\"Float32\" Name=\"density\" NumberOfComponents=\"1\" format=\"ascii\">"
    << std::endl;

  for(int lc1 = 0; lc1 < nrp; ++lc1)
  {
    output_file << "               "
                << std::scientific << std::setw(13) << std::setprecision(6)
                << Real(p_realarray[SoArealData::density][lc1]) << std::endl;
  }

  output_file << "            " << "</DataArray>" << std::endl;

  output_file << "         " << "</PointData>" << std::endl;

  output_file << "         " << "<Points>" << std::endl;

  output_file << "            " << "<DataArray type=\"Float32\" "
              << "Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">"
              << std::endl;

  for(int lc1 = 0; lc1 < nrp; lc1++)
  {
    const ParticleType& part = pstruct[lc1];

    output_file << "               "
                << std::scientific << std::setw(13) << std::setprecision(6)
                << Real(part.pos(0)) << "    "
                << Real(part.pos(1)) << "    "
                << Real(part.pos(2)) << std::endl;
  }

  output_file << "            " << "</DataArray>" << std::endl;
  output_file << "         " << "</Points>" << std::endl;

  // Write tags for data not included (vtp format style)
  output_file << "         " << "<Verts></Verts>" << std::endl;
  output_file << "         " << "<Lines></Lines>" << std::endl;
  output_file << "         " << "<Strips></Strips>" << std::endl;
  output_file << "         " << "<Polys></Polys>" << std::endl;

  output_file << "      " << "</Piece>" << std::endl;
  output_file << "   " << "</PolyData>" << std::endl;
  output_file << "</VTKFile>" << std::endl;

  output_file.close();

  return;
}
