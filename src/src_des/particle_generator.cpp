#include <particle_generator.H>
#include <AMReX_AmrParGDB.H>
#include <param_mod_F.H>
#include <ic_mod_F.H>
#include <discretelement_mod_F.H>
#include <calc_cell_F.H>

#include <limits>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>

#include <MFIX_DEM_Parms.H>

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
                              const amrex::Real dx,
                              const amrex::Real dy,
                              const amrex::Real dz)
{
  const int init_pc(pc);

  const amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();

  int np(0);
  int icv(1);
  int type(1);

  // Get the IC index
  int icv0(1);

  for(; icv0 <= get_dim_ic(); icv0++)
  {
    if(ic_defined_cpp(icv0) and std::abs(get_ic_ep_g(icv0)-1.0) > tolerance)
    {
      // Get the solids type index
      int type0(1);

      for(; type0 <= DEM::solve; type0++)
        if(get_ic_ep_s(icv0,type0) > tolerance)
          break;

      char ic_pack_type[16];
      get_ic_pack_type(icv0, ic_pack_type);

      std::string ic_pack_type_str(ic_pack_type);

      int np0(0);

      if(ic_pack_type_str.compare("HCP") == 0)
        hex_close_pack(icv0, type0, lo, hi, np0, pc, dx, dy, dz);
      else if(ic_pack_type_str.compare("RANDOM") == 0)
        random_fill(icv0, type0, lo, hi, np0, pc, dx, dy, dz, false);
      else if(ic_pack_type_str.compare("PSEUDO_RANDOM") == 0)
        random_fill(icv0, type0, lo, hi, np0, pc, dx, dy, dz, true);
      else if(ic_pack_type_str.compare("ONEPER") == 0)
        one_per_fill(icv0, type0, lo, hi, np0, pc, dx, dy, dz);
      else if(ic_pack_type_str.compare("EIGHTPER") == 0)
        eight_per_fill(icv0, type0, lo, hi, np0, pc, dx, dy, dz);
      else
      {
        amrex::Print() << "Unknown particle generator fill type" << std::endl;
        exit(1000);
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

  amrex::Gpu::ManagedVector<amrex::Real> dp(np, 0);
  amrex::Gpu::ManagedVector<amrex::Real> ro_s(np, 0);

  amrex::Real* p_dp = dp.data();
  amrex::Real* p_ro_s = ro_s.data();

  // Setup particle diameters
  char ic_dp_dist[16];
  get_ic_dp_dist(icv, type, ic_dp_dist);

  std::string ic_dp_dist_str(ic_dp_dist);

  if(ic_dp_dist_str.compare("NORMAL") == 0)
  {
    nor_rno(dp, get_ic_dp_mean(icv,type), get_ic_dp_std(icv,type),
            get_ic_dp_min(icv,type), get_ic_dp_max(icv,type));
  }
  else if(ic_dp_dist_str.compare("UNIFORM") == 0)
  {
    uni_rno(dp, get_ic_dp_min(icv,type), get_ic_dp_max(icv,type));
  }
  else
  {
    const amrex::Real ic_dp_mean = get_ic_dp_mean(icv, type);

    amrex::ParallelFor(np, [p_dp,ic_dp_mean]
      AMREX_GPU_DEVICE (int p) noexcept { p_dp[p] = ic_dp_mean; });
  }

  char ic_ro_s_dist[16];
  get_ic_ro_s_dist(icv, type, ic_ro_s_dist);

  std::string ic_ro_s_dist_str(ic_ro_s_dist);

  if(ic_ro_s_dist_str.compare("NORMAL") == 0)
  {
    nor_rno(ro_s, get_ic_ro_s_mean(icv,type), get_ic_ro_s_std(icv,type),
            get_ic_ro_s_min(icv,type), get_ic_ro_s_max(icv,type));
  }
  else if(ic_ro_s_dist_str.compare("UNIFORM") == 0)
  {
    uni_rno(ro_s, get_ic_ro_s_min(icv,type), get_ic_ro_s_max(icv,type));
  }
  else
  {
    const amrex::Real ic_ro_s_mean = get_ic_ro_s_mean(icv, type);

    amrex::ParallelFor(np, [p_ro_s,ic_ro_s_mean]
      AMREX_GPU_DEVICE (int p) noexcept { p_ro_s[p] = ic_ro_s_mean; });
  }

  pc = init_pc;

  const amrex::Real ic_u_s = get_ic_u_s(icv, type);
  const amrex::Real ic_v_s = get_ic_v_s(icv, type);
  const amrex::Real ic_w_s = get_ic_w_s(icv, type);

  amrex::Real* p_rdata = m_rdata.data();
  int* p_idata = m_idata.data();

  const int local_nr = this->nr;
  const int local_ni = this->ni;

  amrex::ParallelFor(np,
    [p_rdata,p_idata,p_dp,p_ro_s,pc,ic_u_s,ic_v_s,ic_w_s,type,local_nr,local_ni]
    AMREX_GPU_DEVICE (int p) noexcept
  {
    // amrex::Real pvol = (M_PI/6.0) * (dp[p]*dp[p]*dp[p]); // UNUSED_VARIABLE

    const int local_pc = pc + p;

    //< Radius................. 4
    p_rdata[local_pc*local_nr + 3] = 0.5 * p_dp[p];

    //< Density................ 5
    p_rdata[local_pc*local_nr + 4] = p_ro_s[p];

    //< Linear velocity........ 6,7,8
    p_rdata[local_pc*local_nr + 5] = ic_u_s;
    p_rdata[local_pc*local_nr + 6] = ic_v_s;
    p_rdata[local_pc*local_nr + 7] = ic_w_s;

    //< Type................... 1
    p_idata[local_pc*local_ni + 0] = type;
  });

  pc += np;

  amrex::Gpu::synchronize();

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
                                    const amrex::Real dx,
                                    const amrex::Real dy,
                                    const amrex::Real dz)
{
  // indices
  int i_w, i_e, j_s, j_n, k_b, k_t;

  amrex::RealVect ic_dlo, ic_dhi;
  amrex::Real max_dp, max_rp;

  const Real sqrt3 = std::sqrt(3.0);
  const Real sqrt6o3x2 = 2.0*std::sqrt(6.0)/3.0;

  IntVect max_seed, seed_lo, seed_hi, delta_bx;

  calc_cell_ic(dx, dy, dz, get_ic_x_w(icv), get_ic_y_s(icv), get_ic_z_b(icv),
               get_ic_x_e(icv), get_ic_y_n(icv), get_ic_z_t(icv),
               i_w, i_e, j_s, j_n, k_b, k_t);

  // Start/end of IC domain bounds
  ic_dlo[0] = (std::max(lo[0], i_w)) * dx;
  ic_dlo[1] = (std::max(lo[1], j_s)) * dy;
  ic_dlo[2] = (std::max(lo[2], k_b)) * dz;

  ic_dhi[0] = (std::min(hi[0], i_e)+1) * dx;
  ic_dhi[1] = (std::min(hi[1], j_n)+1) * dy;
  ic_dhi[2] = (std::min(hi[2], k_t)+1) * dz;

  // physical volume of IC region
  const amrex::Real ic_vol = (get_ic_x_e(icv) - get_ic_x_w(icv)) *
                             (get_ic_y_n(icv) - get_ic_y_s(icv)) *
                             (get_ic_z_t(icv) - get_ic_z_b(icv));

  // Spacing is based on maximum particle size
  if(is_defined_db_cpp(get_ic_dp_max(icv,type)))
    max_dp = get_ic_dp_max(icv,type);
  else
    max_dp = get_ic_dp_mean(icv,type);

  max_rp = 0.5 * max_dp;

  // Particle count is based on mean particle size
  const amrex::Real mean = get_ic_dp_mean(icv,type);
  const int seed =
    static_cast<int>(ic_vol * get_ic_ep_s(icv,type) / ((M_PI/6.0)*mean*mean*mean));

  // Total to seed over the whole IC region
  max_seed[1] = static_cast<int>((get_ic_y_n(icv) - get_ic_y_s(icv) - max_dp) / max_dp);
  max_seed[2] = static_cast<int>((get_ic_z_t(icv) - get_ic_z_b(icv) - max_dp) / (sqrt3*max_rp));
  max_seed[0] = static_cast<int>(seed / (max_seed[1]*max_seed[2]));

  // local grid seed loop hi/lo
  seed_lo[0] = std::round((ic_dlo[0] - i_w*dx) / ((sqrt6o3x2) * max_rp));
  seed_lo[1] = std::round((ic_dlo[1] - j_s*dy) / max_dp);
  seed_lo[2] = std::round((ic_dlo[2] - k_b*dz) / (sqrt3 * max_rp));

  seed_hi[0] = std::round((ic_dhi[0] - i_w*dx) / ((sqrt6o3x2) * max_rp) - seed_lo[1]*max_dp);
  seed_hi[1] = std::round((ic_dhi[1] - j_s*dy) /  max_dp - seed_lo[1]*max_dp);
  seed_hi[2] = std::round((ic_dhi[2] - k_b*dz) / (sqrt3 * max_rp) - seed_lo[1]*max_dp);

  seed_hi[0] = std::min(max_seed[0], seed_hi[0]-1);
  seed_hi[1] = std::min(max_seed[1], seed_hi[1]-1);
  seed_hi[2] = std::min(max_seed[2], seed_hi[2]-1);

  const Box bx(seed_lo, seed_hi);

  delta_bx[0] = std::max(0, seed_hi[0] - seed_lo[0] + 1);
  delta_bx[1] = std::max(0, seed_hi[1] - seed_lo[1] + 1);
  delta_bx[2] = std::max(0, seed_hi[2] - seed_lo[2] + 1);

  np = delta_bx[0] * delta_bx[1] * delta_bx[2];
  grow_pdata(pc + np);

  amrex::Real* p_rdata = m_rdata.data();

  const int local_nr = this->nr;

  amrex::ParallelFor(bx,
    [p_rdata,seed_lo,delta_bx,max_rp,i_w,j_s,k_b,local_nr,pc,sqrt6o3x2,sqrt3,dx,dy,dz]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int local_i = i - seed_lo[0];
      const int local_j = j - seed_lo[1];
      const int local_k = k - seed_lo[2];

      const int local_pc =
        pc + (local_j + local_k*delta_bx[1] + local_i*delta_bx[1]*delta_bx[2]);

      p_rdata[local_pc*local_nr + 0] = i_w*dx + max_rp*(1. + i*sqrt6o3x2);
      p_rdata[local_pc*local_nr + 1] = j_s*dy + max_rp*(1. + 2.*j + ((i+k)%2));
      p_rdata[local_pc*local_nr + 2] = k_b*dz + max_rp*(1. + sqrt3*(k+((i%2)/3.)));
    });

  pc += np;

  return;
}

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: one_per_fill                                            !
//                                                                      !
//  Purpose: Generate initial solids packing based on putting one       !
//           per cell                                                   !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void
ParticlesGenerator::one_per_fill (const int icv,
                                  const int type,
                                  const IntVect& lo,
                                  const IntVect& hi,
                                  int& np,
                                  int& pc,
                                  const amrex::Real dx,
                                  const amrex::Real dy,
                                  const amrex::Real dz)
{
  // indices
  int i_w, i_e, j_s, j_n, k_b, k_t;

  amrex::IntVect seed_lo, seed_hi, delta_bx;

  calc_cell_ic(dx, dy, dz, get_ic_x_w(icv), get_ic_y_s(icv), get_ic_z_b(icv),
               get_ic_x_e(icv), get_ic_y_n(icv), get_ic_z_t(icv),
               i_w, i_e, j_s, j_n, k_b, k_t);

  // local grid seed loop hi/lo
  seed_lo[0] = std::max<int>(lo[0], std::ceil(get_ic_x_w(icv)/dx - .5));
  seed_lo[1] = std::max<int>(lo[1], std::ceil(get_ic_y_s(icv)/dy - .5));
  seed_lo[2] = std::max<int>(lo[2], std::ceil(get_ic_z_b(icv)/dz - .5));

  seed_hi[0] = std::min<int>(hi[0], std::floor(get_ic_x_e(icv)/dx - .5));
  seed_hi[1] = std::min<int>(hi[1], std::floor(get_ic_y_n(icv)/dy - .5));
  seed_hi[2] = std::min<int>(hi[2], std::floor(get_ic_z_t(icv)/dz - .5));

  const Box bx(seed_lo, seed_hi);

  delta_bx[0] = std::max(0, seed_hi[0] - seed_lo[0] + 1);
  delta_bx[1] = std::max(0, seed_hi[1] - seed_lo[1] + 1);
  delta_bx[2] = std::max(0, seed_hi[2] - seed_lo[2] + 1);

  np = delta_bx[0] * delta_bx[1] * delta_bx[2];
  grow_pdata(pc + np);

  amrex::Real* p_rdata = m_rdata.data();

  const int local_nr = this->nr;

  amrex::ParallelFor(bx, [p_rdata,seed_lo,delta_bx,local_nr,dx,dy,dz,pc]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int local_i = i - seed_lo[0];
      const int local_j = j - seed_lo[1];
      const int local_k = k - seed_lo[2];

      const int local_pc =
        pc + (local_i + local_k*delta_bx[0] + local_j*delta_bx[0]*delta_bx[2]);

      p_rdata[local_pc*local_nr + 0] = (i + 0.5) * dx;
      p_rdata[local_pc*local_nr + 1] = (j + 0.5) * dy;
      p_rdata[local_pc*local_nr + 2] = (k + 0.5) * dz;
    });

  pc += np;

  return;
}

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: eight_per_fill                                          !
//                                                                      !
//  Purpose: Generate initial solids packing based on putting eight     !
//           per cell                                                   !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void
ParticlesGenerator::eight_per_fill (const int icv,
                                    const int type,
                                    const IntVect& lo,
                                    const IntVect& hi,
                                    int& np,
                                    int& pc,
                                    const amrex::Real dx,
                                    const amrex::Real dy,
                                    const amrex::Real dz)
{
  // indices
  amrex::IntVect seed_lo, seed_hi, delta_bx;

  seed_lo[0] = std::max<int>(2*lo[0], std::ceil(get_ic_x_w(icv)*(2/dx) - .5));
  seed_lo[1] = std::max<int>(2*lo[1], std::ceil(get_ic_y_s(icv)*(2/dy) - .5));
  seed_lo[2] = std::max<int>(2*lo[2], std::ceil(get_ic_z_b(icv)*(2/dz) - .5));

  seed_hi[0] = std::min<int>(2*hi[0]+1, std::floor(get_ic_x_e(icv)*(2/dx) - .5));
  seed_hi[1] = std::min<int>(2*hi[1]+1, std::floor(get_ic_y_n(icv)*(2/dy) - .5));
  seed_hi[2] = std::min<int>(2*hi[2]+1, std::floor(get_ic_z_t(icv)*(2/dz) - .5));

  const Box bx(seed_lo, seed_hi);

  delta_bx[0] = std::max(0, seed_hi[0] - seed_lo[0] + 1);
  delta_bx[1] = std::max(0, seed_hi[1] - seed_lo[1] + 1);
  delta_bx[2] = std::max(0, seed_hi[2] - seed_lo[2] + 1);

  np = delta_bx[0] * delta_bx[1] * delta_bx[2];
  grow_pdata(pc + np);

  amrex::Real* p_rdata = m_rdata.data();

  const int local_nr = this->nr;

  amrex::ParallelFor(bx, [p_rdata,seed_lo,delta_bx,pc,dx,dy,dz,local_nr]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int local_i = i - seed_lo[0];
      const int local_j = j - seed_lo[1];
      const int local_k = k - seed_lo[2];

      const int local_pc =
        pc + (local_i + local_k*delta_bx[0] + local_j*delta_bx[0]*delta_bx[2]);

      p_rdata[local_pc*local_nr + 0] = (i + 0.5)*dx/2.;
      p_rdata[local_pc*local_nr + 1] = (j + 0.5)*dy/2.;
      p_rdata[local_pc*local_nr + 2] = (k + 0.5)*dz/2.;
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
ParticlesGenerator::random_fill (const int icv,
                                 const int type,
                                 const IntVect& lo,
                                 const IntVect& hi,
                                 int& np,
                                 int& pc,
                                 const amrex::Real dx,
                                 const amrex::Real dy,
                                 const amrex::Real dz,
                                 const bool fix_seed)
{
    // indices
  int i_w, i_e, j_s, j_n, k_b, k_t;

  const int maxfails = 1000;

  amrex::RealVect ic_dlo, ic_dhi, ic_len;
  amrex::Real max_dp, max_rp;
  
  calc_cell_ic(dx, dy, dz, get_ic_x_w(icv), get_ic_y_s(icv), get_ic_z_b(icv),
               get_ic_x_e(icv), get_ic_y_n(icv), get_ic_z_t(icv),
               i_w, i_e, j_s, j_n, k_b, k_t);

  // Start/end of IC domain bounds
  ic_dlo[0] = (std::max(lo[0], i_w)) * dx;
  ic_dlo[1] = (std::max(lo[1], j_s)) * dy;
  ic_dlo[2] = (std::max(lo[2], k_b)) * dz;

  ic_dhi[0] = (std::min(hi[0], i_e)+1) * dx;
  ic_dhi[1] = (std::min(hi[1], j_n)+1) * dy;
  ic_dhi[2] = (std::min(hi[2], k_t)+1) * dz;

  // physical volume of IC region intersecting this grid
  const amrex::Real ic_vol = (ic_dhi[0] - ic_dlo[0]) *
                             (ic_dhi[1] - ic_dlo[1]) *
                             (ic_dhi[2] - ic_dlo[2]);

  // Spacing is based on maximum particle size
  if(is_defined_db_cpp(get_ic_dp_max(icv,type)))
     max_dp = get_ic_dp_max(icv,type);
  else
     max_dp = get_ic_dp_mean(icv,type);

  // Particle count is based on mean particle size
  const amrex::Real mean = get_ic_dp_mean(icv,type);
  const int seed =
    static_cast<int>(ic_vol * get_ic_ep_s(icv,type) / ((M_PI/6.0)*mean*mean*mean));

  max_rp = 0.5 * max_dp;

  ic_len[0] = ic_dhi[0] - ic_dlo[0] - max_dp;
  ic_len[1] = ic_dhi[1] - ic_dlo[1] - max_dp;
  ic_len[2] = ic_dhi[2] - ic_dlo[2] - max_dp;

  ic_dlo[0] += max_rp;
  ic_dlo[1] += max_rp;
  ic_dlo[2] += max_rp;
  
  const amrex::Real mindist = (1.01*max_dp)*(1.01*max_dp);

  const RealVect Oodx = {1/dx, 1/dy, 1/dz};

  int nb(8);

  IntVect delta_bx;
  delta_bx[0] = hi[0] - lo[0] + 1;
  delta_bx[1] = hi[1] - lo[1] + 1;
  delta_bx[2] = hi[2] - lo[2] + 1;

  amrex::Vector<int> pinc(delta_bx[0]*delta_bx[1]*delta_bx[2], 0);
  amrex::Vector<int> pbin(delta_bx[0]*delta_bx[1]*delta_bx[2]*nb, 0);

  IntVect max_fitting_parts;
  max_fitting_parts[0] = static_cast<int>(std::ceil(ic_len[0]/max_dp));
  max_fitting_parts[1] = static_cast<int>(std::ceil(ic_len[1]/max_dp));
  max_fitting_parts[2] = static_cast<int>(std::ceil(ic_len[2]/max_dp));
  
  const int max_np = max_fitting_parts[0]*max_fitting_parts[1]*max_fitting_parts[2];
  grow_pdata(max_np);

  amrex::ResetRandomSeed(ParallelDescriptor::MyProc()+1);

  amrex::Real* p_rdata = m_rdata.data();

  np = 0;
  int iterations(0);

  while(iterations < seed)
  {
    bool stop = false;
    int fails(0);

    while(fails < maxfails and (not stop)) {
      RealVect pos;
      pos[0] = ic_dlo[0] + ic_len[0]*amrex::Random();
      pos[1] = ic_dlo[1] + ic_len[1]*amrex::Random();
      pos[2] = ic_dlo[2] + ic_len[2]*amrex::Random();

      // Grid containing the new particle
      IntVect idx;
      idx[0] = std::floor(pos[0]*Oodx[0]);
      idx[1] = std::floor(pos[1]*Oodx[1]);
      idx[2] = std::floor(pos[2]*Oodx[2]);

      // Local grid search for collisions.
      int overlaps = 0;

      IntVect seed_lo, seed_hi;
      seed_lo[0] = std::max(lo[0], idx[0]-1);
      seed_lo[1] = std::max(lo[1], idx[1]-1);
      seed_lo[2] = std::max(lo[2], idx[2]-1);

      seed_hi[0] = std::min(hi[0], idx[0]+1);
      seed_hi[1] = std::min(hi[1], idx[1]+1);
      seed_hi[2] = std::min(hi[2], idx[2]+1);

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
              const int local_pc =
                pbin[local_i + local_j*delta_bx[0] +
                     local_k*delta_bx[0]*delta_bx[1] +
                     p*delta_bx[0]*delta_bx[1]*delta_bx[2]];

              const Real dx = p_rdata[local_pc*nr + 0] - pos[0];
              const Real dy = p_rdata[local_pc*nr + 1] - pos[1];
              const Real dz = p_rdata[local_pc*nr + 2] - pos[2];

              const amrex::Real dist = dx*dx + dy*dy + dz*dz;

              if(dist < mindist)
                overlaps++;
            }
          }
        }
      }

      if(overlaps == 0)
      {
        p_rdata[np*nr + 0] = pos[0];
        p_rdata[np*nr + 1] = pos[1];
        p_rdata[np*nr + 2] = pos[2];

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

  amrex::Real* p_rdata = m_rdata.data();
  int* p_idata = m_idata.data();

  const int local_nr = this->nr;
  const int local_ni = this->ni;

  amrex::ParallelFor(nrp, [pstruct,p_rdata,p_idata,local_nr,local_ni]
    AMREX_GPU_DEVICE (int p) noexcept
  {
    ParticleType& part = pstruct[p];

    part.pos(0) = p_rdata[p*local_nr + 0];
    part.pos(1) = p_rdata[p*local_nr + 1];
    part.pos(2) = p_rdata[p*local_nr + 2];

    amrex::Real rad = p_rdata[p*local_nr + 3];
    amrex::Real rho = p_rdata[p*local_nr + 4];

    part.rdata(realData::velx) = p_rdata[p*local_nr + 5];
    part.rdata(realData::vely) = p_rdata[p*local_nr + 6];
    part.rdata(realData::velz) = p_rdata[p*local_nr + 7];

    amrex::Real vol  = (4.0/3.0)*M_PI*rad*rad*rad;
    amrex::Real mass = vol * rho;
    amrex::Real omoi = 2.5/(mass * rad*rad);

    part.rdata(realData::radius) = rad;
    part.rdata(realData::density) = rho;

    part.rdata(realData::volume) = vol;
    part.rdata(realData::mass) = mass;
    part.rdata(realData::oneOverI) = omoi;

    part.rdata(realData::omegax) = 0.0;
    part.rdata(realData::omegay) = 0.0;
    part.rdata(realData::omegaz) = 0.0;

    part.rdata(realData::dragx) = 0.0;
    part.rdata(realData::dragy) = 0.0;
    part.rdata(realData::dragz) = 0.0;

    part.idata(intData::phase) = p_idata[p*local_ni + 0];
    part.idata(intData::state) = 1;
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
ParticlesGenerator::nor_rno (amrex::Gpu::ManagedVector<amrex::Real>& dp,
                             const amrex::Real mean,
                             const amrex::Real sigma,
                             const amrex::Real dp_min,
                             const amrex::Real dp_max)
{
  // Local variables
  //-----------------------------------------------
  bool debug = false;
  //-----------------------------------------------

  const unsigned int nsize = dp.size();

  const int nsize_half = std::ceil(nsize/2.);

  const amrex::Real tolerance =
    std::numeric_limits<amrex::Real>::epsilon();

  amrex::Real* p_dp = dp.data();

  const int maxfails = 10000;
  int fails(0);

#ifdef AMREX_USE_CUDA
  Gpu::DeviceScalar<int> fails_gpu(fails);
  int* p_fails = fails_gpu.dataPtr();
#endif

  amrex::ParallelFor(nsize_half,
    [p_dp,dp_min,dp_max,maxfails,tolerance,sigma,mean,nsize,
#ifdef AMREX_USE_CUDA
     p_fails]
#else
     &fails]
#endif
    AMREX_GPU_DEVICE (int i) noexcept
    {
      amrex::Real x(0);
      amrex::Real y(0);

      amrex::Real dp1 = dp_min - 1;
      amrex::Real dp2 = dp_min - 1;

      int iterations(0);

      while(not(dp1 >= dp_min and dp1 <= dp_max and dp2 >= dp_min and dp2 <= dp_max) and
            iterations < maxfails)
      {
        amrex::Real w(1.1);
        while(w > 1 or std::abs(w-1) < tolerance)
        {
          x = 2*amrex::Random() - 1;
          y = 2*amrex::Random() - 1;
          w = x*x + y*y;
        }

        w = std::sqrt((-2.0 * std::log(w)) / w);

        amrex::Real dp1 = x * w * sigma + mean;
        amrex::Real dp2 = y * w * sigma + mean;

        p_dp[2*i] = dp1;

        if((2*i+1) < nsize)
          p_dp[2*i+1] = dp2;

        iterations++;
      }
      
      if(not(iterations < maxfails))
      {
#ifdef AMREX_USE_CUDA
        Gpu::Atomic::Add(p_fails, 1);
#else
        fails++;
#endif
      }
    });

  Gpu::synchronize();

#ifdef AMREX_USE_CUDA
  fails = fails_gpu.dataValue();
#endif

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fails == 0,
      "In ParticlesGenerator: nor_rno routine failed.");

  if(debug)
  {
    amrex::Real lmean(0), lvariance(0), lsigma(0);
    std::accumulate(dp.begin(), dp.end(), lmean);

    lmean /= nsize;

    for(int i = 0; i < nsize; i++)
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
ParticlesGenerator::uni_rno (amrex::Gpu::ManagedVector<amrex::Real>& dp,
                             const amrex::Real dp_min,
                             const amrex::Real dp_max)
{
  amrex::Real lscale = dp_max - dp_min;

  const unsigned int nsize = dp.size();

  amrex::Real* p_dp = dp.data();

  amrex::ParallelFor(nsize, [p_dp,dp_min,lscale]
    AMREX_GPU_DEVICE (int lc) noexcept { p_dp[lc] = dp_min + lscale*amrex::Random(); });

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
  const int r_gsize = gsize*nr;
  const int i_gsize = gsize*ni;

  // Increase real data
  if(m_rdata.size() == 0)
  {
    const int size = std::max<int>(r_gsize, 1024*nr);

    m_rdata.resize(size, 0);
  }
  else if(r_gsize >= m_rdata.size())
  {
    const int size = std::max<int>(2*m_rdata.size(), r_gsize);

    m_rdata.resize(size, 0);
  }

  // Increase integer data
  if(m_idata.size() == 0)
  {
    const int size = std::max<int>(i_gsize, 1024*ni);

    m_idata.resize(size, 0);
  }
  else if(i_gsize >= m_idata.size())
  {
    const int size = std::max<int>(2*m_idata.size(), i_gsize);

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
                           const ParticleTileType& particles,
                           const int nstep) const
{
  const auto& aos = particles.GetArrayOfStructs();
  const ParticleType* pstruct = aos().dataPtr();

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
    const ParticleType& part = pstruct[lc1];

    output_file << "               "
                << std::scientific << std::setw(13) << std::setprecision(6)
                << std::setfill(' ') << amrex::Real(part.rdata(realData::radius)) << std::endl;
  }

  output_file << "            " << "</DataArray>" << std::endl;

  output_file << "            "
    << "<DataArray type=\"Float32\" Name=\"density\" NumberOfComponents=\"1\" format=\"ascii\">"
    << std::endl;

  for(int lc1 = 0; lc1 < nrp; ++lc1)
  {
    const ParticleType& part = pstruct[lc1];

    output_file << "               "
                << std::scientific << std::setw(13) << std::setprecision(6)
                << amrex::Real(part.rdata(realData::density)) << std::endl;
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
                << amrex::Real(part.pos(0)) << "    "
                << amrex::Real(part.pos(1)) << "    "
                << amrex::Real(part.pos(2)) << std::endl;
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
