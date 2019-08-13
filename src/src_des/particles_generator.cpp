#include <particles_generator.hpp>
#include <AMReX_AmrParGDB.H>
#include <param_mod_F.H>
#include <ic_mod_F.H>
#include <discretelement_mod_F.H>
#include <calc_cell_F.H>
#include <random_nb_mod_F.H>

#include <limits>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>

namespace generator_aux {

struct ParticleType {
  amrex::Real pos[3];   // < Position -- components 1,2,3
  amrex::Real radius;   // < Radius   -- component  4
  amrex::Real volume;   // < Volume   -- component  5
  amrex::Real mass;     // < Mass     -- component  6
  amrex::Real density;  // < Density  -- component  7
  amrex::Real omoi;     // < One over momentum of inertia -- component 8
  amrex::Real vel[3];   // < Linear velocity              -- components 9,10,11
  amrex::Real omega[3]; // < Angular velocity             -- components 12,13,14
  amrex::Real drag[3];  // < Drag                         -- components 15,16,17
  int id;
  int cpu;
  int phase;
  int state;
};

static const amrex::Real sqrt3 = std::sqrt(3.0);
static const amrex::Real sqrt6o3x2 = 2.0*std::sqrt(6.0)/3.0;

} // end namespace generator_aux

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
using namespace generator_aux;

ParticlesGenerator::ParticlesGenerator()
  : m_rdata(0)
  , m_idata(0)
{}

ParticlesGenerator::~ParticlesGenerator()
{}

void
ParticlesGenerator::generate(int& pc,
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
  for(; icv0 <= get_dim_ic(); icv0++) {

    if(ic_defined_cpp(icv0) and std::abs(get_ic_ep_g(icv0)-1.0) > tolerance) {

      // Get the solids type index
      int type0(1);
      for(; type0 <= get_particle_types(); type0++)
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
      if(np0 > 0){

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

  amrex::Cuda::ManagedVector<amrex::Real> dp(np, 0);
  amrex::Cuda::ManagedVector<amrex::Real> ro_s(np, 0);

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

    AMREX_HOST_DEVICE_FOR_1D(np, p, { p_dp[p] = ic_dp_mean; });
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

    AMREX_HOST_DEVICE_FOR_1D(np, p, { p_ro_s[p] = ic_ro_s_mean; });
  }

  pc = init_pc;

  const amrex::Real ic_u_s = get_ic_u_s(icv, type);
  const amrex::Real ic_v_s = get_ic_v_s(icv, type);
  const amrex::Real ic_w_s = get_ic_w_s(icv, type);

  amrex::Real* p_rdata = m_rdata.data();
  int* p_idata = m_idata.data();

  AMREX_HOST_DEVICE_FOR_1D(np, p,
  {
    // amrex::Real pvol = (M_PI/6.0) * (dp[p]*dp[p]*dp[p]); // UNUSED_VARIABLE

    const int local_pc = pc + p;

    //< Radius................. 4
    p_rdata[local_pc*nr + 3] = 0.5 * p_dp[p];

    //< Density................ 5
    p_rdata[local_pc*nr + 4] = p_ro_s[p];

    //< Linear velocity........ 6,7,8
    p_rdata[local_pc*nr + 5] = ic_u_s;
    p_rdata[local_pc*nr + 6] = ic_v_s;
    p_rdata[local_pc*nr + 7] = ic_w_s;

    //< Type................... 1
    p_idata[local_pc*ni + 0] = type;
  });

  pc += np;

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
ParticlesGenerator::hex_close_pack(const int icv,
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

  amrex::Real ic_dlo[3], ic_dhi[3];
  amrex::Real max_dp, max_rp;

  int seed, max_seed[3], seed_lo[3], seed_hi[3];

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
  seed = int(ic_vol * get_ic_ep_s(icv,type) / ((M_PI/6.0)*mean*mean*mean));

  // Total to seed over the whole IC region
  max_seed[1] = int((get_ic_y_n(icv) - get_ic_y_s(icv) - max_dp) / max_dp);
  max_seed[2] = int((get_ic_z_t(icv) - get_ic_z_b(icv) - max_dp) / (sqrt3*max_rp));
  max_seed[0] = int(seed / (max_seed[1]*max_seed[2]));

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

  const Box bx(IntVect(seed_lo[0], seed_lo[1], seed_lo[2]),
               IntVect(seed_hi[0], seed_hi[1], seed_hi[2]));

  const int delta_bx_x = std::max(0, seed_hi[0] - seed_lo[0] + 1);
  const int delta_bx_y = std::max(0, seed_hi[1] - seed_lo[1] + 1);
  const int delta_bx_z = std::max(0, seed_hi[2] - seed_lo[2] + 1);

  np = delta_bx_x * delta_bx_y * delta_bx_z;
  grow_pdata(pc + np);

  amrex::Real* p_rdata = m_rdata.data();

  const amrex::Real sqrt3_copy(sqrt3);
  const amrex::Real sqrt6o3x2_copy(sqrt6o3x2);

  const amrex::Real seed_lo_x = seed_lo[0];
  const amrex::Real seed_lo_y = seed_lo[1];
  const amrex::Real seed_lo_z = seed_lo[2];

  AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
  {
    const int local_i = i - seed_lo_x;
    const int local_j = j - seed_lo_y;
    const int local_k = k - seed_lo_z;

    const int local_pc =
      pc + (local_j + local_k*delta_bx_y + local_i*delta_bx_y*delta_bx_z);

    p_rdata[local_pc*nr + 0] = i_w*dx + max_rp*(1. + i*sqrt6o3x2_copy);
    p_rdata[local_pc*nr + 1] = j_s*dy + max_rp*(1. + 2.*j + ((i+k)%2));
    p_rdata[local_pc*nr + 2] = k_b*dz + max_rp*(1. + sqrt3_copy*(k+((i%2)/3.)));
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
ParticlesGenerator::one_per_fill(const int icv,
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

  amrex::Real ic_dlo[3], ic_dhi[3];
  amrex::Real max_dp, max_rp;

  int seed, max_seed[3], seed_lo[3], seed_hi[3];

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
  seed = int(ic_vol * get_ic_ep_s(icv,type) / ((M_PI/6.0)*mean*mean*mean));

  // Total to seed over the whole IC region
  max_seed[0] = int((get_ic_x_e(icv) - get_ic_x_w(icv) - max_dp) / max_dp);
  max_seed[2] = int((get_ic_z_t(icv) - get_ic_z_b(icv) - max_dp) / (sqrt3*max_rp));
  max_seed[1] = int(seed / (max_seed[0]*max_seed[2]));

  // local grid seed loop hi/lo
  seed_lo[0] = std::round((ic_dlo[0] - i_w*dx) / max_dp);
  seed_lo[1] = std::round((ic_dlo[1] - j_s*dy) / ((sqrt6o3x2) * max_rp));
  seed_lo[2] = std::round((ic_dlo[2] - k_b*dz) / (sqrt3 * max_rp));

  seed_hi[0] = std::round((ic_dhi[0] - i_w*dx) /  max_dp - seed_lo[0]*max_dp);
  seed_hi[1] = std::round((ic_dhi[1] - j_s*dy) / ((sqrt6o3x2) * max_rp) - seed_lo[0]*max_dp);
  seed_hi[2] = std::round((ic_dhi[2] - k_b*dz) / (sqrt3 * max_rp) - seed_lo[0]*max_dp);

  seed_hi[0] = std::min(max_seed[0], seed_hi[0]-1);
  seed_hi[1] = std::min(max_seed[1], seed_hi[1]-1);
  seed_hi[2] = std::min(max_seed[2], seed_hi[2]-1);

  const int lo_x = std::max<int>(lo[0], std::ceil(get_ic_x_w(icv)/dx - .5));
  const int lo_y = std::max<int>(lo[1], std::ceil(get_ic_y_s(icv)/dy - .5));
  const int lo_z = std::max<int>(lo[2], std::ceil(get_ic_z_b(icv)/dz - .5));

  const int hi_x = std::min<int>(hi[0], std::floor(get_ic_x_e(icv)/dx - .5));
  const int hi_y = std::min<int>(hi[1], std::floor(get_ic_y_n(icv)/dy - .5));
  const int hi_z = std::min<int>(hi[2], std::floor(get_ic_z_t(icv)/dz - .5));

  const Box bx(IntVect(lo_x, lo_y, lo_z), IntVect(hi_x, hi_y, hi_z));

  const int delta_bx_x = std::max(0, hi_x - lo_x + 1);
  const int delta_bx_y = std::max(0, hi_y - lo_y + 1);
  const int delta_bx_z = std::max(0, hi_z - lo_z + 1);

  np = delta_bx_x * delta_bx_y * delta_bx_z;
  grow_pdata(pc + np);

  amrex::Real* p_rdata = m_rdata.data();

  AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
  {
    const int local_i = i - lo_x;
    const int local_j = j - lo_y;
    const int local_k = k - lo_z;

    const int local_pc =
      pc + (local_i + local_k*delta_bx_x + local_j*delta_bx_x*delta_bx_z);

    p_rdata[local_pc*nr + 0] = (i + 0.5) * dx;
    p_rdata[local_pc*nr + 1] = (j + 0.5) * dy;
    p_rdata[local_pc*nr + 2] = (k + 0.5) * dz;
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
ParticlesGenerator::eight_per_fill(const int icv,
                                   const int type,
                                   const IntVect& lo,
                                   const IntVect& hi,
                                   int& np,
                                   int& pc,
                                   const amrex::Real dx,
                                   const amrex::Real dy,
                                   const amrex::Real dz)
{
  const int lo_x = std::max<int>(2*lo[0], std::ceil(get_ic_x_w(icv)*(2/dx) - .5));
  const int lo_y = std::max<int>(2*lo[1], std::ceil(get_ic_y_s(icv)*(2/dy) - .5));
  const int lo_z = std::max<int>(2*lo[2], std::ceil(get_ic_z_b(icv)*(2/dz) - .5));

  const int hi_x = std::min<int>(2*hi[0]+1, std::floor(get_ic_x_e(icv)*(2/dx) - .5));
  const int hi_y = std::min<int>(2*hi[1]+1, std::floor(get_ic_y_n(icv)*(2/dy) - .5));
  const int hi_z = std::min<int>(2*hi[2]+1, std::floor(get_ic_z_t(icv)*(2/dz) - .5));

  const Box bx(IntVect(lo_x, lo_y, lo_z), IntVect(hi_x, hi_y, hi_z));

  const int delta_bx_x = std::max(0, hi_x - lo_x + 1);
  const int delta_bx_y = std::max(0, hi_y - lo_y + 1);
  const int delta_bx_z = std::max(0, hi_z - lo_z + 1);

  np = delta_bx_x * delta_bx_y * delta_bx_z;
  grow_pdata(pc + np);

  amrex::Real* p_rdata = m_rdata.data();

  AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
  {
    const int local_i = i - lo_x;
    const int local_j = j - lo_y;
    const int local_k = k - lo_z;

    const int local_pc =
      pc + (local_i + local_k*delta_bx_x + local_j*delta_bx_x*delta_bx_z);

    p_rdata[local_pc*nr + 0] = (i + 0.5)*dx/2.;
    p_rdata[local_pc*nr + 1] = (j + 0.5)*dy/2.;
    p_rdata[local_pc*nr + 2] = (k + 0.5)*dz/2.;
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
ParticlesGenerator::random_fill(const int icv,
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

  amrex::Real ic_dlo[3], ic_dhi[3], ic_len[3];
  amrex::Real max_dp, max_rp;

  amrex::Vector<amrex::Real> pos(3, -1.e20);

  const amrex::Real Oodx[3] = {1/dx, 1/dy, 1/dz};

  int seed;
  // int max_seed[3], seed_lo[3], seed_hi[3]; // UNUSED_VARIABLE

  const int maxfails = 1000;

  int overlaps(0);

  int fails(0), nb(8);

  const unsigned x_dim = hi[0] - lo[0] + 1;
  const unsigned y_dim = hi[1] - lo[1] + 1;
  const unsigned z_dim = hi[2] - lo[2] + 1;

  amrex::Vector<int> pinc(x_dim*y_dim*z_dim, 0);
  amrex::Vector<int> pbin(x_dim*y_dim*z_dim*nb, 0);

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

  max_rp = 0.5 * max_dp;

  // Particle count is based on mean particle size
  const amrex::Real mean = get_ic_dp_mean(icv,type);
  seed = int(ic_vol * get_ic_ep_s(icv,type) / ((M_PI/6.0)*mean*mean*mean));

  ic_len[0] = ic_dhi[0] - ic_dlo[0] - max_dp;
  ic_dlo[0] += max_rp;
  ic_len[1] = ic_dhi[1] - ic_dlo[1] - max_dp;
  ic_dlo[1] += max_rp;
  ic_len[2] = ic_dhi[2] - ic_dlo[2] - max_dp;
  ic_dlo[2] += max_rp;

  np = 0;

  const amrex::Real mindist = (1.01*max_dp)*(1.01*max_dp);

  while(np < seed and fails < maxfails)
  {
    bool stop = false;

    do {
      pos[0] = ic_dlo[0] + ic_len[0]*get_random();
      pos[1] = ic_dlo[1] + ic_len[1]*get_random();
      pos[2] = ic_dlo[2] + ic_len[2]*get_random();

      // Grid containing the new particle
      const int idx_x = std::floor(pos[0]*Oodx[0]);
      const int idx_y = std::floor(pos[1]*Oodx[1]);
      const int idx_z = std::floor(pos[2]*Oodx[2]);

      // Local grid search for collisions.
      overlaps = 0;

      const int lo_x = std::max(lo[0], idx_x-1);
      const int lo_y = std::max(lo[1], idx_y-1);
      const int lo_z = std::max(lo[2], idx_z-1);

      const int hi_x = std::min(hi[0], idx_x+1);
      const int hi_y = std::min(hi[1], idx_y+1);
      const int hi_z = std::min(hi[2], idx_z+1);

      for(int k(lo_z); k <= hi_z; ++k)
        for(int j(lo_y); j <= hi_y; ++j)
          for(int i(lo_x); i <= hi_x; ++i)
          {
            const int local_i = i - lo[0];
            const int local_j = j - lo[1];
            const int local_k = k - lo[2];

            const int upper_limit =
              pinc[local_i + local_j*x_dim + local_k*x_dim*y_dim];

            for(int l = 0; l < upper_limit; l++)
            {
              const int local_pc =
                pbin[local_i + local_j*x_dim + local_k*x_dim*y_dim + l*x_dim*y_dim*z_dim];

              amrex::Real x = m_rdata[local_pc*nr + 0] - pos[0];
              amrex::Real y = m_rdata[local_pc*nr + 1] - pos[1];
              amrex::Real z = m_rdata[local_pc*nr + 2] - pos[2];

              const amrex::Real dist = x*x + y*y + z*z;

              if(dist < mindist)
                overlaps++;
            }
          }

      if(overlaps == 0)
      {
        np++; // local to call
        pc++; // local to grid

        grow_pdata(pc);

        const int np_idx = np-1;
        const int pc_idx = pc-1;

        m_rdata[pc_idx*nr + 0] = pos[0];
        m_rdata[pc_idx*nr + 1] = pos[1];
        m_rdata[pc_idx*nr + 2] = pos[2];

        const int local_idx_x = idx_x - lo[0];
        const int local_idx_y = idx_y - lo[1];
        const int local_idx_z = idx_z - lo[2];

        const int local_idx =
          local_idx_x + local_idx_y*x_dim + local_idx_z*x_dim*y_dim;

        const int pinc_value = pinc[local_idx];
        pbin[local_idx + pinc_value*x_dim*y_dim*z_dim] = np_idx;

        pinc[local_idx] = pinc_value + 1;

        if(pinc_value + 2 >= nb)
        {
          nb += 2;
          pbin.resize(x_dim*y_dim*z_dim*nb, 0);
        }

        fails = 0;
        stop = true;
      }
      else
      {
        fails++;
      }
    } while(not stop);

  //  if((mod(np, seed/10) == 0 .and. np < seed*0.95) .or. np==seed) &
  //       write(*,"(2x,'Seeded: ',I9,3x,'(',f5.0,'%)')") np,100*dble(np)/seed

  }

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
ParticlesGenerator::generate_prop(const int nrp,
                                  void* particles)
{
  generator_aux::ParticleType* parts =
    static_cast<generator_aux::ParticleType*>(particles);

  amrex::Real* p_rdata = m_rdata.data();
  int* p_idata = m_idata.data();

  AMREX_HOST_DEVICE_FOR_1D(nrp, p,
  {
    parts[p].pos[0] = p_rdata[p*nr + 0];
    parts[p].pos[1] = p_rdata[p*nr + 1];
    parts[p].pos[2] = p_rdata[p*nr + 2];

    amrex::Real rad = p_rdata[p*nr + 3];
    amrex::Real rho = p_rdata[p*nr + 4];

    parts[p].vel[0] = p_rdata[p*nr + 5];
    parts[p].vel[1] = p_rdata[p*nr + 6];
    parts[p].vel[2] = p_rdata[p*nr + 7];

    amrex::Real vol  = (4.0/3.0)*M_PI*rad*rad*rad;
    amrex::Real mass = vol * rho;
    amrex::Real omoi = 2.5/(mass * rad*rad);

    parts[p].radius  = rad;
    parts[p].density = rho;

    parts[p].volume  = vol;
    parts[p].mass    = mass;
    parts[p].omoi    = omoi;

    parts[p].omega[0] = 0.0;
    parts[p].omega[1] = 0.0;
    parts[p].omega[2] = 0.0;

    parts[p].drag[0]  = 0.0;
    parts[p].drag[1]  = 0.0;
    parts[p].drag[2]  = 0.0;

    parts[p].phase = p_idata[p*ni + 0];
    parts[p].state = 1;
  });

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
ParticlesGenerator::nor_rno(amrex::Cuda::ManagedVector<amrex::Real>& dp,
                            const amrex::Real mean,
                            const amrex::Real sigma,
                            const amrex::Real dp_min,
                            const amrex::Real dp_max)
{
  // Local variables
  //-----------------------------------------------
  amrex::Real x[2];
  bool debug = false;
  //-----------------------------------------------

  const unsigned int nsize = dp.size();
  int i(0);

  amrex::Real* p_dp = dp.data();

  const amrex::Real tolerance =
    std::numeric_limits<amrex::Real>::epsilon();

  while(i < std::ceil(nsize/2.))
  {
    amrex::Real w(1);
    while(w > 1 or std::abs(w-1) < tolerance)
    {
      x[0] = 2*amrex::Random() - 1;
      x[1] = 2*amrex::Random() - 1;
      w = x[0]*x[0] + x[1]*x[1];
    }

    w = std::sqrt((-2.0 * std::log(w)) / w);

    amrex::Real dp1 = x[0] * w * sigma + mean;
    amrex::Real dp2 = x[1] * w * sigma + mean;

    if(dp1 >= dp_min and dp2 >= dp_min and dp1 <= dp_max and dp2 <= dp_max)
    {
      p_dp[2*i] = dp1;

      if((2*i+1) < nsize)
        p_dp[2*i+1] = dp2;

      i++;
    }
  }

  if(debug)
  {
    amrex::Real lmean(0), lvariance(0), lsigma(0);
    // TODO: parallelization
    std::accumulate(dp.begin(), dp.end(), lmean);

    lmean /= nsize;

    for(i = 0; i < nsize; i++)
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
ParticlesGenerator::uni_rno(amrex::Cuda::ManagedVector<amrex::Real>& dp,
                            const amrex::Real dp_min,
                            const amrex::Real dp_max)
{
  // TODO: can we do this on GPU? Is there a CURand generator in AMReX?
  std::fill(dp.begin(), dp.end(), amrex::Random());

  amrex::Real* p_dp = dp.data();

  amrex::Real lscale = dp_max - dp_min;

  const unsigned int nsize = dp.size();

  AMREX_HOST_DEVICE_FOR_1D(nsize, lc,
  {
    p_dp[lc] = dp_min + lscale*p_dp[lc];
  });

  return;
}

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
//                                                                     !
//                                                                     !
//                                                                     !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
void
ParticlesGenerator::grow_pdata(const int gsize)
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
ParticlesGenerator::write(const int nrp,
                          void* particles) const
{
  const generator_aux::ParticleType* parts =
    static_cast<const generator_aux::ParticleType*>(particles);

  std::ofstream output_file;
  output_file.open("test.vtp");

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

  // TODO openMP
  for(int lc1 = 0; lc1 < nrp; ++lc1)
  {
    output_file << "               "
                << std::scientific << std::setw(13) << std::setprecision(6)
                << std::setfill(' ') << amrex::Real(parts[lc1].radius) << std::endl;
  }

  output_file << "            " << "</DataArray>" << std::endl;

  output_file << "            "
    << "<DataArray type=\"Float32\" Name=\"density\" NumberOfComponents=\"1\" format=\"ascii\">"
    << std::endl;

  // TODO openMP
  for(int lc1 = 0; lc1 < nrp; ++lc1)
  {
    output_file << "               "
                << std::scientific << std::setw(13) << std::setprecision(6)
                << amrex::Real(parts[lc1].density) << std::endl;
  }

  output_file << "            " << "</DataArray>" << std::endl;

  output_file << "         " << "</PointData>" << std::endl;

  output_file << "         " << "<Points>" << std::endl;

  output_file << "            " << "<DataArray type=\"Float32\" "
              << "Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">"
              << std::endl;

  // TODO openMP
  for(int lc1 = 0; lc1 < nrp; lc1++)
  {
    output_file << "               "
                << std::scientific << std::setw(13) << std::setprecision(6)
                << amrex::Real(parts[lc1].pos[0]) << "    "
                << amrex::Real(parts[lc1].pos[1]) << "    "
                << amrex::Real(parts[lc1].pos[2]) << std::endl;
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
