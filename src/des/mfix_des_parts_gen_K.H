#ifndef _MFIX_DES_PARTS_GEN_K_H_
#define _MFIX_DES_PARTS_GEN_K_H_

#include <mfix_calc_cell.H>

#include <limits>

#include <mfix_ic.H>

#include <AMReX_IntVect.H>
#include <AMReX_RealVect.H>
#include <AMReX_PODVector.H>

// This is a very small number we add to the length to catch
// any particles that are "on the edge"
#define small_number 10.*std::numeric_limits<amrex::Real>::epsilon()

#define loc_tolerance std::numeric_limits<amrex::Real>::epsilon()

#define sqrt6x2o3 std::sqrt(6.0)*(2.0/3.0)
#define sqrt3o3 std::sqrt(3.0)/3.0
#define sqrt3 std::sqrt(3.0)

#define DEBUG_HCP 0

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
class Hex_ClosePack
{
  public:
    // Constructor
    AMREX_GPU_HOST_DEVICE
    Hex_ClosePack ()
      : m_is_initialized(0)
      , m_plo(0.)
      , m_dx(0.)
      , m_np(0)
      , m_p1r1_np(0)
      , m_p2r1_np(0)
      , m_p1_np(0)
      , m_p1p2_np(0)
      , m_r1r2_np(0)
      , m_klo(0)
      , m_p1_jlo(0)
      , m_p1r1_ilo(0)
      , m_p1r2_ilo(0)
      , m_p2_jlo(0)
      , m_p2r1_ilo(0)
      , m_p2r2_ilo(0)
      , m_ic_dlo(0.)
      , m_eff_rad(0)
      , m_rbx_lo(0.)
      , m_rbx_hi(0.)
      , m_verbose(0)
    {}


    AMREX_GPU_HOST_DEVICE
    Hex_ClosePack (const amrex::RealVect& plo,
                   const amrex::RealVect& dx)
      : m_is_initialized(1)
      , m_plo(plo)
      , m_dx(dx)
      , m_np(0)
      , m_p1r1_np(0)
      , m_p2r1_np(0)
      , m_p1_np(0)
      , m_p1p2_np(0)
      , m_r1r2_np(0)
      , m_klo(0)
      , m_p1_jlo(0)
      , m_p1r1_ilo(0)
      , m_p1r2_ilo(0)
      , m_p2_jlo(0)
      , m_p2r1_ilo(0)
      , m_p2r2_ilo(0)
      , m_ic_dlo(0.)
      , m_eff_rad(0)
      , m_rbx_lo(0.)
      , m_rbx_hi(0.)
      , m_verbose(0)
    {}

    AMREX_GPU_HOST_DEVICE
    void initialize (const amrex::RealVect& plo,
                     const amrex::RealVect& dx)
    {
      m_is_initialized = 1;

      for (int dir(0); dir < 3; ++dir) {
        m_plo[dir] = plo[dir];
        m_dx[dir] = dx[dir];
      }
    }

    AMREX_GPU_HOST_DEVICE
    void setup (const amrex::IntVect& lo,
                const amrex::IntVect& hi,
                const amrex::RealBox& region,
                const amrex::Real mean,
                const amrex::Real eps,
                const amrex::Real eps_tolerance = 1.e-15)
    {
#if 0
      amrex::ParmParse pp("ic");

      std::vector<std::string> input_regions;
      pp.queryarr("regions", input_regions);

      m_verbose = (input_regions[icv] == "init-bed-FR") ? 100 : 0;
#endif


      // physical volume of IC region
      const amrex::Real vol = region.volume();

      m_np = 0;

      if (!m_is_initialized) {
        return;
      }

      if (!(eps > eps_tolerance && eps < 1.)) {
        return;
      }

      // Total particle count is based on mean particle size
      const long seed =
        static_cast<long>(vol * eps / ((M_PI/6.0)*mean*mean*mean));

      if (seed == 0) {
        return;
      }

      m_ic_dlo = amrex::RealVect(region.lo(0), region.lo(1), region.lo(2));

      amrex::RealVect ic_dhi(region.hi(0), region.hi(1), region.hi(2));

      const amrex::Real xlen = ic_dhi[0] - m_ic_dlo[0];
      const amrex::Real ylen = ic_dhi[1] - m_ic_dlo[1];
      const amrex::Real zlen = ic_dhi[2] - m_ic_dlo[2];

      amrex::RealVect rbx_lo(m_plo[0] + m_dx[0]*static_cast<amrex::Real>(lo[0]),
                             m_plo[1] + m_dx[1]*static_cast<amrex::Real>(lo[1]),
                             m_plo[2] + m_dx[2]*static_cast<amrex::Real>(lo[2]));

      amrex::RealVect rbx_hi(m_plo[0] + m_dx[0]*static_cast<amrex::Real>(hi[0]+1),
                             m_plo[1] + m_dx[1]*static_cast<amrex::Real>(hi[1]+1),
                             m_plo[2] + m_dx[2]*static_cast<amrex::Real>(hi[2]+1));

#if 0
      amrex::PrintToFile("setup", amrex::Print::AllProcs) << "\n"
                     << " m_ic_dlo " << m_ic_dlo[0] << "  " << m_ic_dlo[1] << "  " << m_ic_dlo[2] << "\n"
                     << " m_ic_dhi " <<   ic_dhi[0] << "  " <<   ic_dhi[1] << "  " <<   ic_dhi[2] << "\n";

      amrex::PrintToFile("setup", amrex::Print::AllProcs) << "\n"
                     << " rbx_lo   " << rbx_lo[0] << "  " << rbx_lo[1] << "  " << rbx_lo[2] << "\n"
                     << " rbx_hi   " << rbx_hi[0] << "  " << rbx_hi[1] << "  " << rbx_hi[2] << "\n";
#endif

      const amrex::Real mean_rad = 0.5*mean;

      // Verify that the IC region overlaps with the IC region.
      if (ic_dhi[0] <= rbx_lo[0] + mean_rad || rbx_hi[0] <= m_ic_dlo[0] + mean_rad ||
          ic_dhi[1] <= rbx_lo[1] + mean_rad || rbx_hi[1] <= m_ic_dlo[1] + mean_rad ||
          ic_dhi[2] <= rbx_lo[2] + mean_rad || rbx_hi[2] <= m_ic_dlo[2] + mean_rad) {

        return;
      }

      // This is the radius we need for spacing particles to get the
      // desired volume fraction and particle count.
      m_eff_rad = hcp_rad(mean_rad, xlen, ylen, zlen, seed);
      //printf ("\nPacking radius: %21.18f\n", m_eff_rad);

#if 0
      amrex::PrintToFile("setup", amrex::Print::AllProcs)
                     << "\neff_rad " << m_eff_rad << "\n";
#endif

      // Pull back just less than one particle radius to keep particles
      // fully inside the IC region.
      ic_dhi -= (m_eff_rad - small_number);

      rbx_hi[0] = amrex::min(rbx_hi[0], ic_dhi[0]);
      rbx_hi[1] = amrex::min(rbx_hi[1], ic_dhi[1]);
      rbx_hi[2] = amrex::min(rbx_hi[2], ic_dhi[2]);

      m_rbx_hi = rbx_hi;
      m_rbx_lo = rbx_lo;

      // Transition into relative coordinates with origin at (0,0,0)
      rbx_lo -= m_ic_dlo;
      rbx_hi -= m_ic_dlo;

#if 0
      amrex::PrintToFile("setup", amrex::Print::AllProcs) << "\n"
                     << "shifted rbx: \n"
                     << " rbx_lo   " << rbx_lo[0] << "  " << rbx_lo[1] << "  " << rbx_lo[2] << "\n"
                     << " rbx_hi   " << rbx_hi[0] << "  " << rbx_hi[1] << "  " << rbx_hi[2] << "\n";
#endif

      // Compute the first plane.

      m_klo = hcp_k(m_eff_rad, rbx_lo[2]);        // k index of first plane
      const int p_cnt = hcp_k(m_eff_rad, rbx_hi[2], m_klo) - m_klo; // Total number of planes


      const int p2_cnt = (p_cnt - (p_cnt%2))/2;

      const int p1_cnt = p_cnt - p2_cnt;

      m_p1_jlo = hcp_j(m_eff_rad, rbx_lo[1], (m_klo%2));
      const int p1_rows  = hcp_j(m_eff_rad, rbx_hi[1], (m_klo%2), m_p1_jlo) - m_p1_jlo;


      const int p1r2_cnt = (p1_rows - (p1_rows%2))/2;
      const int p1r1_cnt = p1_rows - p1r2_cnt;

      m_p1r1_ilo = hcp_i(m_eff_rad, rbx_lo[0], ((m_klo+m_p1_jlo  )%2));
      m_p1r2_ilo = hcp_i(m_eff_rad, rbx_lo[0], ((m_klo+m_p1_jlo+1)%2));

      m_p1r1_np  = hcp_i(m_eff_rad, rbx_hi[0], ((m_klo+m_p1_jlo  )%2), m_p1r1_ilo) - m_p1r1_ilo;
      const int p1r2_np  = hcp_i(m_eff_rad, rbx_hi[0], ((m_klo+m_p1_jlo+1)%2), m_p1r2_ilo) - m_p1r2_ilo;

      m_p1_np = m_p1r1_np*p1r1_cnt + p1r2_np*p1r2_cnt;

      // Build out plane 2
      m_p2_jlo = hcp_j(m_eff_rad, rbx_lo[1], ((m_klo+1)%2));
      const int p2_rows  = hcp_j(m_eff_rad, rbx_hi[1], ((m_klo+1)%2), m_p2_jlo) - m_p2_jlo;


      const int p2r2_cnt = (p2_rows - (p2_rows%2))/2;
      const int p2r1_cnt = p2_rows - p2r2_cnt;

      m_p2r1_ilo = hcp_i(m_eff_rad, rbx_lo[0], ((m_klo+1+m_p2_jlo  )%2));
      m_p2r2_ilo = hcp_i(m_eff_rad, rbx_lo[0], ((m_klo+1+m_p2_jlo+1)%2));

      m_p2r1_np  = hcp_i(m_eff_rad, rbx_hi[0], ((m_klo+1+m_p2_jlo  )%2), m_p2r1_ilo) - m_p2r1_ilo;
      const int p2r2_np  = hcp_i(m_eff_rad, rbx_hi[0], ((m_klo+1+m_p2_jlo+1)%2), m_p2r2_ilo) - m_p2r2_ilo;

      const int p2_np = m_p2r1_np*p2r1_cnt + p2r2_np*p2r2_cnt;

      m_np = m_p1_np*p1_cnt + p2_np*p2_cnt;

      m_p1p2_np = m_p1_np + p2_np;

      m_r1r2_np = m_p1r1_np + p1r2_np;

#ifndef AMREX_USE_GPU
      if(m_verbose > 100) {
        amrex::AllPrint() << "rad: " << m_eff_rad << "   seed " << seed << "   m_np " << m_np << "\n";
      }
#endif


#if 0
      amrex::PrintToFile("setup", amrex::Print::AllProcs) << "\n"
                     << "m_klo           " << m_klo << "\n"
                     << "p_cnt           " << p_cnt << "\n"
                     << "p1_cnt          " << p1_cnt << "\n"
                     << "p2_cnt          " << p2_cnt << "\n"
                     << "m_p1_jlo        " << m_p1_jlo << "\n"
                     << "p1_rows         " << p1_rows << "\n"
                     << "p1r1_cnt        " << p1r1_cnt << "\n"
                     << "p1r2_cnt        " << p1r2_cnt << "\n"
                     << "m_p1r1_ilo      " << m_p1r1_ilo << "\n"
                     << "m_p1r2_ilo      " << m_p1r2_ilo << "\n"
                     << "m_p1_np         " << m_p1_np << "\n"
                     << "\n"
                     << "m_p2_jlo        " << m_p2_jlo << "\n"
                     << "p2_rows         " << p2_rows << "\n"
                     << "p2r1_cnt        " << p2r1_cnt << "\n"
                     << "p2r2_cnt        " << p2r2_cnt << "\n"
                     << "m_p2r1_ilo      " << m_p2r1_ilo << "\n"
                     << "m_p2r2_ilo      " << m_p2r2_ilo << "\n"
                     << "m_p2r1_np       " << m_p2r1_np << "\n"
                     << "  p2r2_np       " <<   p2r2_np << "\n"
                     << "p2_np           " << p2_np << "\n"
                     << "\n"
                     << "m_np            " << m_np      << "\n"
                     << "m_p1p2_np       " << m_p1p2_np << "\n"
                     << "m_r1r2_np       " << m_r1r2_np << "\n";
#endif

    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int get_particles_number () const
    { return m_np; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_particles_radius () const
    { return m_eff_rad; }

    template <amrex::RunOn >
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::RealVect operator() (const int p_id,
                                amrex::RandomEngine const&) const
    {

      // Compute the grouped (type-1 + type-2) plane containing the particle.
      const int lk0 = p_id/m_p1p2_np; // (int/int) no need for floor()
      // Compute the particle index in the grouped plane.
      const int lk1 = p_id - lk0*m_p1p2_np;
      // 0: particle is in type-1 plane
      // 1: particle is in type-2 plane
      const int lk2 = (lk1 < m_p1_np ? 0 : 1);

      // Particle index in the hcp on this grid.
      const int local_k = lk2 + 2*lk0;
      // Particle index in the full hcp.
      const int k = m_klo + local_k;

      // Particle index in the k-th plane.
      const int local_ko2 = local_k/2; // (int/int) no need for floor()
      const int lj0 = p_id - ((local_k%2)*m_p1_np + local_ko2*m_p1p2_np);
      // Index of grouped (type-1 + type-2) row in the k-th plane.
      const int lj1 = lj0/m_r1r2_np; // (int/int) not need for floor()
      // Particle index within the grouped row
      const int lj2 = lj0 - lj1*m_r1r2_np;
      // Number of particles in first row or this plane
      const int lj3 = (1-(local_k%2))*m_p1r1_np + (local_k%2)*m_p2r1_np;

      // Particle index in the hcp on this grid.
      // 0: particle is in the first row of the group
      // 1: particle is in the second row of the group
      const int local_j = (lj2 < lj3 ? 0 : 1) + 2*lj1;
      // Particle index in the full hcp.
      const int j = local_j + (1-local_k%2)*m_p1_jlo
                            + (  local_k%2)*m_p2_jlo;


      // Particle index in the hcp on this grid.
      const int li0 = (lj2 < lj3 ? 0 : 1)*lj3;
      const int local_i = lj2 - li0;
      // Particle index in the full hcp.
      const int i = local_i + (1-(local_k%2))*((1-(local_j%2))*m_p1r1_ilo  +
                                               (  (local_j%2))*m_p1r2_ilo) +
                              (  (local_k%2))*((1-(local_j%2))*m_p2r1_ilo  +
                                               (  (local_j%2))*m_p2r2_ilo);

      const amrex::Real x = m_ic_dlo[0] + m_eff_rad * (1.0 + static_cast<amrex::Real>(2*i + (j+k)%2));
      const amrex::Real y = m_ic_dlo[1] + m_eff_rad * (1.0 + sqrt3o3*(static_cast<amrex::Real>(3*j + k%2)));
      const amrex::Real z = m_ic_dlo[2] + m_eff_rad * (1.0 + sqrt6x2o3*(static_cast<amrex::Real>(k)));

#if 0
      {
         amrex::PrintToFile("generate", amrex::Print::AllProcs) << "p_id: " << p_id
                        << "  local: " << amrex::IntVect(local_i, local_j, local_k)
                        << "  global: " << amrex::IntVect(i,j,k)
                        << "  pos: " << x << "  " << y << "  " << z
                        << "  rbx: " << m_rbx_lo[0] << "  " << m_rbx_lo[1] << "  " << m_rbx_lo[2]
                        << "   hi: " << m_rbx_hi[0] << "  " << m_rbx_hi[1] << "  " << m_rbx_hi[2]
                        //<< "  ??  " << (1-((local_j+local_k)%2)) << "  " << m_lo0[0]
                        //<< "  OR  " << (  ((local_j+local_k)%2)) << "  " << m_lo1[0]
                        << "\n";
      }

#endif

      AMREX_ALWAYS_ASSERT(m_rbx_lo[0] <= x && x < m_rbx_hi[0]);
      AMREX_ALWAYS_ASSERT(m_rbx_lo[1] <= y && y < m_rbx_hi[1]);
      AMREX_ALWAYS_ASSERT(m_rbx_lo[2] <= z && z < m_rbx_hi[2]);


      return amrex::RealVect(x, y, z);
    }

  private:
    int m_is_initialized;

    amrex::RealVect m_plo;

    amrex::RealVect m_dx;

    int m_np;

    // Number of particles in plane type 1; row type 1
    int m_p1r1_np;

    // Number of particles in plane type 2; row type 1
    int m_p2r1_np;

    // Number of particles in plane type 1
    int m_p1_np;

    // Number of particles in combined type 1 + type 2 plane
    int m_p1p2_np;

    // Number of particles in combined type 1 + type 2 plane
    int m_r1r2_np;

    int m_klo;

    int m_p1_jlo;

    int m_p1r1_ilo;
    int m_p1r2_ilo;

    int m_p2_jlo;
    int m_p2r1_ilo;
    int m_p2r2_ilo;

    amrex::RealVect m_ic_dlo;

    amrex::Real m_eff_rad;

    amrex::RealVect m_rbx_lo;
    amrex::RealVect m_rbx_hi;

    int m_verbose = 0;

    // Auxiliary functions
    AMREX_GPU_HOST_DEVICE
    int hcp_i (amrex::Real rad,
               amrex::Real xend,
               int offset_int,
               int start_i = 0) const
    {
      const amrex::Real offset = static_cast<amrex::Real>(offset_int);
      int i(start_i);
      while ( loc_tolerance < xend - (1.0 + 2.0*static_cast<amrex::Real>(i) + offset)*rad ) {
        i++;
      }
      return i;
    }

    AMREX_GPU_HOST_DEVICE
    int hcp_j (amrex::Real rad,
               amrex::Real yend,
               int offset_int,
               int start_j = 0) const
    {
      const amrex::Real offset = static_cast<amrex::Real>(offset_int);
      int j(start_j);
      while ( loc_tolerance < yend - (1.0 + sqrt3o3 * (3.0 * static_cast<amrex::Real>(j) + offset))*rad ) {
        j++;
      }
      return j;
    }


    AMREX_GPU_HOST_DEVICE
    int hcp_k (amrex::Real rad,
               amrex::Real zend,
               int start_k = 0) const
    {
      int k(start_k);
      while ( loc_tolerance < zend - (1.0 + sqrt6x2o3*static_cast<amrex::Real>(k))*rad ) {
        k++;
      }
      return k;
    }

    AMREX_GPU_HOST_DEVICE
    long hcp_count (amrex::Real rad,
                    amrex::Real xlen,
                    amrex::Real ylen,
                    amrex::Real zlen) const
    {


      // Compute the first plane.

      const long p_cnt = static_cast<long>(hcp_k(rad, zlen)); // Total number of planes

      const long p2_cnt = (p_cnt - (p_cnt%2))/2;
      const long p1_cnt = p_cnt - p2_cnt;

      const long p1_rows = static_cast<long>(hcp_j(rad, ylen, 0));

      const long p1r2_cnt = (p1_rows - (p1_rows%2))/2;
      const long p1r1_cnt = p1_rows - p1r2_cnt;

      const long p1r1_np = static_cast<long>(hcp_i(rad, xlen, 0));
      const long p1r2_np = static_cast<long>(hcp_i(rad, xlen, 1));

      const long p1_np = p1r1_np*p1r1_cnt + p1r2_np*p1r2_cnt;

      // Build out plane 2
      const long p2_rows = static_cast<long>(hcp_j(rad, ylen, 1));

      const long p2r2_cnt = (p2_rows - (p2_rows%2))/2;
      const long p2r1_cnt = p2_rows - p2r2_cnt;

      const long p2r1_np  = static_cast<long>(hcp_i(rad, xlen, 1));
      const long p2r2_np  = static_cast<long>(hcp_i(rad, xlen, 0));

      const int p2_np = p2r1_np*p2r1_cnt + p2r2_np*p2r2_cnt;

      // Total number of particles
      const long np = p1_np*p1_cnt + p2_np*p2_cnt;

#ifndef AMREX_USE_GPU
      if(m_verbose >= 100) {
        amrex::PrintToFile("hcp_rad", amrex::Print::AllProcs)
          << p1_np << "  "
          << p1_cnt << "  "
          << p2_np << "  "
          << p2_cnt << "  "
          << "\n";
      }
#endif

      return np;

    }

    AMREX_GPU_HOST_DEVICE
    amrex::Real hcp_rad (amrex::Real rad_in,
                         amrex::Real xlen,
                         amrex::Real ylen,
                         amrex::Real zlen,
                         long target) const
    {
      amrex::Real high = 1.0;
      amrex::Real low  = 1.0;

      amrex::Real radius(rad_in);

      const int max_iter = 15;
      int iter(0);

#ifndef AMREX_USE_GPU
      if(m_verbose >= 100) {
        amrex::PrintToFile("hcp_rad", amrex::Print::AllProcs)
          << "\nTarget particle count: " << target << "\n";
      }
#endif

      long pc0 = hcp_count(rad_in, xlen, ylen, zlen);
      while ((pc0 > target) && (iter < max_iter)) {
        low = high;
        high *= 2.;
        radius = high*rad_in;
        pc0 = hcp_count(radius, xlen, ylen, zlen);

#ifndef AMREX_USE_GPU
        if(m_verbose >= 100) {
          amrex::PrintToFile("hcp_rad", amrex::Print::AllProcs)
            << "iter " << iter << "   low: " << low
            << "   high: " << high << "   counts: " << pc0 << "\n";
        }
#endif
        iter++;
      }

      amrex::Real scale=0.5*(low + high);

      iter = 0;
      int keep_searching(1);
      while (keep_searching) {

        radius = scale*rad_in;
        long pc = hcp_count(radius, xlen, ylen, zlen);

#ifndef AMREX_USE_GPU
        if(m_verbose >= 100) {
          amrex::PrintToFile("hcp_rad", amrex::Print::AllProcs)
            << "Iter " << iter << "   Current count: " << pc
            << "   previous count: " << pc0 << "  scale: " << scale << "\n";
        }
#endif

        if(pc == target || (pc == pc0 && iter > 0)) {
          scale = 1.0;
#ifndef AMREX_USE_GPU
          if(m_verbose >= 100) {
            amrex::PrintToFile("hcp_rad", amrex::Print::AllProcs) << "          "
              << "We got exactly want we wanted!" << pc << "\n\n";
          }
#endif

        } else if (pc < target) {
          high  = scale;
          scale = 0.5*(low + scale);
#ifndef AMREX_USE_GPU
          if(m_verbose >= 100) {
            amrex::PrintToFile("hcp_rad", amrex::Print::AllProcs) << "          "
              << "We need more particles!  scale: " << scale << "\n\n";
          }
#endif

        } else {
          high *= ( (scale >= high) ? 2. : 1.);
          low   = amrex::max(scale, 1.0);
          scale = 0.5*(high + scale);

#ifndef AMREX_USE_GPU
          if(m_verbose >= 100) {
            amrex::PrintToFile("hcp_rad", amrex::Print::AllProcs) << "         "
              << "We need fewer particles!  scale: " << scale << "\n\n";
          }
#endif
        }

        pc0 = pc;
        iter++;
        keep_searching = ((scale != 1.0) && (iter < max_iter)) ? 1 : 0;
      }

#ifndef AMREX_USE_GPU
      if(m_verbose >= 100) {
        amrex::PrintToFile("hcp_rad", amrex::Print::AllProcs)
          << "Done: Iterations: " << iter << "  radius: " << radius << "\n\n";
      }
#endif

      return radius;
    }
};


//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: n_cube_per_fill                                         !
//                                                                      !
//  Purpose: Generate initial solids packing based on putting one       !
//           per cell                                                   !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
class nCubePer_Fill
{
  public:
    // Constructor
    nCubePer_Fill (const int base,
                   const amrex::RealVect& plo,
                   const amrex::RealVect& dx)
      : m_base(base)
      , m_plo(plo)
      , m_dx(dx)
      , m_np(0)
    {}

    void setup (MFIXInitialConditions& initial_conditions,
                const amrex::IntVect& lo,
                const amrex::IntVect& hi,
                const int icv,
                const int phase)
    {
      // This assertion should not be needed.
      AMREX_ALWAYS_ASSERT(amrex::Math::abs(m_dx[0]-m_dx[1]) < 1.e-15);
      AMREX_ALWAYS_ASSERT(amrex::Math::abs(m_dx[0]-m_dx[2]) < 1.e-15);

      int i_w, i_e, j_s, j_n, k_b, k_t;
      calc_cell_ic(m_dx[0], m_dx[1], m_dx[2],
                   initial_conditions.ic(icv).region->lo(),
                   initial_conditions.ic(icv).region->hi(),
                   m_plo.dataPtr(),
                   i_w, i_e, j_s, j_n, k_b, k_t);

      const amrex::IntVect ic_lo(AMREX_D_DECL(i_w, j_s, k_b));
      const amrex::IntVect ic_hi(AMREX_D_DECL(i_e, j_n, k_t));

      const amrex::Box mfi_bx(lo, hi);
      const amrex::Box ic_bx(ic_lo, ic_hi);

      if(!ic_bx.intersects(mfi_bx)){
        amrex::Print() << "Nothing to do. No intersection!\n";
        return;
      }

      // Intersection of mfi box and ic box
      const amrex::Box bx_int = mfi_bx&ic_bx;

      const auto ic_solid = initial_conditions.ic(icv).get_solid(phase);

      if (ic_solid == nullptr) {

        amrex::Abort("Wrong initial conditions setup");

      } else {

        const auto& diameter = ic_solid->diameter;

        const amrex::Real mean = diameter.get_mean();
        amrex::Real max_dp(mean);

        // Spacing is based on maximum particle size
        if(diameter.get_max() > 0.0)
          max_dp = diameter.get_max();

        // Scaled m_dx[0] for particle spacing
        m_dxn = m_dx[0]/static_cast<amrex::Real>(m_base);

        if (max_dp > m_dxn) {

          std::vector<std::string> regions;
          amrex::ParmParse pp("ic");
          pp.queryarr("regions", regions);
          amrex::Print() << "\n\n";
          amrex::Print() << "**************************************************************\n";
          amrex::Print() << "  Invalid particle initialization. The spacing for " << m_base << "-cube\n";
          amrex::Print() << "  packing is less than the maximum particle diameter. \n";
          amrex::Print() << "  Initial Condition Name: " << regions[icv] << "\n";
          amrex::Print() << "  Fix the inputs file.\n";
          amrex::Print() << "**************************************************************\n";
          amrex::Print() << "\n\n";
          amrex::Abort("Fix the inputs file.");
        }

        // indices
        amrex::IntVect seed_hi;

        m_seed_lo[0] = bx_int.smallEnd(0)*m_base;
        m_seed_lo[1] = bx_int.smallEnd(1)*m_base;
        m_seed_lo[2] = bx_int.smallEnd(2)*m_base;

        seed_hi[0] = (bx_int.bigEnd(0) + 1)*m_base - 1;
        seed_hi[1] = (bx_int.bigEnd(1) + 1)*m_base - 1;
        seed_hi[2] = (bx_int.bigEnd(2) + 1)*m_base - 1;

        const amrex::Box bx(m_seed_lo, seed_hi);

        m_delta_bx[0] = amrex::max(0, seed_hi[0] - m_seed_lo[0] + 1);
        m_delta_bx[1] = amrex::max(0, seed_hi[1] - m_seed_lo[1] + 1);
        m_delta_bx[2] = amrex::max(0, seed_hi[2] - m_seed_lo[2] + 1);

        m_np = m_delta_bx[0] * m_delta_bx[1] * m_delta_bx[2];
      }
    }

    int get_particles_number () const
    { return m_np; }

    template <amrex::RunOn >
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::RealVect operator() (const int p_id,
                                amrex::RandomEngine const&) const
    {
      const int local_j = int(amrex::Math::floor(amrex::Real(p_id) / amrex::Real(m_delta_bx[0]*m_delta_bx[2])));
      const int local_k = int(amrex::Math::floor(amrex::Real(p_id - local_j*m_delta_bx[0]*m_delta_bx[2]) / amrex::Real(m_delta_bx[0])));
      const int local_i = p_id - local_j*m_delta_bx[0]*m_delta_bx[2] - local_k*m_delta_bx[0];

      const int i = local_i + m_seed_lo[0];
      const int j = local_j + m_seed_lo[1];
      const int k = local_k + m_seed_lo[2];

      const amrex::Real x = m_plo[0] + (i + 0.5) * m_dxn;
      const amrex::Real y = m_plo[1] + (j + 0.5) * m_dxn;
      const amrex::Real z = m_plo[2] + (k + 0.5) * m_dxn;

      return amrex::RealVect(x, y, z);
    }

  private:
    const int m_base;

    const amrex::RealVect m_plo;

    const amrex::RealVect m_dx;

    int m_np;

    amrex::IntVect m_delta_bx;

    amrex::IntVect m_seed_lo;

    amrex::Real m_dxn;
};


//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: random fill                                             !
//                                                                      !
//  Purpose: Generate initial solids packing based on randomly placing  !
//           particles in the ic region.                                !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
class RandomFill_DEM
{
  public:
    // Constructor
    RandomFill_DEM (const amrex::RealVect& plo,
                    const amrex::RealVect& dx)
      : m_plo(plo)
      , m_dx(dx)
      , m_np(0)
    {}

    void setup (MFIXInitialConditions& initial_conditions,
                MFIXDEM& dem,
                const amrex::IntVect& lo,
                const amrex::IntVect& hi,
                const int icv,
                const int phase,
                amrex::Gpu::HostVector<amrex::Real>& h_data,
                amrex::Gpu::DeviceVector<amrex::Real>& d_data,
                const bool /*fix_seed*/)
    {
      // indices
      int i_w, i_e, j_s, j_n, k_b, k_t;

      const int maxfails = 1000;

      amrex::RealVect ic_dlo, ic_dhi, ic_len;
      amrex::Real max_dp, max_rp;

      calc_cell_ic(m_dx[0], m_dx[1], m_dx[2],
                   initial_conditions.ic(icv).region->lo(),
                   initial_conditions.ic(icv).region->hi(),
                   m_plo.dataPtr(),
                   i_w, i_e, j_s, j_n, k_b, k_t);

      // Start/end of IC domain bounds
      ic_dlo[0] = (amrex::max(lo[0], i_w)) * m_dx[0];
      ic_dlo[1] = (amrex::max(lo[1], j_s)) * m_dx[1];
      ic_dlo[2] = (amrex::max(lo[2], k_b)) * m_dx[2];

      ic_dhi[0] = (amrex::min(hi[0], i_e)+1) * m_dx[0];
      ic_dhi[1] = (amrex::min(hi[1], j_n)+1) * m_dx[1];
      ic_dhi[2] = (amrex::min(hi[2], k_t)+1) * m_dx[2];

      // physical volume of IC region intersecting this grid
      const amrex::Real ic_vol = (ic_dhi[0] - ic_dlo[0]) *
                                 (ic_dhi[1] - ic_dlo[1]) *
                                 (ic_dhi[2] - ic_dlo[2]);

      const auto ic_solid = initial_conditions.ic(icv).get_solid(phase);

      if (ic_solid == nullptr) {

        amrex::Abort("Wrong initial conditions setup");

      } else {

        const auto& diameter = ic_solid->diameter;

        amrex::Real mean = diameter.get_mean();

        // Spacing is based on maximum particle size
        if(diameter.get_max() > 0.)
          max_dp = diameter.get_max();
        else
          max_dp = diameter.get_mean();

        // If coarse-grain DEM is activated
        if (dem.cg_dem())
        {
          amrex::Real statwt = ic_solid->statwt;
          mean = std::cbrt(statwt) * mean;
          max_dp = std::cbrt(statwt) * max_dp;
        }

        // Particle count is based on mean particle size
        const int seed =
          static_cast<int>(ic_vol * ic_solid->volfrac / ((M_PI/6.0)*mean*mean*mean));

        max_rp = 0.5 * max_dp;

        ic_len[0] = ic_dhi[0] - ic_dlo[0] - max_dp;
        ic_len[1] = ic_dhi[1] - ic_dlo[1] - max_dp;
        ic_len[2] = ic_dhi[2] - ic_dlo[2] - max_dp;

        ic_dlo[0] += max_rp;
        ic_dlo[1] += max_rp;
        ic_dlo[2] += max_rp;

        const amrex::Real mindist = (1.01*max_dp)*(1.01*max_dp);

        const amrex::RealVect Oodx = {1/m_dx[0], 1/m_dx[1], 1/m_dx[2]};

        int nb(8);

        amrex::IntVect delta_bx;
        delta_bx[0] = hi[0] - lo[0] + 1;
        delta_bx[1] = hi[1] - lo[1] + 1;
        delta_bx[2] = hi[2] - lo[2] + 1;

        amrex::Vector<int> pinc(delta_bx[0]*delta_bx[1]*delta_bx[2], 0);
        amrex::Vector<int> pbin(delta_bx[0]*delta_bx[1]*delta_bx[2]*nb, 0);

        amrex::IntVect max_fitting_parts;
        max_fitting_parts[0] = static_cast<int>(amrex::Math::ceil(ic_len[0]/max_dp));
        max_fitting_parts[1] = static_cast<int>(amrex::Math::ceil(ic_len[1]/max_dp));
        max_fitting_parts[2] = static_cast<int>(amrex::Math::ceil(ic_len[2]/max_dp));

        const int max_np = max_fitting_parts[0]*max_fitting_parts[1]*max_fitting_parts[2];

        amrex::ResetRandomSeed(amrex::ParallelDescriptor::MyProc()+1);

        h_data.reserve(max_np);

        m_np = 0;
        int iterations(0);

        while(iterations < seed) {

          bool stop = false;
          int fails(0);

          while(fails < maxfails && (!stop)) {

            amrex::RealVect pos;
            pos[0] = ic_dlo[0] + ic_len[0]*amrex::Random();
            pos[1] = ic_dlo[1] + ic_len[1]*amrex::Random();
            pos[2] = ic_dlo[2] + ic_len[2]*amrex::Random();

            // Grid containing the new particle
            amrex::IntVect idx;
            idx[0] = static_cast<int>(amrex::Math::floor(pos[0]*Oodx[0]));
            idx[1] = static_cast<int>(amrex::Math::floor(pos[1]*Oodx[1]));
            idx[2] = static_cast<int>(amrex::Math::floor(pos[2]*Oodx[2]));

            // Local grid search for collisions.
            int overlaps = 0;

            amrex::IntVect seed_lo, seed_hi;
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
                    const int local_pc =
                      pbin[local_i + local_j*delta_bx[0] +
                           local_k*delta_bx[0]*delta_bx[1] +
                           p*delta_bx[0]*delta_bx[1]*delta_bx[2]];

                    const amrex::Real dist_x = h_data[local_pc*3 + 0] - pos[0] - m_plo[0];
                    const amrex::Real dist_y = h_data[local_pc*3 + 1] - pos[1] - m_plo[1];
                    const amrex::Real dist_z = h_data[local_pc*3 + 2] - pos[2] - m_plo[2];

                    const amrex::Real dist = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;

                    if(dist < mindist)
                      overlaps++;
                  }
                }
              }
            }

            if(overlaps == 0)
            {
              h_data.push_back(pos[0] + m_plo[0]);
              h_data.push_back(pos[1] + m_plo[1]);
              h_data.push_back(pos[2] + m_plo[2]);

              const int local_idx_x = idx[0] - lo[0];
              const int local_idx_y = idx[1] - lo[1];
              const int local_idx_z = idx[2] - lo[2];

              const int local_idx =
                local_idx_x + local_idx_y*delta_bx[0] +
                local_idx_z*delta_bx[0]*delta_bx[1];

              const int pinc_value = pinc[local_idx];
              pbin[local_idx + pinc_value*delta_bx[0]*delta_bx[1]*delta_bx[2]] = m_np;

              pinc[local_idx] = pinc_value + 1;

              if(pinc_value + 2 >= nb)
              {
                nb += 2;
                pbin.resize(delta_bx[0]*delta_bx[1]*delta_bx[2]*nb, 0);
              }

              m_np++;

              stop = true;
            }
            else
            {
              fails++;
            }
          }

          iterations++;
        }

        d_data.resize(3*m_np);
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_data.begin(), h_data.end(), d_data.begin());

        m_h_data_ptr = h_data.data();
        m_d_data_ptr = d_data.data();
      }
    }

    int get_particles_number () const
    { return m_np; }

    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::RealVect operator() (const int p_id,
                                amrex::RandomEngine const&) const
    {
      amrex::Real x;
      amrex::Real y;
      amrex::Real z;

      if (run_on == amrex::RunOn::Gpu) {
        x = m_d_data_ptr[p_id*3+0];
        y = m_d_data_ptr[p_id*3+1];
        z = m_d_data_ptr[p_id*3+2];
      } else if (run_on == amrex::RunOn::Cpu) {
        x = m_h_data_ptr[p_id*3+0];
        y = m_h_data_ptr[p_id*3+1];
        z = m_h_data_ptr[p_id*3+2];
      }

      return amrex::RealVect(x, y, z);
    }

  private:
    const amrex::RealVect m_plo;

    const amrex::RealVect m_dx;

    int m_np;

    amrex::Real* m_h_data_ptr;
    amrex::Real* m_d_data_ptr;
};


//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: random fill                                             !
//                                                                      !
//  Purpose: Generate initial solids packing based on randomly placing  !
//           particles in the ic region.                                !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
class RandomFill_PIC
{
  public:
    // Constructor
    RandomFill_PIC (const amrex::RealVect& plo,
                    const amrex::RealVect& dx)
      : m_plo(plo)
      , m_dx(dx)
      , m_np(0)
    {}

    void setup (MFIXInitialConditions& initial_conditions,
                const amrex::IntVect& lo,
                const amrex::IntVect& hi,
                const int icv,
                const int phase,
                const bool /*fix_seed*/)
    {
      const auto& ics = initial_conditions.ic();
      const auto& ic = ics[icv];
      const auto& ic_region = ic.region;

      const amrex::Real* ic_rlo = ic_region->lo();
      const amrex::Real* ic_rhi = ic_region->hi();

      amrex::IntVect seed_hi;

      m_seed_lo[0] = amrex::max<int>(lo[0], static_cast<int>(amrex::Math::ceil((ic_rlo[0]-m_plo[0])/m_dx[0] - .5)));
      m_seed_lo[1] = amrex::max<int>(lo[1], static_cast<int>(amrex::Math::ceil((ic_rlo[1]-m_plo[1])/m_dx[1] - .5)));
      m_seed_lo[2] = amrex::max<int>(lo[2], static_cast<int>(amrex::Math::ceil((ic_rlo[2]-m_plo[2])/m_dx[2] - .5)));

      seed_hi[0] = amrex::min<int>(hi[0], static_cast<int>(amrex::Math::floor((ic_rhi[0]-m_plo[0])/m_dx[0] - .5)));
      seed_hi[1] = amrex::min<int>(hi[1], static_cast<int>(amrex::Math::floor((ic_rhi[1]-m_plo[1])/m_dx[1] - .5)));
      seed_hi[2] = amrex::min<int>(hi[2], static_cast<int>(amrex::Math::floor((ic_rhi[2]-m_plo[2])/m_dx[2] - .5)));

      const amrex::Box bx(m_seed_lo, seed_hi);

      m_delta_bx[0] = amrex::max(0, seed_hi[0] - m_seed_lo[0] + 1);
      m_delta_bx[1] = amrex::max(0, seed_hi[1] - m_seed_lo[1] + 1);
      m_delta_bx[2] = amrex::max(0, seed_hi[2] - m_seed_lo[2] + 1);

      // Number of cells in IC region
      const amrex::Real ic_cells = m_delta_bx[0] * m_delta_bx[1] * m_delta_bx[2];

      const auto ic_solid = ic.get_solid(phase);

      if (ic_solid == nullptr) {

        amrex::Abort("Wrong initial conditions setup");

      } else {

        // Particle properties: Mean diameter, stat weight, volume
        const amrex::Real mean = ic_solid->diameter.get_mean();
        const amrex::Real statwt(ic_solid->statwt);
        const amrex::Real parcel_volume = statwt*((M_PI/6.0)*mean*mean*mean);

        // Parcels per cell to satisfy the IC condition
        const amrex::Real parcels_per_cell =
            (m_dx[0] * m_dx[1] * m_dx[2] * ic_solid->volfrac / parcel_volume);

        const int whole_parcels_per_cell = static_cast<int>(amrex::Math::floor(parcels_per_cell));
        m_wppc = whole_parcels_per_cell;

        // This is the total number of particles that will be generated.
        m_np = static_cast<int>(ic_cells) * whole_parcels_per_cell;

        amrex::ResetRandomSeed(amrex::ParallelDescriptor::MyProc()+1);
      }
    }

    int get_particles_number () const
    { return m_np; }

    template <amrex::RunOn >
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::RealVect operator() (const int p_id,
                                amrex::RandomEngine const& engine) const
    {
      const int local_j = int(amrex::Math::floor(amrex::Real(p_id) / amrex::Real(m_wppc*m_delta_bx[0]*m_delta_bx[2])));

      const int local_k = int(amrex::Math::floor(amrex::Real(p_id - local_j*m_wppc*m_delta_bx[0]*m_delta_bx[2]) / amrex::Real(m_wppc*m_delta_bx[0])));

      const int local_i = int(amrex::Math::floor(amrex::Real(p_id - m_wppc*(local_k*m_delta_bx[0] + local_j*m_delta_bx[0]*m_delta_bx[2])) / amrex::Real(m_wppc)));

      const int i = local_i + m_seed_lo[0];
      const int j = local_j + m_seed_lo[1];
      const int k = local_k + m_seed_lo[2];

      const amrex::Real x = m_plo[0] + (i * m_dx[0]) + m_dx[0]*amrex::Random(engine);
      const amrex::Real y = m_plo[1] + (j * m_dx[1]) + m_dx[1]*amrex::Random(engine);
      const amrex::Real z = m_plo[2] + (k * m_dx[2]) + m_dx[2]*amrex::Random(engine);

      return amrex::RealVect(x, y, z);
    }

  private:
    const amrex::RealVect m_plo;

    const amrex::RealVect m_dx;

    int m_np;

    int m_wppc;

    amrex::IntVect m_delta_bx;

    amrex::IntVect m_seed_lo;
};


  // TODO: At this point, we have seeded the bulk of the parcels. However there
  // are two oversites --
  // 1) This will get us close to the requested solids volume fraction, but
  //    not exact. We need to deal with the remainder.
  // 2) If the IC region is dilute (with less than one parcel per cell), this
  //    will not create any parcles. This could be addressed by adding a
  //    stride to the above logic so we skip past the stride cells.

#endif
