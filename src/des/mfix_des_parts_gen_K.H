#ifndef _MFIX_DES_PARTS_GEN_K_H_
#define _MFIX_DES_PARTS_GEN_K_H_

//#include <mfix_particle_generator.H>
#include <mfix_calc_cell.H>

#include <limits>
//#include <string>
//#include <iostream>
//#include <iomanip>
//#include <fstream>
//#include <sstream>

#include <mfix_ic_parms.H>
//#include <mfix_solids_parms.H>
//#include <mfix_dem_parms.H>
//#include <mfix_pic_parms.H>

#include <AMReX_IntVect.H>
#include <AMReX_RealVect.H>
#include <AMReX_PODVector.H>

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
    Hex_ClosePack (const amrex::GpuArray<amrex::Real,3>& plo,
                   const amrex::GpuArray<amrex::Real,3>& dx)
      : m_plo(plo)
      , m_dx(dx)
    {}

    void setup (const amrex::IntVect& lo,
                const amrex::IntVect& hi,
                const int icv,
                const int type)
    {
      // physical volume of IC region
      const amrex::Real ic_vol = IC::ic[icv].region->volume();

      const amrex::Real mean = IC::ic[icv].solids[type].diameter.mean;

      // Total particle count is based on mean particle size
      const long seed =
        static_cast<long>(ic_vol * IC::ic[icv].solids[type].volfrac / ((M_PI/6.0)*mean*mean*mean));

      const amrex::Real xlen = IC::ic[icv].region->hi(0) - IC::ic[icv].region->lo(0);
      const amrex::Real ylen = IC::ic[icv].region->hi(1) - IC::ic[icv].region->lo(1);
      const amrex::Real zlen = IC::ic[icv].region->hi(2) - IC::ic[icv].region->lo(2);

      // This is the radius we need for spacing particles to get the
      // desired volume fraction and particle count.
      amrex::Real eff_rad = hcp_rad(0.5*mean, xlen, ylen, zlen, seed);
      //printf ("\nPacking radius: %9.6f\n", eff_rad);

      amrex::RealVect rbx_lo(m_plo[0] + m_dx[0]*static_cast<amrex::Real>(lo[0]),
                             m_plo[1] + m_dx[1]*static_cast<amrex::Real>(lo[1]),
                             m_plo[2] + m_dx[2]*static_cast<amrex::Real>(lo[2]));

      amrex::RealVect rbx_hi(m_plo[0] + m_dx[0]*static_cast<amrex::Real>(hi[0]+1),
                             m_plo[1] + m_dx[1]*static_cast<amrex::Real>(hi[1]+1),
                             m_plo[2] + m_dx[2]*static_cast<amrex::Real>(hi[2]+1));

      m_ic_dlo = amrex::RealVect(IC::ic[icv].region->lo(0),
                                 IC::ic[icv].region->lo(1),
                                 IC::ic[icv].region->lo(2));

      amrex::RealVect ic_dhi(IC::ic[icv].region->hi(0),
                             IC::ic[icv].region->hi(1),
                             IC::ic[icv].region->hi(2));

      // Verify that the IC region overlaps with the IC region.
      if (ic_dhi[0] < rbx_lo[0] || rbx_hi[0] < m_ic_dlo[0] ||
          ic_dhi[1] < rbx_lo[1] || rbx_hi[1] < m_ic_dlo[1] ||
          ic_dhi[2] < rbx_lo[2] || rbx_hi[2] < m_ic_dlo[2]) {

        //amrex::Print() << "Nothing to do in this MFIter -- bail out\n";
        return;
      }


      const int kmax = hcp_k(eff_rad, zlen);

      int klo(0);

      for(;klo<kmax;klo++) {
        if(m_ic_dlo[2] + (1+klo*sqrt6x2o3)*eff_rad > rbx_lo[2])
          break;
      }

      int khi(klo);
      for(;khi<kmax;khi++) {
        if(m_ic_dlo[2] + (1+khi*sqrt6x2o3)*eff_rad >= rbx_hi[2])
          break;
      }
      khi--;

      const int jmax = std::max(hcp_j(eff_rad, ylen, 0),
                                hcp_j(eff_rad, ylen, 1));

      int jlo(0);
      for(; jlo<jmax; jlo++) {
        if(m_ic_dlo[1] + (1+jlo*sqrt3)*eff_rad > rbx_lo[1])
          break;
      }

      int jhi(jlo);
      const amrex::Real third(1.0/3.0);
      for(; jhi<jmax; jhi++) {
        if(m_ic_dlo[1] + (1+(jhi + third)*sqrt3)*eff_rad >= rbx_hi[1])
          break;
      }
      jhi--;


      const int imax = std::max(hcp_i(eff_rad, xlen, 0),
                                hcp_i(eff_rad, xlen, 1));

      int ilo(0);
      for(; ilo<imax; ilo++) {
        if(m_ic_dlo[0] + (2*ilo+1)*eff_rad > rbx_lo[0])
          break;
      }

      int ihi(ilo);
      for(; ihi<imax; ihi++) {
        if(m_ic_dlo[0] + (2*ihi+1)*eff_rad >= rbx_hi[0])
          break;
      }
      ihi--;

      m_seed_lo = amrex::IntVect(ilo, jlo, klo);
      amrex::IntVect seed_hi(ihi, jhi, khi);

      m_delta_bx = amrex::IntVect(amrex::max(0, seed_hi[0] - m_seed_lo[0] + 1),
                                  amrex::max(0, seed_hi[1] - m_seed_lo[1] + 1),
                                  amrex::max(0, seed_hi[2] - m_seed_lo[2] + 1));

      const amrex::Box bx(m_seed_lo, seed_hi);
      m_max_rp = eff_rad;

      m_np = m_delta_bx[0] * m_delta_bx[1] * m_delta_bx[2];
    }

    int get_particles_number () const
    { return m_np; }

    template <amrex::RunOn >
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::RealVect operator() (const int p_id,
                                amrex::RandomEngine const&) const
    {
      const int local_i = int(amrex::Math::floor(amrex::Real(p_id) / amrex::Real(m_delta_bx[1]*m_delta_bx[2])));
      const int local_k = int(amrex::Math::floor(amrex::Real(p_id - local_i*m_delta_bx[1]*m_delta_bx[2]) / amrex::Real(m_delta_bx[1])));
      const int local_j = p_id - local_i*m_delta_bx[1]*m_delta_bx[2] - local_k*m_delta_bx[1];

      const int i = local_i + m_seed_lo[0];
      const int j = local_j + m_seed_lo[1];
      const int k = local_k + m_seed_lo[2];

      const amrex::Real x = m_ic_dlo[0] + m_max_rp * (1.0 + static_cast<amrex::Real>(2*i + (j+k)%2));
      const amrex::Real y = m_ic_dlo[1] + m_max_rp * (1.0 + sqrt3o3*(static_cast<amrex::Real>(3*j + k%2)));
      const amrex::Real z = m_ic_dlo[2] + m_max_rp * (1.0 + sqrt6x2o3*(static_cast<amrex::Real>(k)));

      return amrex::RealVect(x, y, z);
    }

  private:
    const amrex::Real sqrt6x2o3 = std::sqrt(6.0)*(2.0/3.0);
    const amrex::Real sqrt3o3 = std::sqrt(3.0)/3.0;
    const amrex::Real sqrt3 = std::sqrt(3.0);

    const amrex::GpuArray<amrex::Real,3>& m_plo;

    const amrex::GpuArray<amrex::Real,3>& m_dx;

    int m_np;

    amrex::IntVect m_delta_bx;

    amrex::IntVect m_seed_lo;

    amrex::RealVect m_ic_dlo;

    amrex::Real m_max_rp;

    // Auxiliary functions
    int hcp_i (amrex::Real rad,
               amrex::Real xlen,
               int offset_int) const 
    {
      int i(0);
      const amrex::Real offset = static_cast<amrex::Real>(offset_int);
      while ((1.0 + 2.0*static_cast<amrex::Real>(i) + offset)*rad <= xlen-rad) {
        i++;
      }
      return i;
    }

    int hcp_j (amrex::Real rad,
               amrex::Real ylen,
               int offset_int) const
    {
      int j(0);
      const amrex::Real third = amrex::Real(3.0);
      const amrex::Real offset = static_cast<amrex::Real>(offset_int);
      while ((1.0 + sqrt3*static_cast<amrex::Real>(j) + offset*third)*rad <=ylen-rad) {
        j++;
      }
      return j;
    }

    int hcp_k (amrex::Real rad,
               amrex::Real zlen) const
    {
      int k(0);
      while ((1.0 + sqrt6x2o3*static_cast<amrex::Real>(k))*rad <= zlen-rad) {
        k++;
      }
      return k;
    }

    long hcp_count (amrex::Real rad,
                    amrex::Real xlen,
                    amrex::Real ylen,
                    amrex::Real zlen) const
    {
      const int a1 = hcp_i(rad, xlen, 0);
      const int a2 = hcp_i(rad, xlen, 1);

      const int b1 = hcp_j(rad, ylen, 0);
      const int b2 = hcp_j(rad, ylen, 1);

      const int b1r = b1%2;
      const int b2r = b2%2;

      const int c1 = hcp_k(rad, zlen);
      const int c1r = c1%2;

      return ((a1*((b1 - b1r)/2 + b1r)) + (a2*((b1 - b1r)/2)))*( (c1 - c1r)/2 + c1r) +
             ((a2*((b2 - b2r)/2 + b1r)) + (a1*((b2 - b2r)/2)))*( (c1 - c1r)/2      );
    }

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

      //amrex::Print() << "\nTarget particle count: " << target << "\n";

      long pc0 = hcp_count(rad_in, xlen, ylen, zlen);
      while ((pc0 > target) && (iter < max_iter)) {
        low = high;
        high *= 2.;
        radius = high*rad_in;
        pc0 = hcp_count(radius, xlen, ylen, zlen);

        //amrex::Print() << "iter " << iter << "   low: " << low
        //    << "   high: " << high << "   counts: " << pc0 << "\n";
        iter++;
      }

      amrex::Real scale=0.5*(low + high);

      iter = 0;
      int keep_searching(1);
      while (keep_searching) {

        radius = scale*rad_in;
        long pc = hcp_count(radius, xlen, ylen, zlen);

        //amrex::Print() << "iter " << iter << "   current: " << pc
        //    << "   previous: " << pc0 << "  scale: " << scale << "\n";

        if(pc == target || (pc == pc0 && iter > 0)) {
          //amrex::Print() << "We got exactly want we wanted! " << pc << "\n";
          scale = 1.0;

        } else if (pc < target) {
          //amrex::Print() << "We need more particles!  scale: " << scale << "\n";

          high  = scale;
          scale = 0.5*(low + scale);
          //scale = amrex::max(0.5*(low + scale),1.0);

        } else {
          //amrex::Print() << "We need fewer particles!  scale: " << scale << "\n";
          if( scale >= high)
            high *= 2.;

          low   = amrex::max(scale, 1.0);
          scale = 0.5*(high + scale);
        }

        pc0 = pc;
        iter++;
        keep_searching = ((scale != 1.0) && (iter < max_iter)) ? 1 : 0;
      }

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
                   const amrex::GpuArray<amrex::Real,3>& plo,
                   const amrex::GpuArray<amrex::Real,3>& dx)
      : m_base(base)
      , m_plo(plo)
      , m_dx(dx)
    {}
    
    void setup (const amrex::IntVect& lo,
                const amrex::IntVect& hi,
                const int icv,
                const int type)
    {
      // This assertion should not be needed.
      AMREX_ALWAYS_ASSERT(amrex::Math::abs(m_dx[0]-m_dx[1]) < 1.e-15);
      AMREX_ALWAYS_ASSERT(amrex::Math::abs(m_dx[0]-m_dx[2]) < 1.e-15);

      int i_w, i_e, j_s, j_n, k_b, k_t;
      calc_cell_ic(m_dx[0], m_dx[1], m_dx[2],
                   IC::ic[icv].region->lo(),
                   IC::ic[icv].region->hi(),
                   m_plo.data(),
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

      const amrex::Real mean = IC::ic[icv].solids[type].diameter.mean;
      amrex::Real max_dp(mean);

      // Spacing is based on maximum particle size
      if(IC::ic[icv].solids[type].diameter.max > 0.0)
        max_dp = IC::ic[icv].solids[type].diameter.max;

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

    int get_particles_number () const
    { return m_np; }

    template <amrex::RunOn >
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
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

    const amrex::GpuArray<amrex::Real,3>& m_plo;

    const amrex::GpuArray<amrex::Real,3>& m_dx;

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
    RandomFill_DEM (const amrex::GpuArray<amrex::Real,3>& plo,
                    const amrex::GpuArray<amrex::Real,3>& dx)
      : m_plo(plo)
      , m_dx(dx)
    {}

    void setup (const amrex::IntVect& lo,
                const amrex::IntVect& hi,
                const int icv,
                const int type,
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
                   IC::ic[icv].region->lo(),
                   IC::ic[icv].region->hi(),
                   m_plo.data(),
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

      amrex::Real mean = IC::ic[icv].solids[type].diameter.mean;

      // Spacing is based on maximum particle size
      if(IC::ic[icv].solids[type].diameter.max > 0.)
        max_dp = IC::ic[icv].solids[type].diameter.max;
      else
        max_dp = IC::ic[icv].solids[type].diameter.mean;

      // If coarse-grain DEM is activated
      if (DEM::cg_dem)
      {
        amrex::Real statwt = IC::ic[icv].solids[type].statwt;
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
                  const amrex::Real* plo_ptr = m_plo.data();
    
                  const int local_pc =
                    pbin[local_i + local_j*delta_bx[0] +
                         local_k*delta_bx[0]*delta_bx[1] +
                         p*delta_bx[0]*delta_bx[1]*delta_bx[2]];
    
                  const amrex::Real dist_x = h_data[local_pc*3 + 0] - pos[0] - plo_ptr[0];
                  const amrex::Real dist_y = h_data[local_pc*3 + 1] - pos[1] - plo_ptr[1];
                  const amrex::Real dist_z = h_data[local_pc*3 + 2] - pos[2] - plo_ptr[2];
    
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

    int get_particles_number () const
    { return m_np; }

    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
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
    const amrex::GpuArray<amrex::Real,3>& m_plo;

    const amrex::GpuArray<amrex::Real,3>& m_dx;

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
    RandomFill_PIC (const amrex::GpuArray<amrex::Real,3>& plo,
                    const amrex::GpuArray<amrex::Real,3>& dx)
      : m_plo(plo)
      , m_dx(dx)
    {}

    void setup (const amrex::IntVect& lo,
                const amrex::IntVect& hi,
                const int icv,
                const int type,
                const bool /*fix_seed*/)
    {
      const amrex::Real* ic_rlo = IC::ic[icv].region->lo();
      const amrex::Real* ic_rhi = IC::ic[icv].region->hi();

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

      // Particle properties: Mean diameter, stat weight, volume
      const amrex::Real mean = IC::ic[icv].solids[type].diameter.mean;
      const amrex::Real statwt(IC::ic[icv].solids[type].statwt);
      const amrex::Real parcel_volume = statwt*((M_PI/6.0)*mean*mean*mean);

      // Parcels per cell to satisfy the IC condition
      const amrex::Real parcels_per_cell =
          (m_dx[0] * m_dx[1] * m_dx[2] * IC::ic[icv].solids[type].volfrac / parcel_volume);

      const int whole_parcels_per_cell = static_cast<int>(amrex::Math::floor(parcels_per_cell));
      m_wppc = whole_parcels_per_cell;

      // This is the total number of particles that will be generated.
      m_np = static_cast<int>(ic_cells) * whole_parcels_per_cell;

      amrex::ResetRandomSeed(amrex::ParallelDescriptor::MyProc()+1);
    }

    int get_particles_number () const
    { return m_np; }

    template <amrex::RunOn >
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
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
    const amrex::GpuArray<amrex::Real,3>& m_plo;

    const amrex::GpuArray<amrex::Real,3>& m_dx;

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
