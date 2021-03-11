#include <mfix_init_fluid.H>

#include <mfix_calc_fluid_coeffs.H>
#include <mfix_calc_cell.H>
#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>
#include <mfix_ic_parms.H>
#include <Redistribution.H>

using namespace amrex;

// Forward declarations
void set_ic_vel (const Box& sbx, const Box& domain,
                 const Real dx, const Real dy, const Real dz,
                 const GpuArray<Real, 3>& plo, FArrayBox& vel_g_fab);

void set_ic_temp (const Box& sbx, const Box& domain,
                  const Real dx, const Real dy, const Real dz,
                  const GpuArray<Real, 3>& plo, FArrayBox& T_g_fab);

void set_ic_species_g (const Box& sbx, const Box& domain,
                       const Real dx, const Real dy, const Real dz,
                       const GpuArray<Real, 3>& plo, FArrayBox& X_gk_fab);

void init_helix (const Box& bx, const Box& domain, FArrayBox& vel_g_fab,
                 const Real dx, const Real dy, const Real dz);

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: init_fluid                                              !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

void init_fluid (const Box& sbx,
                 const Box& bx,
                 const Box& domain,
                 const MFIter& mfi,
                 LevelData& ld,
                 const Real dx,
                 const Real dy,
                 const Real dz,
                 const Real xlength,
                 const Real ylength,
                 const Real zlength,
                 const GpuArray<Real, 3>& plo,
                 bool test_tracer_conservation,
                 const int advect_enthalpy,
                 const int advect_fluid_species)
{
  // Set user specified initial conditions (IC)
  // Set initial fluid velocity
  set_ic_vel(sbx, domain, dx, dy, dz, plo, (*ld.vel_g)[mfi]);

  // init_periodic_vortices (bx, domain, vel_g_fab, dx, dy, dz);
  // init_helix (bx, domain, vel_g_fab, dx, dy, dz);

  // Set the initial fluid density
  Array4<Real> const& ro_g = ld.ro_g->array(mfi);

  const Real ro_g0  = FLUID::ro_g0;

  ParallelFor(sbx, [ro_g,ro_g0] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  { ro_g(i,j,k) = ro_g0; });

  // Set the initial fluid tracer
  Array4<Real> const& trac = ld.trac->array(mfi);

  const Real trac_0 = FLUID::trac_0;

  ParallelFor(sbx, [trac,trac_0] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  { trac(i,j,k) = trac_0; });

  // Set the initial fluid temperature
  if (advect_enthalpy) {
    set_ic_temp(sbx, domain, dx, dy, dz, plo, (*ld.T_g)[mfi]);
  }

  // Fluid SPECIES Initialization
  if (advect_fluid_species) {
    // Set the initial fluid species mass fractions
    set_ic_species_g(sbx, domain, dx, dy, dz, plo, (*ld.X_gk)[mfi]);
  }

  // Initialize all the fluid and fluid species parameters
  init_fluid_parameters(bx, mfi, ld, advect_enthalpy, advect_fluid_species);

  if (test_tracer_conservation)
    init_periodic_tracer(bx, domain, (*ld.vel_g)[mfi], (*ld.trac)[mfi], dx, dy, dz);
}

void init_helix (const Box& bx,
                 const Box& domain,
                 FArrayBox& vel_g_fab,
                 const Real dx,
                 const Real dy,
                 const Real dz)
{
  Array4<Real> const& velocity = vel_g_fab.array();

  const int plane = 3;

  const amrex::Real fac = 0.01;

  switch (plane)
  {
    case 1:  // around x-axis
      amrex::ParallelFor(bx,[fac,dy,dz,velocity]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real y = (Real(j) + .5) * dy - .0016;
        Real z = (Real(k) + .5) * dz - .0016;
        Real r = std::sqrt(y*y + z*z);

        velocity(i,j,k,0) = 0.0;
        velocity(i,j,k,1) =  fac * (z/r);
        velocity(i,j,k,2) = -fac * (y/r);
      });
      break;

    case 2:  // around y-axis
      amrex::ParallelFor(bx,[fac,dx,dz,velocity]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real x = (Real(i) + .5) * dx - .0016;
        Real z = (Real(k) + .5) * dz - .0016;
        Real r = std::sqrt(x*x + z*z);

        velocity(i,j,k,0) =  fac * (z/r);
        velocity(i,j,k,1) = 0.0;
        velocity(i,j,k,2) = -fac * (x/r);
      });
      break;

    case 3:  // around z-axis
      amrex::ParallelFor(bx,[fac,dx,dy,velocity]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real x = (Real(i) + .5) * dx;
        Real y = (Real(j) + .5) * dy;
        Real r = std::sqrt(x*x + y*y);

        velocity(i,j,k,0) =  fac * (y/r);
        velocity(i,j,k,1) = -fac * (x/r);
        velocity(i,j,k,2) = 0.0;
      });
      break;

    default:
      amrex::Abort("Error: wrong plane number");
      break;
  }
}

void init_periodic_vortices (const Box& bx,
                             const Box& domain,
                             FArrayBox& vel_g_fab,
                             const Real dx,
                             const Real dy,
                             const Real dz)
{
  const amrex::Real twopi = 2. * M_PI;

  Array4<Real> const& velocity = vel_g_fab.array();

  const int plane = 1;

  switch (plane)
  {
    case 1:  // x-y plane
      // x-direction
      amrex::ParallelFor(bx,[twopi,dx,dy,velocity]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real x = (Real(i) + .5) * dx;
        Real y = (Real(j) + .5) * dy;

        velocity(i,j,k,0) = std::tanh(30*(.25 - amrex::Math::abs(y-.5)));
        velocity(i,j,k,1) = .05 * std::sin(twopi*x);
        velocity(i,j,k,2) = 0.;
      });
      break;

    case 2:  // x-z plane
      // x-direction
      amrex::ParallelFor(bx,[twopi,dx,dz,velocity]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real x = (Real(i) + .5) * dx;
        Real z = (Real(k) + .5) * dz;

        velocity(i,j,k,0) = std::tanh(30*(.25 - amrex::Math::abs(z-.5)));
        velocity(i,j,k,1) = 0.;
        velocity(i,j,k,2) = .05 * std::sin(twopi*x);
      });
      break;

    case 3:  // y-z plane
      // x-direction
      amrex::ParallelFor(bx,[twopi,dy,dz,velocity]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real y = (Real(j) + .5) * dy;
        Real z = (Real(k) + .5) * dz;

        velocity(i,j,k,0) = 0.;
        velocity(i,j,k,1) = std::tanh(30*(.25 - amrex::Math::abs(z-.5)));
        velocity(i,j,k,2) = .05 * std::sin(twopi*y);
      }); break;

    default:
      amrex::Abort("Error: wrong plane number");
      break;
  }
}


void init_periodic_tracer (const Box& bx,
                           const Box& domain,
                           FArrayBox& vel_g_fab,
                           FArrayBox& trac_fab,
                           const Real dx,
                           const Real dy,
                           const Real dz)
{
    const amrex::Real twopi(2. * M_PI);

    Array4<Real> const& trac = trac_fab.array();
    Array4<Real> const&  vel = vel_g_fab.array();

    int dir(0);
    const amrex::Real A(1.0);

    amrex::Real L(0); // domain size
    amrex::Real C(0); // sin coefficient

    dir += (domain.bigEnd(1) == domain.bigEnd(2)) ? 1 : 0;
    dir += (domain.bigEnd(0) == domain.bigEnd(2)) ? 2 : 0;
    dir += (domain.bigEnd(0) == domain.bigEnd(1)) ? 4 : 0;

    switch (dir)
      {
      case 1:  // x-direction

        L = Real(domain.bigEnd(0)+1) * dx;
        C = twopi / L;
        amrex::ParallelFor(bx,[A,C,dx,dy,dz,trac,vel]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            Real x = (Real(i) + .5) * dx - .00037;
            Real y = (Real(j) + .5) * dy - .00073;
            Real z = (Real(k) + .5) * dz - .00123;

            trac(i,j,k) = A*( std::sin(C*(y+z) - 0.00042) + 1.0) * exp(x);

            vel(i,j,k,1) += 0.1*( std::sin(C*(x+z) - 0.00042) + 1.0) * exp(y);
            vel(i,j,k,2) += 0.1*( std::sin(C*(x+y) - 0.00042) + 1.0) * exp(z);
        });
        break;

    case 2: // y-direction

        L = Real(domain.bigEnd(1)+1) * dy;
        C = twopi / L;
        amrex::ParallelFor(bx,[A,C,dx,dy,dz,trac,vel]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            Real y = (Real(j) + .5) * dy - .00037;

            Real x = (Real(i) + .5) * dx - .00073;
            Real z = (Real(k) + .5) * dz - .00123;

            trac(i,j,k) = A*( std::sin(C*(x+z) - 0.00042) + 1.0) * exp(y);

            vel(i,j,k,0) += 0.1*( std::sin(C*(y+z) - 0.00042) + 1.0) * exp(x);
            vel(i,j,k,2) += 0.1*( std::sin(C*(y+x) - 0.00042) + 1.0) * exp(z);

        });
        break;

    case 4: // z-direction

        L = Real(domain.bigEnd(2)+1) * dz;
        C = twopi / L;
        amrex::ParallelFor(bx,[A,C,dx,dy,dz,trac,vel]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real z = (Real(k) + .5) * dz - .00037;

            Real x = (Real(i) + .5) * dx - .00123;
            Real y = (Real(j) + .5) * dy - .00073;

            trac(i,j,k) = A*( std::sin(C*(x+y) - 0.00042) + 1.0) * exp(z);

            vel(i,j,k,0) += 0.1*( std::sin(C*(z+y) - 0.00042) + 1.0) * exp(x);
            vel(i,j,k,1) += 0.1*( std::sin(C*(z+x) - 0.00042) + 1.0) * exp(y);
        });
        break;

    default:
        amrex::Abort("Error: wrong direction number");
        break;
    }
}

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: init_fluid_parameters                                   !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void init_fluid_parameters (const Box& bx,
                            const MFIter& mfi,
                            LevelData& ld,
                            const int advect_enthalpy,
                            const int advect_fluid_species)
{
  Array4<Real> const& MW_g = ld.MW_g->array(mfi);

  // Initialize mu_g
  calc_mu_g(bx, (*ld.mu_g)[mfi]);

  // Initialize D_gk
  if (advect_fluid_species) {
    calc_D_gk(bx, (*ld.D_gk)[mfi]);
  }

  // Initialize k_g
  if (advect_enthalpy)
    calc_k_g(bx, (*ld.k_g)[mfi], (*ld.T_g)[mfi]);

  // Initialize cp_gk and h_gk
  if (advect_enthalpy and advect_fluid_species) {
    calc_cp_gk(bx, (*ld.cp_gk)[mfi], (*ld.T_g)[mfi]);
    calc_h_gk(bx, (*ld.h_gk)[mfi], (*ld.cp_gk)[mfi], (*ld.T_g)[mfi]);
  }

  // Initialize MW_g0, cp_g, h_g, cp_gk, and h_gk
  if (not FLUID::is_a_mixture)
  {
    const Real MW_g0  = FLUID::MW_g0;

    ParallelFor(bx, [MW_g,MW_g0] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    { MW_g(i,j,k) = MW_g0; });

    // Initialize Cp_g and k_g
    if (advect_enthalpy) {
      calc_cp_g(bx, (*ld.cp_g)[mfi], (*ld.T_g)[mfi]);
      calc_h_g(bx, (*ld.h_g)[mfi], (*ld.cp_g)[mfi], (*ld.T_g)[mfi]);
    }
  }
  else // fluid_is_a_mixture == true
  {
    const int nspecies_g = FLUID::nspecies;

    Array4<const Real> const& X_gk = ld.X_gk->const_array(mfi);

    Gpu::DeviceVector< Real > MW_gk0_d(nspecies_g);
    Gpu::copyAsync(Gpu::hostToDevice, FLUID::MW_gk0.begin(), FLUID::MW_gk0.end(), MW_gk0_d.begin());

    Real* p_MW_gk0 = MW_gk0_d.data();

    // Set initial species molecular weights and fluid molecular weight
    ParallelFor(bx, [nspecies_g,X_gk,MW_g,p_MW_gk0]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      Real MW_g_sum(0);

      for (int n(0); n < nspecies_g; n++)
        MW_g_sum += X_gk(i,j,k,n) / p_MW_gk0[n];

      MW_g(i,j,k) = 1./MW_g_sum;
    });

    if (advect_enthalpy)
    {
      Array4<Real> const& cp_g  = ld.cp_g->array(mfi);
      Array4<Real> const& cp_gk = ld.cp_gk->array(mfi);
      Array4<Real> const& h_g   = ld.h_g->array(mfi);
      Array4<Real> const& h_gk  = ld.h_gk->array(mfi);

      // Update cp_g and h_g from species quantities
      ParallelFor(bx, [nspecies_g,X_gk,cp_gk,h_gk,cp_g,h_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real cp_g_sum(0);
        Real h_g_sum(0);

        for (int n(0); n < nspecies_g; n++) {
          cp_g_sum += X_gk(i,j,k,n) * cp_gk(i,j,k,n);
          h_g_sum += X_gk(i,j,k,n) * h_gk(i,j,k,n);
        }

        cp_g(i,j,k) = cp_g_sum;
        h_g(i,j,k) = h_g_sum;
      });

    }

    Gpu::synchronize();

  }
}

//!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//!                                                                      !
//!  Subroutine: SET_IC_VEL                                              !
//!                                                                      !
//!  Purpose: Set velocity initial conditions.                           !
//!                                                                      !
//!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void set_ic_vel (const Box& sbx,
                 const Box& domain,
                 const Real dx,
                 const Real dy,
                 const Real dz,
                 const GpuArray<Real, 3>& plo,
                 FArrayBox& vel_g_fab)
{
  const IntVect slo(sbx.loVect());
  const IntVect shi(sbx.hiVect());

  const IntVect domlo(domain.loVect());
  const IntVect domhi(domain.hiVect());

  Array4<Real> const& velocity = vel_g_fab.array();

  // Set the initial conditions.
  for(int icv(0); icv < IC::ic.size(); ++icv)
  {

    int i_w(0), j_s(0), k_b(0);
    int i_e(0), j_n(0), k_t(0);

    calc_cell_ic(dx, dy, dz,
                 IC::ic[icv].region->lo(),
                 IC::ic[icv].region->hi(),
                 plo.data(),
                 i_w, i_e, j_s, j_n, k_b, k_t);

    // Use the volume fraction already calculated from particle data
    const Real ugx = IC::ic[icv].fluid.velocity[0];
    const Real vgx = IC::ic[icv].fluid.velocity[1];
    const Real wgx = IC::ic[icv].fluid.velocity[2];

    const int istart = amrex::max(slo[0], i_w);
    const int jstart = amrex::max(slo[1], j_s);
    const int kstart = amrex::max(slo[2], k_b);
    const int iend   = amrex::min(shi[0], i_e);
    const int jend   = amrex::min(shi[1], j_n);
    const int kend   = amrex::min(shi[2], k_t);

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      amrex::ParallelFor(box1, [velocity, ugx]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          { velocity(i,j,k,0) = ugx; });

      if(slo[0] < domlo[0] and domlo[0] == istart)
      {
        const IntVect low2(slo[0], jstart, kstart), hi2(istart-1, jend, kend);
        const Box box2(low2, hi2);
        amrex::ParallelFor(box2, [velocity, ugx]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            { velocity(i,j,k,0) = ugx; });
      }

      if(shi[0] > domhi[0] and domhi[0] == iend)
      {
        const IntVect low3(iend+1, jstart, kstart), hi3(shi[0], jend, kend);
        const Box box3(low3, hi3);
        amrex::ParallelFor(box3, [velocity, ugx]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            { velocity(i,j,k,0) = ugx; });
      }
    }

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      amrex::ParallelFor(box1, [velocity, vgx]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          { velocity(i,j,k,1) = vgx; });

      if (slo[1] < domlo[1] and domlo[1] == jstart)
      {
        const IntVect low2(istart, slo[1], kstart), hi2(iend, jstart-1, kend);
        const Box box2(low2, hi2);
        amrex::ParallelFor(box2, [velocity, vgx]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            { velocity(i,j,k,1) = vgx; });
      }

      if (shi[1] > domhi[1] and domhi[1] == jend)
      {
        const IntVect low3(istart, jend+1, kstart), hi3(iend, shi[1], kend);
        const Box box3(low3, hi3);
        amrex::ParallelFor(box3, [velocity, vgx]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            { velocity(i,j,k,1) = vgx; });
      }
    }

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);
      amrex::ParallelFor(box1, [velocity, wgx]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          { velocity(i,j,k,2) = wgx; });

      if (slo[2] < domlo[2] and domlo[2] == kstart)
      {
        const IntVect low2(istart, jstart, slo[2]), hi2(iend, jend, kstart-1);
        const Box box2(low2, hi2);
        amrex::ParallelFor(box2, [velocity, wgx]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            { velocity(i,j,k,2) = wgx; });
      }

      if (shi[2] > domhi[2] and domhi[2] == kend)
      {
        const IntVect low3(istart, jstart, kend+1), hi3(iend, jend, shi[2]);
        const Box box3(low3, hi3);
        amrex::ParallelFor(box3, [velocity, wgx]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            { velocity(i,j,k,2) = wgx; });
      }
    }
  }
}

//!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//!                                                                      !
//!  Subroutine: SET_IC_TEMP                                             !
//!                                                                      !
//!  Purpose: Set fluid temperature initial conditions.                  !
//!                                                                      !
//!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void set_ic_temp (const Box& sbx,
                  const Box& domain,
                  const Real dx,
                  const Real dy,
                  const Real dz,
                  const GpuArray<Real, 3>& plo,
                  FArrayBox& T_g_fab)
{
  const IntVect slo(sbx.loVect());
  const IntVect shi(sbx.hiVect());

  const IntVect domlo(domain.loVect());
  const IntVect domhi(domain.hiVect());

  Array4<Real> const& T_g = T_g_fab.array();

  // Set the initial conditions.
  for(int icv(0); icv < IC::ic.size(); ++icv)
  {

    int i_w(0), j_s(0), k_b(0);
    int i_e(0), j_n(0), k_t(0);

    calc_cell_ic(dx, dy, dz,
                 IC::ic[icv].region->lo(),
                 IC::ic[icv].region->hi(),
                 plo.data(),
                 i_w, i_e, j_s, j_n, k_b, k_t);

    // Use the volume fraction already calculated from particle data
    const Real temperature = IC::ic[icv].fluid.temperature;

    const int istart = amrex::max(slo[0], i_w);
    const int jstart = amrex::max(slo[1], j_s);
    const int kstart = amrex::max(slo[2], k_b);
    const int iend   = amrex::min(shi[0], i_e);
    const int jend   = amrex::min(shi[1], j_n);
    const int kend   = amrex::min(shi[2], k_t);

    {
      {
        const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
        const Box box1(low1, hi1);

        ParallelFor(box1, [T_g,temperature]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { T_g(i,j,k) = temperature; });

        if(slo[0] < domlo[0] and domlo[0] == istart)
        {
          const IntVect low2(slo[0], jstart, kstart), hi2(istart-1, jend, kend);
          const Box box2(low2, hi2);
          ParallelFor(box2, [T_g,temperature]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          { T_g(i,j,k) = temperature; });
        }

        if(shi[0] > domhi[0] and domhi[0] == iend)
        {
          const IntVect low3(iend+1, jstart, kstart), hi3(shi[0], jend, kend);
          const Box box3(low3, hi3);
          ParallelFor(box3, [T_g,temperature]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          { T_g(i,j,k) = temperature; });
        }
      }

      {
        const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
        const Box box1(low1, hi1);

        ParallelFor(box1, [T_g,temperature]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { T_g(i,j,k) = temperature; });

        if (slo[1] < domlo[1] and domlo[1] == jstart)
        {
          const IntVect low2(istart, slo[1], kstart), hi2(iend, jstart-1, kend);
          const Box box2(low2, hi2);
          ParallelFor(box2, [T_g,temperature]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          { T_g(i,j,k) = temperature; });
        }

        if (shi[1] > domhi[1] and domhi[1] == jend)
        {
          const IntVect low3(istart, jend+1, kstart), hi3(iend, shi[1], kend);
          const Box box3(low3, hi3);
          ParallelFor(box3, [T_g,temperature]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          { T_g(i,j,k) = temperature; });
        }
      }

      {
        const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
        const Box box1(low1, hi1);
        ParallelFor(box1, [T_g,temperature]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { T_g(i,j,k) = temperature; });

        if (slo[2] < domlo[2] and domlo[2] == kstart)
        {
          const IntVect low2(istart, jstart, slo[2]), hi2(iend, jend, kstart-1);
          const Box box2(low2, hi2);

          ParallelFor(box2, [T_g,temperature]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          { T_g(i,j,k) = temperature; });
        }

        if (shi[2] > domhi[2] and domhi[2] == kend)
        {
          const IntVect low3(istart, jstart, kend+1), hi3(iend, jend, shi[2]);
          const Box box3(low3, hi3);

          ParallelFor(box3, [T_g,temperature]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          { T_g(i,j,k) = temperature; });
        }
      }
    }
  }
}

//!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//!                                                                      !
//!  Subroutine: SET_IC_SPECIES_G                                        !
//!                                                                      !
//!  Purpose: Set fluid species mass fractions initial conditions.       !
//!                                                                      !
//!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void set_ic_species_g (const Box& sbx,
                       const Box& domain,
                       const Real dx,
                       const Real dy,
                       const Real dz,
                       const GpuArray<Real, 3>& plo,
                       FArrayBox& X_gk_fab)
{
  const IntVect slo(sbx.loVect());
  const IntVect shi(sbx.hiVect());

  const IntVect domlo(domain.loVect());
  const IntVect domhi(domain.hiVect());

  Array4<Real> const& X_gk = X_gk_fab.array();

  const int nspecies_g = X_gk_fab.nComp();

  // Set the initial conditions.
  for(int icv(0); icv < IC::ic.size(); ++icv)
  {
    int i_w(0), j_s(0), k_b(0);
    int i_e(0), j_n(0), k_t(0);

    calc_cell_ic(dx, dy, dz,
                 IC::ic[icv].region->lo(), IC::ic[icv].region->hi(),
                 plo.data(),
                 i_w, i_e, j_s, j_n, k_b, k_t);

    // Get the initial condition values
    Gpu::DeviceVector< Real> mass_fractions_d(nspecies_g);
    Gpu::HostVector  < Real> mass_fractions_h(nspecies_g);
    for (int n(0); n < nspecies_g; n++) {
      mass_fractions_h[n] = IC::ic[icv].fluid.species[n].mass_fraction;
    }
    Gpu::copyAsync(Gpu::hostToDevice, mass_fractions_h.begin(), mass_fractions_h.end(),
                   mass_fractions_d.begin());

    Real* p_mass_fractions = mass_fractions_d.data();

    const int istart = std::max(slo[0], i_w);
    const int jstart = std::max(slo[1], j_s);
    const int kstart = std::max(slo[2], k_b);
    const int iend   = std::min(shi[0], i_e);
    const int jend   = std::min(shi[1], j_n);
    const int kend   = std::min(shi[2], k_t);

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      ParallelFor(box1, nspecies_g, [X_gk,p_mass_fractions]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      { X_gk(i,j,k,n) = p_mass_fractions[n]; });

      if(slo[0] < domlo[0] and domlo[0] == istart)
      {
        const IntVect low2(slo[0], jstart, kstart), hi2(istart-1, jend, kend);
        const Box box2(low2, hi2);
        ParallelFor(box2, nspecies_g, [X_gk,p_mass_fractions]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        { X_gk(i,j,k,n) = p_mass_fractions[n]; });
      }

      if(shi[0] > domhi[0] and domhi[0] == iend)
      {
        const IntVect low3(iend+1, jstart, kstart), hi3(shi[0], jend, kend);
        const Box box3(low3, hi3);
        ParallelFor(box3, nspecies_g, [X_gk,p_mass_fractions]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        { X_gk(i,j,k,n) = p_mass_fractions[n]; });
      }
    }

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      ParallelFor(box1, nspecies_g, [X_gk,p_mass_fractions]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      { X_gk(i,j,k,n) = p_mass_fractions[n]; });

      if (slo[1] < domlo[1] and domlo[1] == jstart)
      {
        const IntVect low2(istart, slo[1], kstart), hi2(iend, jstart-1, kend);
        const Box box2(low2, hi2);
        ParallelFor(box2, nspecies_g, [X_gk,p_mass_fractions]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        { X_gk(i,j,k,n) = p_mass_fractions[n]; });
      }

      if (shi[1] > domhi[1] and domhi[1] == jend)
      {
        const IntVect low3(istart, jend+1, kstart), hi3(iend, shi[1], kend);
        const Box box3(low3, hi3);
        ParallelFor(box3, nspecies_g, [X_gk,p_mass_fractions]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        { X_gk(i,j,k,n) = p_mass_fractions[n]; });
      }
    }

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);
      ParallelFor(box1, nspecies_g, [X_gk,p_mass_fractions]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      { X_gk(i,j,k,n) = p_mass_fractions[n]; });

      if (slo[2] < domlo[2] and domlo[2] == kstart)
      {
        const IntVect low2(istart, jstart, slo[2]), hi2(iend, jend, kstart-1);
        const Box box2(low2, hi2);

        ParallelFor(box2, nspecies_g, [X_gk,p_mass_fractions]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        { X_gk(i,j,k,n) = p_mass_fractions[n]; });
      }

      if (shi[2] > domhi[2] and domhi[2] == kend)
      {
        const IntVect low3(istart, jstart, kend+1), hi3(iend, jend, shi[2]);
        const Box box3(low3, hi3);

        ParallelFor(box3, nspecies_g, [X_gk,p_mass_fractions]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        { X_gk(i,j,k,n) = p_mass_fractions[n]; });
      }
    }

    Gpu::synchronize();

  }
}
