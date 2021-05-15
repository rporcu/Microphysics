#include <mfix_init_fluid.H>

#include <mfix_calc_cell.H>
#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>
#include <mfix_ic_parms.H>

using namespace amrex;

// Forward declarations
void set_ic_vel (const Box& sbx, const Box& domain,
                 const Real dx, const Real dy, const Real dz,
                 const GpuArray<Real, 3>& plo, FArrayBox& vel_g_fab);

void set_ic_temp (const Box& sbx, const Box& domain,
                  const Real dx, const Real dy, const Real dz,
                  const GpuArray<Real, 3>& plo, FArrayBox& T_g_fab,
                  FArrayBox& h_g_fab, FArrayBox& X_gk_fab,
                  const int nspecies_g, const int is_mixture,
                  FluidPhase& fluid);

void set_ic_species_g (const Box& sbx, const Box& domain,
                       const Real dx, const Real dy, const Real dz,
                       const GpuArray<Real, 3>& plo, FArrayBox& X_gk_fab);

void set_ic_pressure_g (const Box& sbx, const Box& domain,
                        const Real dx, const Real dy, const Real dz,
                        const GpuArray<Real, 3>& plo, FArrayBox& pressure_g_fab,
                        FArrayBox& ro_g_fab, FArrayBox& T_g_fab,
                        FArrayBox& X_gk_fab, const Real R,
                        const int nspecies_g, const int is_mixture,
                        FluidPhase& fluid);

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
                 const int advect_fluid_species,
                 const int& idealgas_constraint,
                 FluidPhase& fluid)
{
  // Set user specified initial conditions (IC)
  // Set initial fluid velocity
  set_ic_vel(sbx, domain, dx, dy, dz, plo, (*ld.vel_g)[mfi]);

  // init_periodic_vortices (bx, domain, vel_g_fab, dx, dy, dz);
  // init_helix (bx, domain, vel_g_fab, dx, dy, dz);

  // Set the initial fluid density
  Array4<Real> const& ro_g = ld.ro_g->array(mfi);

  const Real ro_g0  = fluid.ro_g0;

  ParallelFor(sbx, [ro_g,ro_g0] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  { ro_g(i,j,k) = ro_g0; });

  // Set the initial fluid tracer
  Array4<Real> const& trac = ld.trac->array(mfi);

  const Real trac_0 = fluid.trac_0;

  ParallelFor(sbx, [trac,trac_0] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  { trac(i,j,k) = trac_0; });

  // Set the rank of the box
  Array4<Real> const& bx_proc = ld.ba_proc->array(mfi);

  const Real proc = Real(ParallelDescriptor::MyProc());

  ParallelFor(bx, [bx_proc, proc] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  { bx_proc(i,j,k) = proc; });

  // Fluid SPECIES Initialization
  if (advect_fluid_species) {
    // Set the initial fluid species mass fractions
    set_ic_species_g(sbx, domain, dx, dy, dz, plo, (*ld.X_gk)[mfi]);
  }

  // Set the initial fluid temperature
  if (advect_enthalpy) {
    if (fluid.is_a_mixture)
      set_ic_temp(sbx, domain, dx, dy, dz, plo, (*ld.T_g)[mfi], (*ld.h_g)[mfi],
          (*ld.X_gk)[mfi], fluid.nspecies, fluid.is_a_mixture, fluid);
    else {
      FArrayBox empty_fab;
      set_ic_temp(sbx, domain, dx, dy, dz, plo, (*ld.T_g)[mfi], (*ld.h_g)[mfi],
          empty_fab, fluid.nspecies, fluid.is_a_mixture, fluid);
    }
  }

  if (test_tracer_conservation)
    init_periodic_tracer(bx, domain, (*ld.vel_g)[mfi], (*ld.trac)[mfi], dx, dy, dz);

  if (idealgas_constraint == IdealGasConstraint::ClosedSystem) {
    // Set initial thermodynamic pressure according to the equation of state
    // p_g = ro_g * R * T_g / MW_g
    set_ic_pressure_g(sbx, domain, dx, dy, dz, plo, (*ld.pressure_g)[mfi],
        (*ld.ro_g)[mfi], (*ld.T_g)[mfi], (*ld.X_gk)[mfi], fluid.R,
        fluid.nspecies, fluid.is_a_mixture, fluid);
  }

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

  const Real fac = 0.01;

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
  const Real twopi = 2. * M_PI;

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
    const Real twopi(2. * M_PI);

    Array4<Real> const& trac = trac_fab.array();
    Array4<Real> const&  vel = vel_g_fab.array();

    int dir(0);
    const Real A(1.0);

    Real L(0); // domain size
    Real C(0); // sin coefficient

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
                            const int advect_fluid_species,
                            FluidPhase& fluid)
{
  const int fluid_is_a_mixture = fluid.is_a_mixture;
  const int advect_species_enthalpy = advect_fluid_species && advect_enthalpy;

  Array4<Real> nullArray4;

  Array4<Real> const& X_gk = fluid_is_a_mixture ? (ld.X_gk)->array(mfi) : nullArray4;
  Array4<Real> const& T_g  = advect_enthalpy ? (ld.T_g)->array(mfi) : nullArray4;
  Array4<Real> const& h_g  = advect_enthalpy ? (ld.h_g)->array(mfi) : nullArray4;

  const int nspecies_g = fluid.nspecies;

  const Real T_ref = fluid.T_ref;

  auto& fluid_parms = *fluid.parameters;

  // Set the IC values
  amrex::ParallelFor(bx, [nspecies_g,T_g,h_g,X_gk,advect_enthalpy,fluid_parms,
      advect_fluid_species,advect_species_enthalpy,fluid_is_a_mixture,T_ref]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    // set initial fluid enthalpy and  specific enthalpy
    if (advect_enthalpy)
    {
      const Real Tg_loc = T_g(i,j,k);

      if (fluid_is_a_mixture)
      {
        Real h_g_sum(0);
        for (int n(0); n < nspecies_g; n++) {
          const Real h_gk = fluid_parms.calc_h_gk(Tg_loc,n);
          h_g_sum += X_gk(i,j,k,n) * h_gk;
        }

        h_g(i,j,k) = h_g_sum;
      }
      else {
        h_g(i,j,k) = fluid_parms.calc_h_g(Tg_loc);
      }
    }

  });
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

      if(slo[0] < domlo[0] && domlo[0] == istart)
      {
        const IntVect low2(slo[0], jstart, kstart), hi2(istart-1, jend, kend);
        const Box box2(low2, hi2);
        amrex::ParallelFor(box2, [velocity, ugx]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            { velocity(i,j,k,0) = ugx; });
      }

      if(shi[0] > domhi[0] && domhi[0] == iend)
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

      if (slo[1] < domlo[1] && domlo[1] == jstart)
      {
        const IntVect low2(istart, slo[1], kstart), hi2(iend, jstart-1, kend);
        const Box box2(low2, hi2);
        amrex::ParallelFor(box2, [velocity, vgx]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            { velocity(i,j,k,1) = vgx; });
      }

      if (shi[1] > domhi[1] && domhi[1] == jend)
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

      if (slo[2] < domlo[2] && domlo[2] == kstart)
      {
        const IntVect low2(istart, jstart, slo[2]), hi2(iend, jend, kstart-1);
        const Box box2(low2, hi2);
        amrex::ParallelFor(box2, [velocity, wgx]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            { velocity(i,j,k,2) = wgx; });
      }

      if (shi[2] > domhi[2] && domhi[2] == kend)
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
                  FArrayBox& T_g_fab,
                  FArrayBox& h_g_fab,
                  FArrayBox& X_gk_fab,
                  const int nspecies_g,
                  const int is_mixture,
                  FluidPhase& fluid)
{
  const IntVect slo(sbx.loVect());
  const IntVect shi(sbx.hiVect());

  const IntVect domlo(domain.loVect());
  const IntVect domhi(domain.hiVect());

  Array4<Real      > const& T_g  = T_g_fab.array();
  Array4<Real      > const& h_g  = h_g_fab.array();
  Array4<Real const> const& X_gk = is_mixture ?
    X_gk_fab.array() : Array4<const Real>();

  auto& fluid_parms = *fluid.parameters;

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

        ParallelFor(box1, [T_g,h_g,X_gk,temperature,nspecies_g,is_mixture,fluid_parms]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          T_g(i,j,k) = temperature;

          if (is_mixture) {
            Real h_g_sum(0);
            for (int n(0); n < nspecies_g; n++) {
              Real h_gk = fluid_parms.calc_h_gk(temperature,n);
              h_g_sum += X_gk(i,j,k,n) * h_gk;
            }
            h_g(i,j,k) = h_g_sum;
          } else {
            h_g(i,j,k) = fluid_parms.calc_h_g(temperature);
          }
        });

        if(slo[0] < domlo[0] && domlo[0] == istart)
        {
          const IntVect low2(slo[0], jstart, kstart), hi2(istart-1, jend, kend);
          const Box box2(low2, hi2);
          ParallelFor(box2, [T_g,h_g,X_gk,temperature,nspecies_g,is_mixture,fluid_parms]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            T_g(i,j,k) = temperature;

            if (is_mixture) {
              Real h_g_sum(0);
              for (int n(0); n < nspecies_g; n++) {
                Real h_gk = fluid_parms.calc_h_gk(temperature,n);
                h_g_sum += X_gk(i,j,k,n) * h_gk;
              }
              h_g(i,j,k) = h_g_sum;
            } else {
              h_g(i,j,k) = fluid_parms.calc_h_g(temperature);
            }
          });
        }

        if(shi[0] > domhi[0] && domhi[0] == iend)
        {
          const IntVect low3(iend+1, jstart, kstart), hi3(shi[0], jend, kend);
          const Box box3(low3, hi3);
          ParallelFor(box3, [T_g,h_g,X_gk,temperature,nspecies_g,is_mixture,fluid_parms]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            T_g(i,j,k) = temperature;

            if (is_mixture) {
              Real h_g_sum(0);
              for (int n(0); n < nspecies_g; n++) {
                Real h_gk = fluid_parms.calc_h_gk(temperature,n);
                h_g_sum += X_gk(i,j,k,n) * h_gk;
              }
              h_g(i,j,k) = h_g_sum;
            } else {
              h_g(i,j,k) = fluid_parms.calc_h_g(temperature);
            }
          });
        }
      }

      {
        const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
        const Box box1(low1, hi1);

        ParallelFor(box1, [T_g,h_g,X_gk,temperature,nspecies_g,is_mixture,fluid_parms]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          T_g(i,j,k) = temperature;

          if (is_mixture) {
            Real h_g_sum(0);
            for (int n(0); n < nspecies_g; n++) {
              Real h_gk = fluid_parms.calc_h_gk(temperature,n);
              h_g_sum += X_gk(i,j,k,n) * h_gk;
            }
            h_g(i,j,k) = h_g_sum;
          } else {
            h_g(i,j,k) = fluid_parms.calc_h_g(temperature);
          }
        });

        if (slo[1] < domlo[1] && domlo[1] == jstart)
        {
          const IntVect low2(istart, slo[1], kstart), hi2(iend, jstart-1, kend);
          const Box box2(low2, hi2);
          ParallelFor(box2, [T_g,h_g,X_gk,temperature,nspecies_g,is_mixture,fluid_parms]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            T_g(i,j,k) = temperature;

            if (is_mixture) {
              Real h_g_sum(0);
              for (int n(0); n < nspecies_g; n++) {
                Real h_gk = fluid_parms.calc_h_gk(temperature,n);
                h_g_sum += X_gk(i,j,k,n) * h_gk;
              }
              h_g(i,j,k) = h_g_sum;
            } else {
              h_g(i,j,k) = fluid_parms.calc_h_g(temperature);
            }
          });
        }

        if (shi[1] > domhi[1] && domhi[1] == jend)
        {
          const IntVect low3(istart, jend+1, kstart), hi3(iend, shi[1], kend);
          const Box box3(low3, hi3);
          ParallelFor(box3, [T_g,h_g,X_gk,temperature,nspecies_g,is_mixture,fluid_parms]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            T_g(i,j,k) = temperature;

            if (is_mixture) {
              Real h_g_sum(0);
              for (int n(0); n < nspecies_g; n++) {
                Real h_gk = fluid_parms.calc_h_gk(temperature,n);
                h_g_sum += X_gk(i,j,k,n) * h_gk;
              }
              h_g(i,j,k) = h_g_sum;
            } else {
              h_g(i,j,k) = fluid_parms.calc_h_g(temperature);
            }
          });
        }
      }

      {
        const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
        const Box box1(low1, hi1);
        ParallelFor(box1, [T_g,h_g,X_gk,temperature,nspecies_g,is_mixture,fluid_parms]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          T_g(i,j,k) = temperature;

          if (is_mixture) {
            Real h_g_sum(0);
            for (int n(0); n < nspecies_g; n++) {
              Real h_gk = fluid_parms.calc_h_gk(temperature,n);
              h_g_sum += X_gk(i,j,k,n) * h_gk;
            }
            h_g(i,j,k) = h_g_sum;
          } else {
            h_g(i,j,k) = fluid_parms.calc_h_g(temperature);
          }
        });

        if (slo[2] < domlo[2] && domlo[2] == kstart)
        {
          const IntVect low2(istart, jstart, slo[2]), hi2(iend, jend, kstart-1);
          const Box box2(low2, hi2);

          ParallelFor(box2, [T_g,h_g,X_gk,temperature,nspecies_g,is_mixture,fluid_parms]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            T_g(i,j,k) = temperature;

            if (is_mixture) {
              Real h_g_sum(0);
              for (int n(0); n < nspecies_g; n++) {
                Real h_gk = fluid_parms.calc_h_gk(temperature,n);
                h_g_sum += X_gk(i,j,k,n) * h_gk;
              }
              h_g(i,j,k) = h_g_sum;
            } else {
              h_g(i,j,k) = fluid_parms.calc_h_g(temperature);
            }
          });
        }

        if (shi[2] > domhi[2] && domhi[2] == kend)
        {
          const IntVect low3(istart, jstart, kend+1), hi3(iend, jend, shi[2]);
          const Box box3(low3, hi3);

          ParallelFor(box3, [T_g,h_g,X_gk,temperature,nspecies_g,is_mixture,fluid_parms]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            T_g(i,j,k) = temperature;

            if (is_mixture) {
              Real h_g_sum(0);
              for (int n(0); n < nspecies_g; n++) {
                Real h_gk = fluid_parms.calc_h_gk(temperature,n);
                h_g_sum += X_gk(i,j,k,n) * h_gk;
              }
              h_g(i,j,k) = h_g_sum;
            } else {
              h_g(i,j,k) = fluid_parms.calc_h_g(temperature);
            }
          });
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

      if(slo[0] < domlo[0] && domlo[0] == istart)
      {
        const IntVect low2(slo[0], jstart, kstart), hi2(istart-1, jend, kend);
        const Box box2(low2, hi2);
        ParallelFor(box2, nspecies_g, [X_gk,p_mass_fractions]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        { X_gk(i,j,k,n) = p_mass_fractions[n]; });
      }

      if(shi[0] > domhi[0] && domhi[0] == iend)
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

      if (slo[1] < domlo[1] && domlo[1] == jstart)
      {
        const IntVect low2(istart, slo[1], kstart), hi2(iend, jstart-1, kend);
        const Box box2(low2, hi2);
        ParallelFor(box2, nspecies_g, [X_gk,p_mass_fractions]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        { X_gk(i,j,k,n) = p_mass_fractions[n]; });
      }

      if (shi[1] > domhi[1] && domhi[1] == jend)
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

      if (slo[2] < domlo[2] && domlo[2] == kstart)
      {
        const IntVect low2(istart, jstart, slo[2]), hi2(iend, jend, kstart-1);
        const Box box2(low2, hi2);

        ParallelFor(box2, nspecies_g, [X_gk,p_mass_fractions]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        { X_gk(i,j,k,n) = p_mass_fractions[n]; });
      }

      if (shi[2] > domhi[2] && domhi[2] == kend)
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

//!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//!                                                                      !
//!  Subroutine: SET_IC_PRESSURE_G                                       !
//!                                                                      !
//!  Purpose: Set fluid thermodynamic pressure initial conditions.       !
//!                                                                      !
//!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void set_ic_pressure_g (const Box& sbx,
                        const Box& domain,
                        const Real dx,
                        const Real dy,
                        const Real dz,
                        const GpuArray<Real, 3>& plo,
                        FArrayBox& pressure_g_fab,
                        FArrayBox& ro_g_fab,
                        FArrayBox& T_g_fab,
                        FArrayBox& X_gk_fab,
                        const Real R,
                        const int nspecies_g,
                        const int is_mixture,
                        FluidPhase& fluid)
{
  const IntVect slo(sbx.loVect());
  const IntVect shi(sbx.hiVect());

  const IntVect domlo(domain.loVect());
  const IntVect domhi(domain.hiVect());

  Array4<Real      > const& pressure_g = pressure_g_fab.array();
  Array4<Real const> const& ro_g       = ro_g_fab.array();
  Array4<Real const> const& T_g        = T_g_fab.array();
  Array4<Real const> const& X_gk       = X_gk_fab.array();

  // IC values for MW_g
  const Real MW_g0 = fluid.MW_g0;

  Gpu::DeviceVector<Real> MW_gk0_d(nspecies_g);
  if (is_mixture)
    Gpu::copyAsync(Gpu::hostToDevice, fluid.MW_gk0.begin(), fluid.MW_gk0.end(), MW_gk0_d.begin());
  Real* p_MW_gk0 = is_mixture ? MW_gk0_d.data() : nullptr;

  // Set the initial conditions.
  for(int icv(0); icv < IC::ic.size(); ++icv)
  {
    int i_w(0), j_s(0), k_b(0);
    int i_e(0), j_n(0), k_t(0);

    calc_cell_ic(dx, dy, dz,
                 IC::ic[icv].region->lo(), IC::ic[icv].region->hi(),
                 plo.data(),
                 i_w, i_e, j_s, j_n, k_b, k_t);

    const int istart = std::max(slo[0], i_w);
    const int jstart = std::max(slo[1], j_s);
    const int kstart = std::max(slo[2], k_b);
    const int iend   = std::min(shi[0], i_e);
    const int jend   = std::min(shi[1], j_n);
    const int kend   = std::min(shi[2], k_t);

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      ParallelFor(box1, [pressure_g,ro_g,T_g,X_gk,MW_g0,p_MW_gk0,R,is_mixture,nspecies_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real MW_g(0);

        // set initial fluid molecular weight
        if (is_mixture) {
          for (int n(0); n < nspecies_g; n++)
            MW_g += X_gk(i,j,k,n) / p_MW_gk0[n];
          MW_g = 1. / MW_g;
        }
        else
          MW_g = MW_g0;

        pressure_g(i,j,k) = ro_g(i,j,k) * R * T_g(i,j,k) / MW_g;
      });

      if(slo[0] < domlo[0] && domlo[0] == istart)
      {
        const IntVect low2(slo[0], jstart, kstart), hi2(istart-1, jend, kend);
        const Box box2(low2, hi2);
        ParallelFor(box2, [pressure_g,ro_g,T_g,X_gk,MW_g0,p_MW_gk0,R,is_mixture,nspecies_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          Real MW_g(0);

          // set initial fluid molecular weight
          if (is_mixture) {
            for (int n(0); n < nspecies_g; n++)
              MW_g += X_gk(i,j,k,n) / p_MW_gk0[n];
            MW_g = 1. / MW_g;
          }
          else
            MW_g = MW_g0;

          pressure_g(i,j,k) = ro_g(i,j,k) * R * T_g(i,j,k) / MW_g;
        });
      }

      if(shi[0] > domhi[0] && domhi[0] == iend)
      {
        const IntVect low3(iend+1, jstart, kstart), hi3(shi[0], jend, kend);
        const Box box3(low3, hi3);
        ParallelFor(box3, [pressure_g,ro_g,T_g,X_gk,MW_g0,p_MW_gk0,R,is_mixture,nspecies_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          Real MW_g(0);

          // set initial fluid molecular weight
          if (is_mixture) {
            for (int n(0); n < nspecies_g; n++)
              MW_g += X_gk(i,j,k,n) / p_MW_gk0[n];
            MW_g = 1. / MW_g;
          }
          else
            MW_g = MW_g0;

          pressure_g(i,j,k) = ro_g(i,j,k) * R * T_g(i,j,k) / MW_g;
        });
      }
    }

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      ParallelFor(box1, [pressure_g,ro_g,T_g,X_gk,MW_g0,p_MW_gk0,R,is_mixture,nspecies_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real MW_g(0);

        // set initial fluid molecular weight
        if (is_mixture) {
          for (int n(0); n < nspecies_g; n++)
            MW_g += X_gk(i,j,k,n) / p_MW_gk0[n];
          MW_g = 1. / MW_g;
        }
        else
          MW_g = MW_g0;

        pressure_g(i,j,k) = ro_g(i,j,k) * R * T_g(i,j,k) / MW_g;
      });

      if (slo[1] < domlo[1] && domlo[1] == jstart)
      {
        const IntVect low2(istart, slo[1], kstart), hi2(iend, jstart-1, kend);
        const Box box2(low2, hi2);
        ParallelFor(box2, [pressure_g,ro_g,T_g,X_gk,MW_g0,p_MW_gk0,R,is_mixture,nspecies_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          Real MW_g(0);

          // set initial fluid molecular weight
          if (is_mixture) {
            for (int n(0); n < nspecies_g; n++)
              MW_g += X_gk(i,j,k,n) / p_MW_gk0[n];
            MW_g = 1. / MW_g;
          }
          else
            MW_g = MW_g0;

          pressure_g(i,j,k) = ro_g(i,j,k) * R * T_g(i,j,k) / MW_g;
        });
      }

      if (shi[1] > domhi[1] && domhi[1] == jend)
      {
        const IntVect low3(istart, jend+1, kstart), hi3(iend, shi[1], kend);
        const Box box3(low3, hi3);
        ParallelFor(box3, [pressure_g,ro_g,T_g,X_gk,MW_g0,p_MW_gk0,R,is_mixture,nspecies_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          Real MW_g(0);

          // set initial fluid molecular weight
          if (is_mixture) {
            for (int n(0); n < nspecies_g; n++)
              MW_g += X_gk(i,j,k,n) / p_MW_gk0[n];
            MW_g = 1. / MW_g;
          }
          else
            MW_g = MW_g0;

          pressure_g(i,j,k) = ro_g(i,j,k) * R * T_g(i,j,k) / MW_g;
        });
      }
    }

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);
      ParallelFor(box1, [pressure_g,ro_g,T_g,X_gk,MW_g0,p_MW_gk0,R,is_mixture,nspecies_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real MW_g(0);

        // set initial fluid molecular weight
        if (is_mixture) {
          for (int n(0); n < nspecies_g; n++)
            MW_g += X_gk(i,j,k,n) / p_MW_gk0[n];
          MW_g = 1. / MW_g;
        }
        else
          MW_g = MW_g0;

        pressure_g(i,j,k) = ro_g(i,j,k) * R * T_g(i,j,k) / MW_g;
      });

      if (slo[2] < domlo[2] && domlo[2] == kstart)
      {
        const IntVect low2(istart, jstart, slo[2]), hi2(iend, jend, kstart-1);
        const Box box2(low2, hi2);

        ParallelFor(box2, [pressure_g,ro_g,T_g,X_gk,MW_g0,p_MW_gk0,R,is_mixture,nspecies_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          Real MW_g(0);

          // set initial fluid molecular weight
          if (is_mixture) {
            for (int n(0); n < nspecies_g; n++)
              MW_g += X_gk(i,j,k,n) / p_MW_gk0[n];
            MW_g = 1. / MW_g;
          }
          else
            MW_g = MW_g0;

          pressure_g(i,j,k) = ro_g(i,j,k) * R * T_g(i,j,k) / MW_g;
        });
      }

      if (shi[2] > domhi[2] && domhi[2] == kend)
      {
        const IntVect low3(istart, jstart, kend+1), hi3(iend, jend, shi[2]);
        const Box box3(low3, hi3);

        ParallelFor(box3, [pressure_g,ro_g,T_g,X_gk,MW_g0,p_MW_gk0,R,is_mixture,nspecies_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          Real MW_g(0);

          // set initial fluid molecular weight
          if (is_mixture) {
            for (int n(0); n < nspecies_g; n++)
              MW_g += X_gk(i,j,k,n) / p_MW_gk0[n];
            MW_g = 1. / MW_g;
          }
          else
            MW_g = MW_g0;

          pressure_g(i,j,k) = ro_g(i,j,k) * R * T_g(i,j,k) / MW_g;
        });
      }
    }
  }
}
