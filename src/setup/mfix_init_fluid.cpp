#include <mfix_init_fluid.H>

#include <mfix_calc_cell.H>
#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>
#include <mfix_ic_parms.H>


using namespace amrex;


namespace init_fluid_aux {

// Forward declarations
void set_ic_vel (const Box& sbx, const Box& domain,
                 const Real dx, const Real dy, const Real dz,
                 const GpuArray<Real, 3>& plo, FArrayBox& vel_g_fab);

void set_ic_temp (const Box& sbx, const Box& domain,
                  const Real dx, const Real dy, const Real dz,
                  const GpuArray<Real, 3>& plo, FArrayBox& T_g_fab,
                  FArrayBox* h_g_fab, FArrayBox* X_gk_fab,
                  FluidPhase& fluid);

void set_ic_species_g (const Box& sbx, const Box& domain,
                       const Real dx, const Real dy, const Real dz,
                       const GpuArray<Real, 3>& plo, FArrayBox& X_gk_fab);

void set_ic_ro_g (const Box& sbx, const Box& domain,
                  const Real dx, const Real dy, const Real dz,
                  const GpuArray<Real, 3>& plo, FArrayBox& ro_g_fab);

void set_ic_thermo_p_g (const Box& sbx, const Box& domain,
                        const Real dx, const Real dy, const Real dz,
                        const GpuArray<Real, 3>& plo, FArrayBox& p_g_fab,
                        const FluidPhase& fluid);

void init_helix (const Box& bx, const Box& domain, FArrayBox& vel_g_fab,
                 const Real dx, const Real dy, const Real dz);

} // end namespace init_fluid_aux


using namespace init_fluid_aux;


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
                 const Real /*xlength*/,
                 const Real /*ylength*/,
                 const Real /*zlength*/,
                 const GpuArray<Real, 3>& plo,
                 bool test_tracer_conservation,
                 FluidPhase& fluid)
{
  // Set user specified initial conditions (IC)

  // **************************************************************************
  // Set initial fluid velocity
  // **************************************************************************
  set_ic_vel(sbx, domain, dx, dy, dz, plo, (*ld.vel_g)[mfi]);

  // init_periodic_vortices (bx, domain, vel_g_fab, dx, dy, dz);
  // init_helix (bx, domain, vel_g_fab, dx, dy, dz);

  // **************************************************************************
  // Set initial fluid tracer
  // **************************************************************************
  if (test_tracer_conservation) {

    init_periodic_tracer(bx, domain, (*ld.vel_g)[mfi], (*ld.trac)[mfi], dx, dy, dz);

  } else {

    const Real trac_0 = fluid.trac_0;
    Array4<Real> const& trac = ld.trac->array(mfi);
    ParallelFor(sbx, [trac,trac_0] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { trac(i,j,k) = trac_0; });
  }

  // **************************************************************************
  // Set initial fluid species mass fractions
  // **************************************************************************
  if (fluid.solve_species) {
    // Set the initial fluid species mass fractions
    set_ic_species_g(sbx, domain, dx, dy, dz, plo, (*ld.X_gk)[mfi]);
  }

  // ************************************************************************
  // Set initial fluid density
  // ************************************************************************
  set_ic_ro_g(sbx, domain, dx, dy, dz, plo, (*ld.ro_g)[mfi]);

  // **************************************************************************
  // Set initial fluid temperature
  // **************************************************************************
  if (fluid.solve_enthalpy ||
      (fluid.constraint_type == ConstraintType::IdealGasOpenSystem ||
       fluid.constraint_type == ConstraintType::IdealGasClosedSystem)) {

    FArrayBox* X_gk_fab = fluid.is_a_mixture ? &((*ld.X_gk)[mfi]) : nullptr;
    FArrayBox* h_g_fab = fluid.solve_enthalpy ? &((*ld.h_g)[mfi]) : nullptr;

    set_ic_temp(sbx, domain, dx, dy, dz, plo, (*ld.T_g)[mfi], h_g_fab, X_gk_fab, fluid);

    if (!fluid.solve_enthalpy) {
      ((*ld.T_go)[mfi]).copy((*ld.T_g)[mfi], 0, 0, 1);
    }
  }

  // ************************************************************************
  // Set initial fluid thermodynamic_pressure
  // ************************************************************************
  if (fluid.solve_enthalpy &&
      (fluid.constraint_type == ConstraintType::IdealGasOpenSystem ||
       fluid.constraint_type == ConstraintType::IdealGasClosedSystem)) {
    set_ic_thermo_p_g(sbx, domain, dx, dy, dz, plo, (*ld.thermodynamic_p_g)[mfi], fluid);
  }
}


namespace init_fluid_aux {

void init_helix (const Box& bx,
                 const Box& /*domain*/,
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

} // end namespace init_fluid_aux


void init_periodic_vortices (const Box& bx,
                             const Box& /*domain*/,
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

            Real x = (Real(i) + .5) * dx - .00073;
            Real y = (Real(j) + .5) * dy - .00123;

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
                            FluidPhase& fluid)
{
  const int run_on_device = Gpu::inLaunchRegion() ? 1 : 0;

  const EBFArrayBox& epg_fab = static_cast<EBFArrayBox const&>((*ld.ep_g)[mfi]);
  const EBCellFlagFab& flags = epg_fab.getEBCellFlagFab();

  const int solve_enthalpy = fluid.solve_enthalpy;
  const int fluid_is_a_mixture = fluid.is_a_mixture;

  Array4<Real> dummy_arr;

  Array4<Real> const& X_gk = fluid_is_a_mixture ? (ld.X_gk)->array(mfi) : dummy_arr;
  Array4<Real> const& T_g  = fluid.solve_enthalpy ? (ld.T_g)->array(mfi) : dummy_arr;
  Array4<Real> const& h_g  = fluid.solve_enthalpy ? (ld.h_g)->array(mfi) : dummy_arr;

  auto const& flags_arr = flags.const_array();

  const int nspecies_g = fluid.nspecies;

  auto& fluid_parms = *fluid.parameters;

  // Set the IC values
  amrex::ParallelFor(bx, [nspecies_g,T_g,h_g,X_gk,solve_enthalpy,fluid_parms,
      fluid_is_a_mixture,run_on_device,flags_arr]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

    // set initial fluid enthalpy and  specific enthalpy
    if (solve_enthalpy)
    {
      const Real Tg_loc = T_g(i,j,k);

      if (fluid_is_a_mixture)
      {
        Real h_g_sum(0);
        for (int n(0); n < nspecies_g; n++) {
          const Real h_gk = run_on_device ?
            fluid_parms.calc_h_gk<RunOn::Device>(Tg_loc, n, cell_is_covered) :
            fluid_parms.calc_h_gk<RunOn::Host>(Tg_loc, n, cell_is_covered);

          h_g_sum += X_gk(i,j,k,n) * h_gk;
        }

        h_g(i,j,k) = h_g_sum;
      }
      else {
        h_g(i,j,k) = run_on_device ?
          fluid_parms.calc_h_g<RunOn::Device>(Tg_loc, cell_is_covered) :
          fluid_parms.calc_h_g<RunOn::Host>(Tg_loc, cell_is_covered);
      }
    }

  });
}


namespace init_fluid_aux {

//!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
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
                  FArrayBox* h_g_fab,
                  FArrayBox* X_gk_fab,
                  FluidPhase& fluid)
{
  const int run_on_device = Gpu::inLaunchRegion() ? 1 : 0;

  const EBFArrayBox& Tg_EB_fab = static_cast<EBFArrayBox const&>(T_g_fab);
  const EBCellFlagFab& flags = Tg_EB_fab.getEBCellFlagFab();

  const IntVect slo(sbx.loVect());
  const IntVect shi(sbx.hiVect());

  const IntVect domlo(domain.loVect());
  const IntVect domhi(domain.hiVect());

  const int solve_enthalpy = fluid.solve_enthalpy;
  const int fluid_is_a_mixture = fluid.is_a_mixture;
  const int nspecies_g = fluid.nspecies;

  Array4<Real> dummy_arr;
  Array4<const Real> dummy_const_arr;

  Array4<Real      > const& T_g  = T_g_fab.array();
  Array4<Real      > const& h_g  = solve_enthalpy ? h_g_fab->array() : dummy_arr;
  Array4<Real const> const& X_gk = fluid_is_a_mixture ? X_gk_fab->array() : dummy_const_arr;

  auto const& flags_arr = flags.const_array();

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

    const Real temperature = IC::ic[icv].fluid.temperature;

    const int istart = amrex::max(slo[0], i_w);
    const int jstart = amrex::max(slo[1], j_s);
    const int kstart = amrex::max(slo[2], k_b);
    const int iend   = amrex::min(shi[0], i_e);
    const int jend   = amrex::min(shi[1], j_n);
    const int kend   = amrex::min(shi[2], k_t);

    // Define the function to be used on the different Box-es
    auto set_quantities = [T_g,h_g,X_gk,temperature,nspecies_g,fluid_is_a_mixture,
         fluid_parms,run_on_device,flags_arr,solve_enthalpy]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

      T_g(i,j,k) = temperature;

      if (solve_enthalpy) {
        if (fluid_is_a_mixture) {
          Real h_g_sum(0);
          for (int n(0); n < nspecies_g; n++) {
            Real h_gk = run_on_device ?
              fluid_parms.calc_h_gk<RunOn::Device>(temperature, n, cell_is_covered) :
              fluid_parms.calc_h_gk<RunOn::Host>(temperature, n, cell_is_covered);

            h_g_sum += X_gk(i,j,k,n) * h_gk;
          }
          h_g(i,j,k) = h_g_sum;
        } else {
          h_g(i,j,k) = run_on_device ?
            fluid_parms.calc_h_g<RunOn::Device>(temperature, cell_is_covered) :
            fluid_parms.calc_h_g<RunOn::Host>(temperature, cell_is_covered);
        }
      }
    };

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      ParallelFor(box1, [set_quantities] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { set_quantities(i,j,k); });

      if(slo[0] < domlo[0] && domlo[0] == istart)
      {
        const IntVect low2(slo[0], jstart, kstart), hi2(istart-1, jend, kend);
        const Box box2(low2, hi2);
        ParallelFor(box2, [set_quantities] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_quantities(i,j,k); });
      }

      if(shi[0] > domhi[0] && domhi[0] == iend)
      {
        const IntVect low3(iend+1, jstart, kstart), hi3(shi[0], jend, kend);
        const Box box3(low3, hi3);
        ParallelFor(box3, [set_quantities] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_quantities(i,j,k); });
      }
    }

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      ParallelFor(box1, [set_quantities] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { set_quantities(i,j,k); });

      if (slo[1] < domlo[1] && domlo[1] == jstart)
      {
        const IntVect low2(istart, slo[1], kstart), hi2(iend, jstart-1, kend);
        const Box box2(low2, hi2);
        ParallelFor(box2, [set_quantities] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_quantities(i,j,k); });
      }

      if (shi[1] > domhi[1] && domhi[1] == jend)
      {
        const IntVect low3(istart, jend+1, kstart), hi3(iend, shi[1], kend);
        const Box box3(low3, hi3);
        ParallelFor(box3, [set_quantities] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_quantities(i,j,k); });
      }
    }

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);
      ParallelFor(box1, [set_quantities] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { set_quantities(i,j,k); });

      if (slo[2] < domlo[2] && domlo[2] == kstart)
      {
        const IntVect low2(istart, jstart, slo[2]), hi2(iend, jend, kstart-1);
        const Box box2(low2, hi2);

        ParallelFor(box2, [set_quantities] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_quantities(i,j,k); });
      }

      if (shi[2] > domhi[2] && domhi[2] == kend)
      {
        const IntVect low3(istart, jstart, kend+1), hi3(iend, jend, shi[2]);
        const Box box3(low3, hi3);

        ParallelFor(box3, [set_quantities] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_quantities(i,j,k); });
      }
    }
  }
}

//!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
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
    Gpu::copy(Gpu::hostToDevice, mass_fractions_h.begin(), mass_fractions_h.end(),
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
//!  Purpose: Set fluid thermodynamic pressure initial conditions.       !
//!                                                                      !
//!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void set_ic_thermo_p_g (const Box& sbx,
                        const Box& domain,
                        const Real dx,
                        const Real dy,
                        const Real dz,
                        const GpuArray<Real, 3>& plo,
                        FArrayBox& thermodynamic_pressure_fab,
                        const FluidPhase& fluid)
{
  const IntVect slo(sbx.loVect());
  const IntVect shi(sbx.hiVect());

  const IntVect domlo(domain.loVect());
  const IntVect domhi(domain.hiVect());

  Array4<Real> const& thermo_p_g = thermodynamic_pressure_fab.array();

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

    const Real pressure = fluid.thermodynamic_pressure;

    // Define the function
    auto set_pressure = [thermo_p_g,pressure]
      AMREX_GPU_DEVICE (int i, int j, int k) -> void
    { thermo_p_g(i,j,k) = pressure; };

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      ParallelFor(box1, [set_pressure] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { set_pressure(i,j,k); });

      if(slo[0] < domlo[0] && domlo[0] == istart) {

        const IntVect low2(slo[0], jstart, kstart), hi2(istart-1, jend, kend);
        const Box box2(low2, hi2);

        ParallelFor(box2, [set_pressure] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_pressure(i,j,k); });
      }

      if(shi[0] > domhi[0] && domhi[0] == iend) {

        const IntVect low3(iend+1, jstart, kstart), hi3(shi[0], jend, kend);
        const Box box3(low3, hi3);

        ParallelFor(box3, [set_pressure] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_pressure(i,j,k); });
      }
    }

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      ParallelFor(box1, [set_pressure] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { set_pressure(i,j,k); });

      if (slo[1] < domlo[1] && domlo[1] == jstart) {
        const IntVect low2(istart, slo[1], kstart), hi2(iend, jstart-1, kend);
        const Box box2(low2, hi2);

        ParallelFor(box2, [set_pressure] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_pressure(i,j,k); });
      }

      if (shi[1] > domhi[1] && domhi[1] == jend) {
        const IntVect low3(istart, jend+1, kstart), hi3(iend, shi[1], kend);
        const Box box3(low3, hi3);

        ParallelFor(box3, [set_pressure] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_pressure(i,j,k); });
      }
    }

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      ParallelFor(box1, [set_pressure] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { set_pressure(i,j,k); });

      if (slo[2] < domlo[2] && domlo[2] == kstart) {
        const IntVect low2(istart, jstart, slo[2]), hi2(iend, jend, kstart-1);
        const Box box2(low2, hi2);

        ParallelFor(box2, [set_pressure] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_pressure(i,j,k); });
      }

      if (shi[2] > domhi[2] && domhi[2] == kend) {
        const IntVect low3(istart, jstart, kend+1), hi3(iend, jend, shi[2]);
        const Box box3(low3, hi3);

        ParallelFor(box3, [set_pressure] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_pressure(i,j,k); });
      }
    }
  }

}

//!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//!                                                                      !
//!  Purpose: Set fluid density                                          !
//!                                                                      !
//!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void set_ic_ro_g (const Box& sbx,
                  const Box& domain,
                  const Real dx,
                  const Real dy,
                  const Real dz,
                  const GpuArray<Real, 3>& plo,
                  FArrayBox& ro_g_fab)
{
  const IntVect slo(sbx.loVect());
  const IntVect shi(sbx.hiVect());

  const IntVect domlo(domain.loVect());
  const IntVect domhi(domain.hiVect());

  Array4<Real> const& ro_g = ro_g_fab.array();

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

    const Real density = IC::ic[icv].fluid.density;

    // Define the function
    auto set_density = [ro_g,density]
      AMREX_GPU_DEVICE (int i, int j, int k) -> void
    { ro_g(i,j,k) = density; };

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      ParallelFor(box1, [set_density] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { set_density(i,j,k); });

      if(slo[0] < domlo[0] && domlo[0] == istart) {

        const IntVect low2(slo[0], jstart, kstart), hi2(istart-1, jend, kend);
        const Box box2(low2, hi2);

        ParallelFor(box2, [set_density] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_density(i,j,k); });
      }

      if(shi[0] > domhi[0] && domhi[0] == iend) {

        const IntVect low3(iend+1, jstart, kstart), hi3(shi[0], jend, kend);
        const Box box3(low3, hi3);

        ParallelFor(box3, [set_density] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_density(i,j,k); });
      }
    }

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      ParallelFor(box1, [set_density] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { set_density(i,j,k); });

      if (slo[1] < domlo[1] && domlo[1] == jstart) {
        const IntVect low2(istart, slo[1], kstart), hi2(iend, jstart-1, kend);
        const Box box2(low2, hi2);

        ParallelFor(box2, [set_density] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_density(i,j,k); });
      }

      if (shi[1] > domhi[1] && domhi[1] == jend) {
        const IntVect low3(istart, jend+1, kstart), hi3(iend, shi[1], kend);
        const Box box3(low3, hi3);

        ParallelFor(box3, [set_density] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_density(i,j,k); });
      }
    }

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      ParallelFor(box1, [set_density] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { set_density(i,j,k); });

      if (slo[2] < domlo[2] && domlo[2] == kstart) {
        const IntVect low2(istart, jstart, slo[2]), hi2(iend, jend, kstart-1);
        const Box box2(low2, hi2);

        ParallelFor(box2, [set_density] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_density(i,j,k); });
      }

      if (shi[2] > domhi[2] && domhi[2] == kend) {
        const IntVect low3(istart, jstart, kend+1), hi3(iend, jend, shi[2]);
        const Box box3(low3, hi3);

        ParallelFor(box3, [set_density] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_density(i,j,k); });
      }
    }
  }
}

} // end namespace init_fluid_aux
