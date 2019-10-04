#include <mfix_init_fluid.hpp>
#include <fld_constants_mod_F.H>
#include <mfix_calc_mu_g.hpp>
#include <ic_mod_F.H>
#include <param_mod_F.H>
#include <calc_cell_F.H>

using namespace amrex;

// Forward declarations
void set_ic(const Box& sbx, const Box& domain, const Real dx, const Real dy,
    const Real dz, FArrayBox& vel_g_fab);

void init_helix(const Box& bx, const Box& domain, FArrayBox& vel_g_fab,
    const Real dx, const Real dy, const Real dz);

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: init_fluid                                              !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   
void init_fluid(const Box& sbx,
                const Box& bx,
                const Box& domain,
                const FArrayBox& ep_g_fab, 
                      FArrayBox& ro_g_fab, 
                      FArrayBox& trac_fab, 
                const FArrayBox& p_g, 
                      FArrayBox& vel_g_fab,
                      FArrayBox& mu_g_fab, 
                const Real dx,
                const Real dy,
                const Real dz,
                const Real xlength,
                const Real ylength,
                const Real zlength)
{
      // Set user specified initial conditions (IC)
      set_ic(sbx, domain, dx, dy, dz, vel_g_fab);

      // init_periodic_vortices (bx, domain, vel_g_fab, dx, dy, dz);
      // init_helix (bx, domain, vel_g_fab, dx, dy, dz);

      // Set the initial fluid density and viscosity
      Array4<Real> const& ro_g = ro_g_fab.array();
      Array4<Real> const& trac = trac_fab.array();

      const Real ro_g0  = get_ro_g0();
      const Real trac_0 = get_trac0();

      AMREX_HOST_DEVICE_FOR_3D(sbx, i, j, k, {ro_g(i,j,k) = ro_g0;});
      AMREX_HOST_DEVICE_FOR_3D(sbx, i, j, k, {trac(i,j,k) = trac_0;});

      Gpu::streamSynchronize();

      calc_mu_g(bx, mu_g_fab);
}

void init_helix(const Box& bx,
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
      AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
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
      AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
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
      AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
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

  Gpu::streamSynchronize();
}

void init_periodic_vortices(const Box& bx,
                            const Box& domain,
                            FArrayBox& vel_g_fab,
                            const Real dx,
                            const Real dy,
                            const Real dz)
{
//  const amrex::Real twopi = 8. * std::atan(1);
  const amrex::Real twopi = 2 * M_PI;

  Array4<Real> const& velocity = vel_g_fab.array();

  const int plane = 1;

  switch (plane)
  {
    case 1:  // x-y plane
      // x-direction
      AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
      {
        Real x = (Real(i) + .5) * dx;
        Real y = (Real(j) + .5) * dy;

        velocity(i,j,k,0) = std::tanh(30*(.25 - std::abs(y-.5)));
        velocity(i,j,k,1) = .05 * std::sin(twopi*x);
        velocity(i,j,k,2) = 0.;
      });
      break;

    case 2:  // x-z plane
      // x-direction
      AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
      {
        Real x = (Real(i) + .5) * dx;
        Real z = (Real(k) + .5) * dz;

        velocity(i,j,k,0) = std::tanh(30*(.25 - std::abs(z-.5)));
        velocity(i,j,k,1) = 0.;
        velocity(i,j,k,2) = .05 * std::sin(twopi*x);
      });
      break;

    case 3:  // y-z plane
      // x-direction
      AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
      {
        Real y = (Real(j) + .5) * dy;
        Real z = (Real(k) + .5) * dz;

        velocity(i,j,k,0) = 0.;
        velocity(i,j,k,1) = std::tanh(30*(.25 - std::abs(z-.5)));
        velocity(i,j,k,2) = .05 * std::sin(twopi*y);
      });
      break;

    default:
      amrex::Abort("Error: wrong plane number");
      break;
  }

  Gpu::streamSynchronize();
}

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: init_fluid_restart                                      !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void init_fluid_restart(const Box& bx,
                        FArrayBox& mu_g_fab)
{
  calc_mu_g(bx, mu_g_fab);
}

//!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//!                                                                      !
//!  Subroutine: SET_IC                                                  !
//!                                                                      !
//!  Purpose: This module sets all the initial conditions.               !
//!                                                                      !
//!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void set_ic(const Box& sbx,
            const Box& domain,
            const Real dx,
            const Real dy,
            const Real dz,
            FArrayBox& vel_g_fab)
{
  const IntVect slo(sbx.loVect());
  const IntVect shi(sbx.hiVect());

  const IntVect domlo(domain.loVect());
  const IntVect domhi(domain.hiVect());

  Array4<Real> const& velocity = vel_g_fab.array();

  // Set the initial conditions.
  for(int icv(1); icv <= get_dim_ic(); ++icv)
  {
    if (ic_defined_cpp(icv))
    {
      int i_w(0), j_s(0), k_b(0);
      int i_e(0), j_n(0), k_t(0);

      calc_cell_ic(dx, dy, dz, get_ic_x_w(icv), get_ic_y_s(icv),
                   get_ic_z_b(icv), get_ic_x_e(icv), get_ic_y_n(icv),
                   get_ic_z_t(icv), i_w, i_e, j_s, j_n, k_b, k_t);

      // Use the volume fraction already calculated from particle data
      const Real ugx = get_ic_u_g(icv);
      const Real vgx = get_ic_v_g(icv);
      const Real wgx = get_ic_w_g(icv);

      const int istart = std::max(slo[0], i_w);
      const int jstart = std::max(slo[1], j_s);
      const int kstart = std::max(slo[2], k_b);
      const int iend   = std::min(shi[0], i_e);
      const int jend   = std::min(shi[1], j_n);
      const int kend   = std::min(shi[2], k_t);

      if (is_defined_db_cpp(ugx))
      {
        const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
        const Box box1(low1, hi1);
        AMREX_HOST_DEVICE_FOR_3D(box1, i, j, k, {velocity(i,j,k,0) = ugx;});

        if(slo[0] < domlo[0] and domlo[0] == istart)
        {
          const IntVect low2(slo[0], jstart, kstart), hi2(istart-1, jend, kend);
          const Box box2(low2, hi2);
          AMREX_HOST_DEVICE_FOR_3D(box2, i, j, k, {velocity(i,j,k,0) = ugx;});
        }

        if(shi[0] > domhi[0] and domhi[0] == iend)
        {
          const IntVect low3(iend+1, jstart, kstart), hi3(shi[0], jend, kend);
          const Box box3(low3, hi3);
          AMREX_HOST_DEVICE_FOR_3D(box3, i, j, k, {velocity(i,j,k,0) = ugx;});
        }
      }

      if (is_defined_db_cpp(vgx))
      {
        const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
        const Box box1(low1, hi1);
        AMREX_HOST_DEVICE_FOR_3D(box1, i, j, k, {velocity(i,j,k,1) = vgx;});

        if (slo[1] < domlo[1] and domlo[1] == jstart)
        {
          const IntVect low2(istart, slo[1], kstart), hi2(iend, jstart-1, kend);
          const Box box2(low2, hi2);
          AMREX_HOST_DEVICE_FOR_3D(box2, i, j, k, {velocity(i,j,k,1) = vgx;});
        }

        if (shi[1] > domhi[1] and domhi[1] == jend)
        {
          const IntVect low3(istart, jend+1, kstart), hi3(iend, shi[1], kend);
          const Box box3(low3, hi3);
          AMREX_HOST_DEVICE_FOR_3D(box3, i, j, k, {velocity(i,j,k,1) = vgx;});
        }
      }

      if (is_defined_db_cpp(wgx))
      {
        const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
        const Box box1(low1, hi1);
        AMREX_HOST_DEVICE_FOR_3D(box1, i, j, k, {velocity(i,j,k,2) = wgx;});

        if (slo[2] < domlo[2] and domlo[2] == kstart)
        {
          const IntVect low2(istart, jstart, slo[2]), hi2(iend, jend, kstart-1);
          const Box box2(low2, hi2);
          AMREX_HOST_DEVICE_FOR_3D(box2, i, j, k, {velocity(i,j,k,2) = wgx;});
        }

        if (shi[2] > domhi[2] and domhi[2] == kend)
        {
          const IntVect low3(istart, jstart, kend+1), hi3(iend, jend, shi[2]);
          const Box box3(low3, hi3);
          AMREX_HOST_DEVICE_FOR_3D(box3, i, j, k, {velocity(i,j,k,2) = wgx;});
        }
      }
    }

    Gpu::streamSynchronize();
  }
}
