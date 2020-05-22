#ifndef MFIX_DEM_PARMS_H_
#define MFIX_DEM_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Vector.H>

#include <MFIX_SPECIES_Parms.H>

namespace DEM
{
    constexpr int NMAX = 10;
    extern AMREX_GPU_DEVICE_MANAGED int NPHASE;

    enum COLLISIONMODEL {
        LSD
    };

    extern AMREX_GPU_DEVICE_MANAGED COLLISIONMODEL CollisionModel;

    extern int solve;

    extern AMREX_GPU_DEVICE_MANAGED amrex::Real dtsolid;

    // tangential component of spring constant
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real kt;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real kt_w;

    // normal component of spring constant
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real kn;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real kn_w;

    // friction coefficients for particle-particle and particle-wall collisions
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real mew;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real mew_w;

    // assumed to be a constant factor relating these for linear spring-dashpot
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real kt_fac;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real kt_w_fac;

    // normal and tangential components of the damping coefficients
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real etan[NMAX][NMAX];
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real etan_w[NMAX];

    extern AMREX_GPU_DEVICE_MANAGED amrex::Real etat[NMAX][NMAX];
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real etat_w[NMAX];

    // coefficients of restitution, normal and tangential
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real en[NMAX][NMAX];
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real en_w[NMAX];


    extern AMREX_GPU_DEVICE_MANAGED amrex::Real eta_fac;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real eta_w_fac;


    extern AMREX_GPU_DEVICE_MANAGED amrex::Real small_number;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real large_number;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real eps;

    extern AMREX_GPU_DEVICE_MANAGED amrex::Real neighborhood;

    extern amrex::Vector<std::string> names;

    // Flag to solve species DEM particles equations
    extern int solve_species;

    // Particle species
    extern amrex::Vector<std::string> species_dem;

    // Particle species fractions
    extern amrex::Vector<amrex::Real> spec_frac_dem;

    // Total number of dem species
    extern int AMREX_GPU_DEVICE_MANAGED nspecies_dem;

    // Coarse-grain DEM
    extern int AMREX_GPU_DEVICE_MANAGED cg_dem;

    struct INPUT_DIST_t {
      std::string distribution;
      amrex::Real mean; // Mean value of distriubtion (or constant value)
      amrex::Real std;  // Standard deivation
      amrex::Real min;  // Minimum value of distribution
      amrex::Real max;  // Max value of distribution

      INPUT_DIST_t() :
        mean(-1.0),
        std(-1.0),
        min(-1.0),
        max(-1.0) { }

    };

    struct DEM_t {

      amrex::Real volfrac;
      amrex::Real statwt;
      amrex::Real temperature;
      amrex::Vector<amrex::Real> velocity;

      SPECIES::SPECIES_t species;

      INPUT_DIST_t diameter;
      INPUT_DIST_t density;

      DEM_t() :
        volfrac(0.0),
        statwt(0.0),
        temperature(0.0),
        velocity(0),
        species(DEM::nspecies_dem),
        diameter(),
        density()
      {}
    };

    void Initialize ();

}

#endif
