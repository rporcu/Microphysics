#ifndef MFIX_DEM_PARMS_H_
#define MFIX_DEM_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

#include <mfix_solids_parms.H>

namespace DEM
{
    extern int NPHASE;

    enum COLLISIONMODEL {
        LSD
    };

    extern COLLISIONMODEL CollisionModel;

    extern int solve;

    extern amrex::Real dtsolid;

    // tangential component of spring constant
    extern amrex::Real kt;
    extern amrex::Real kt_w;

    // normal component of spring constant
    extern amrex::Real kn;
    extern amrex::Real kn_w;

    // friction coefficients for particle-particle and particle-wall collisions
    extern amrex::Real mew;
    extern amrex::Real mew_w;

    // assumed to be a constant factor relating these for linear spring-dashpot
    extern amrex::Real kt_fac;
    extern amrex::Real kt_w_fac;

    template <typename ARRAY_TYPE>
    struct DeviceArray
    {
        typedef ARRAY_TYPE array_type;

        template <typename... Is>
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::Real& operator() (Is... index) const noexcept {
            return (*a)(index...);
        }

        void alloc () {
#ifdef AMREX_USE_GPU
            auto p = (ARRAY_TYPE*)amrex::The_Device_Arena()->alloc(sizeof(ARRAY_TYPE));
            amrex::single_task(amrex::Gpu::gpuStream(), [p] AMREX_GPU_DEVICE
            {
                new (p) ARRAY_TYPE;
            });
            amrex::Gpu::synchronize();
            a = p;
#else
            a = new ARRAY_TYPE;
#endif
        }

        void free () {
#ifdef AMREX_USE_GPU
            amrex::single_task(amrex::Gpu::gpuStream(), [p=a] AMREX_GPU_DEVICE
            {
                p->~ARRAY_TYPE();
            });
            amrex::Gpu::synchronize();
            amrex::The_Device_Arena()->free(a);
#else
            delete a;
#endif
        }

        ARRAY_TYPE* arrayPtr () { return a; }

    private:
        ARRAY_TYPE * a = nullptr;
    };

    using A2D = DeviceArray<amrex::Array2D<amrex::Real,0,SOLIDS::NMAX-1,0,SOLIDS::NMAX-1,amrex::Order::C> >;
    using A1D = DeviceArray<amrex::Array1D<amrex::Real,0,SOLIDS::NMAX-1> >;

    // normal and tangential components of the damping coefficients
    extern A2D etan;
    extern A1D etan_w;

    extern A2D etat;
    extern A1D etat_w;

    // coefficients of restitution, normal and tangential
    extern A2D en;
    extern A1D en_w;

    extern amrex::Real eta_fac;
    extern amrex::Real eta_w_fac;

    extern amrex::Real small_number;
    extern amrex::Real large_number;
    extern amrex::Real eps;

    extern amrex::Real neighborhood;

    // Flag to solve species DEM particles equations
    extern int solve_species;

    // Particle species
    extern amrex::Vector<std::string> species_dem;

    // Total number of dem species
    extern int nspecies_dem;

    // Coarse-grain DEM
    extern int cg_dem;


    void Initialize ();
    void Finalize ();

}

#endif
