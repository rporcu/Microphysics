#ifndef actual_network_H
#define actual_network_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <network_properties.H>

using namespace amrex;

AMREX_INLINE
void actual_network_init() {}

const std::string network_name = "droplet_evaporation";

namespace Rates {
    enum NetworkRates {
        evaporation = 1,
        NumRates = evaporation
    };
}

//namespace RHS {
//
//    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//    constexpr rhs_t rhs_data (int rate)
//    {
//        using namespace Species;
//        using namespace Rates;
//
//        rhs_t data;
//
//        switch (rate) {
//
//        case evaporation:
//            data.species_A = Liquid;
//            data.species_D = Vapor;
//
//            data.number_A = 1.0_rt;
//            data.number_D = 1.0_rt;
//
//            data.exponent_A = 1;
//            data.exponent_D = 1;
//
//            data.apply_identical_particle_factor = 0;
//            break;
//        }
//
//        return data;
//    }
//}

#endif
