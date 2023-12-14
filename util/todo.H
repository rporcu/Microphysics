#ifndef TODO_H
#define TODO_H

#include <vector>
#include <string>

#include <AMReX_REAL.H>

static constexpr bool scale_system = false;
static constexpr bool call_eos_in_rhs = false;
static constexpr bool do_species_clip = false;
static constexpr bool use_number_densities = false;
static constexpr bool renormalize_abundances = false;
static constexpr bool integrate_energy = false;
static constexpr bool use_jacobian_caching = false;
static constexpr bool subtract_internal_energy = false;
static constexpr bool burner_verbose = false;

static constexpr int NumSpec = 1;
static constexpr int NumAux = 0;
static constexpr int ode_max_steps = 0;
static constexpr int jacobian = 0;

static constexpr amrex::Real small_x = 0;
static constexpr amrex::Real SMALL_X_SAFE = 0;
static constexpr amrex::Real react_boost = 0;
static constexpr amrex::Real MAX_TEMP = 0;
static constexpr amrex::Real X_reject_buffer = 0;
static constexpr amrex::Real atol_spec = 0;
static constexpr amrex::Real atol_enuc = 0;
static constexpr amrex::Real rtol_spec = 0;
static constexpr amrex::Real rtol_enuc = 0;
static constexpr amrex::Real ode_max_dt = 0;

static const std::vector<std::string> spec_names_cxx = {};
static const std::vector<std::string> short_spec_names_cxx = {};
static const std::vector<amrex::Real> zion = {};
static const std::vector<amrex::Real> aion = {};
static const std::vector<amrex::Real> aion_inv = {};

#endif
