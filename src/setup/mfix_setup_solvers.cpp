#include <mfix.H>

#include <hydro_NodalProjector.H>
#include <mfix_diffusion_op.H>
#include <mfix_bc.H>

void
mfix::mfix_init_solvers ()
{
    BL_PROFILE("mfix::mfix_init_solvers");

    diffusion_op = std::make_unique<DiffusionOp>(this, amrex::GetVecOfConstPtrs(ebfactory), m_embedded_boundaries, fluid,
                                       m_boundary_conditions.diff_vel_lobc(),         m_boundary_conditions.diff_vel_hibc(),
                                       m_boundary_conditions.diff_scal_lobc(),        m_boundary_conditions.diff_scal_hibc(),
                                       m_boundary_conditions.diff_temperature_lobc(), m_boundary_conditions.diff_temperature_hibc(),
                                       m_boundary_conditions.diff_species_lobc(),     m_boundary_conditions.diff_species_hibc(),
                                       nghost_state());

    macproj = std::make_unique<Hydro::MacProjector>(Geom(0,finest_level),
                                   MLMG::Location::FaceCentroid,  // Location of mac_vec
                                   MLMG::Location::FaceCentroid,  // Location of beta
                                   MLMG::Location::CellCenter,    // Location of solution variable phi
                                   MLMG::Location::CellCentroid);// Location of MAC RHS
}

void
mfix::mfix_setup_solvers ()
{
    BL_PROFILE("mfix::mfix_setup_solvers");

    diffusion_op->setup(this, amrex::GetVecOfConstPtrs(ebfactory));
}
