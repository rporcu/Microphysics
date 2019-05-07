#include <AMReX_REAL.H>
#include <MFIX_DEM_Parms.H>

struct ComputeDragUser
{
    /**
     *  \brief Returns: the calculated drag coefficient.
     *
     *  Inputs:
     *    EPg     - gas volume fraction
     *    Mug     - gas laminar viscosity
     *    ROpg    - gas density * EP_g
     *    vrel    - magnitude of gas-solids relative velocity
     *    DPM     - particle diamater of solids phase M
     *    DPA     - average particle diameter
     *    PHIS    - solids volume fraction of solids phases
	 *    fvelx   - x component of the fluid velocity at the particle position
	 *    fvely   - y component of the fluid velocity at the particle position
	 *    fvelz   - z component of the fluid velocity at the particle position
	 *    i, j, k - particle cell indices
	 *    pid     - particle id number
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real operator() (amrex::Real EPg, amrex::Real Mug, amrex::Real ROPg, amrex::Real vrel, 
                            amrex::Real DPM, amrex::Real DPA, amrex::Real PHIS, 
							amrex::Real fvelx, amrex::Real fvely, amrex::Real fvelz,
							int i, int j, int k, int pid) const
    {
		amrex::Real ROg = ROPg / EPg;
		amrex::Real RE = (Mug > 0.0) ? DPM*vrel*ROg/Mug : DEMParams::large_number;

		amrex::Real Cd = 0.0;
		if (RE > DEMParams::eps) Cd = (24.0/RE)*(1.0 + 0.15*std::pow(RE, 0.687));

		return 0.75*(ROg*vrel/DPM)*Cd;
    }
};
