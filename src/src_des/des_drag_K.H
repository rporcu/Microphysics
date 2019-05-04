#include <AMReX_REAL.H>
#include <MFIX_DEM_Parms.H>
/**
   
   Subroutine: DRAG_BVK2                                               
   Purpose: Calculate the gas-solids drag coefficient                  
   
   Literature/Document References:                                     
       Beetstra, van der Hoef, Kuipers, Chem. Eng. Science 62           
	   (Jan 2007)                                                       
   
   Tang, Peters, Kuipers, Kriebitzsch, van der Hoef, AIChEJ,        
       61(2) (Feb 2015)                                                 
                                                                      
 */
struct ComputeDragBVK2
{
    /**
     *  \brief Function for computing drag according to the BVK drag model.
     *  
     *  Returns: the calculated drag coefficient.
     *
     *  Inputs:
     *    EPg  - gas volume fraction
     *    Mug  - gas laminar viscosity
     *    ROpg - gas density * EP_g
     *    vrel - magnitude of gas-solids relative velocity
     *    DPM  - particle diamater of solids phase M
     *    DPA  - average particle diameter
     *    PHIS - solids volume fraction of solids phases
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real operator() (amrex::Real EPg, amrex::Real Mug, amrex::Real ROPg, amrex::Real vrel,
                            amrex::Real DPM, amrex::Real DPA, amrex::Real PHIS) const
    {
        amrex::Real RE = (Mug > 0.0) ? DPA*vrel*ROPg/Mug : DEMParams::large_number;
    
        if (RE > DEMParams::eps)
        {
			amrex::Real oEPgfour = 1.0/EPg/EPg/EPg/EPg;

			// eq(9) BVK J. fluid. Mech. 528, 2005
			// (this F_Stokes is /= of Koch_Hill by a factor of ep_g)
            amrex::Real F_Stokes = 18.0*Mug*EPg/DPM/DPM;

            amrex::Real F = 10.0*PHIS/EPg/EPg + EPg*EPg*(1.0 + 1.5*sqrt(PHIS));

			F += RE*(0.11*PHIS*(1.0+PHIS) - 4.56e-3*oEPgfour + 
					 std::pow(RE, -0.343)*(0.169*EPg + 6.44e-2*oEPgfour));

            // F += 0.413*RE/(24.0*EPg*EPg) *
            //     (1.0/EPg + 3.0*EPg*PHIS + 8.4/std::pow(RE, 0.343)) /
            //     (1.0 + std::pow(10.0, 3.0*PHIS)/std::pow(RE, 0.5 + 2.0*PHIS));
            
            return F*F_Stokes;
        }
        else
        {
            return 0.0;
        }        
    }
};

struct ComputeDragKochHill
{

    /**
     *  \brief Function for computing drag according to the Koch-Hill drag model.
     *  
     *  Returns: the calculated drag coefficient.
     *
     *  Inputs:
     *    EPg  - gas volume fraction
     *    Mug  - gas laminar viscosity
     *    ROpg - gas density * EP_g
     *    vrel - magnitude of gas-solids relative velocity
     *    DPM  - particle diamater of solids phase M
     *    DPA  - average particle diameter
     *    PHIS - solids volume fraction of solids phases
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real operator() (amrex::Real EPg, amrex::Real Mug, amrex::Real ROPg, amrex::Real vrel, 
                            amrex::Real DPM, amrex::Real DPA, amrex::Real PHIS) const
    {
        amrex::Abort("Not implemented");

        return 0.0;
    }
};
