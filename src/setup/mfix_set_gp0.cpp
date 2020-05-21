#include <mfix.H>
#include <climits>

#include <MFIX_FLUID_Parms.H>
#include <MFIX_BC_Parms.H>
#include <MFIX_IC_Parms.H>

#include <string>


void
mfix::set_gp0 (const int lev,
               const Box& domain)
{
  const Real tolerance = std::numeric_limits<Real>::epsilon();

  const IntVect dom_lo(domain.loVect());
  const IntVect dom_hi(domain.hiVect());

  Array4<const int> const& bct_ilo = bc_ilo[lev]->array();
  Array4<const int> const& bct_ihi = bc_ihi[lev]->array();
  Array4<const int> const& bct_jlo = bc_jlo[lev]->array();
  Array4<const int> const& bct_jhi = bc_jhi[lev]->array();
  Array4<const int> const& bct_klo = bc_klo[lev]->array();
  Array4<const int> const& bct_khi = bc_khi[lev]->array();

  int delp_dir_loc(BC::delp_dir);

  Real delp_x = BC::delp[0];
  Real delp_y = BC::delp[1];
  Real delp_z = BC::delp[2];

  // ---------------------------------------------------------------->>>
  //     If the bc's are pressure inflow/outflow then be sure to
  //     capture that in p0 and gp0.
  // ---------------------------------------------------------------->>>

  if ( ( (bct_ilo(dom_lo[0]-1,dom_lo[1],dom_lo[2],0) == bc_list.get_pinf()) and
         (bct_ihi(dom_hi[0]+1,dom_lo[1],dom_lo[2],0) == bc_list.get_pout()) )
       or
       ( (bct_ihi(dom_hi[0]+1,dom_lo[1],dom_lo[2],0) == bc_list.get_pinf()) and
         (bct_ilo(dom_lo[0]-1,dom_lo[1],dom_lo[2],0) == bc_list.get_pout()) ) )
    {

      const int bcv_lo = bct_ilo(dom_lo[0]-1,dom_lo[1],dom_lo[2],1);
      const int bcv_hi = bct_ihi(dom_hi[0]+1,dom_lo[1],dom_lo[2],1);

      const Real  p_lo = m_bc_p_g[bcv_lo];
      const Real  p_hi = m_bc_p_g[bcv_hi];

      delp_dir_loc = 0;
      delp_x = p_lo - p_hi;
    }


  else if( ( (bct_jlo(dom_lo[0],dom_lo[1]-1,dom_lo[2],0) == bc_list.get_pinf()) and
             (bct_jhi(dom_lo[0],dom_hi[1]+1,dom_lo[2],0) == bc_list.get_pout()) )
           or
           ( (bct_jhi(dom_lo[0],dom_hi[1]+1,dom_lo[2],0) == bc_list.get_pinf()) and
             (bct_jlo(dom_lo[0],dom_lo[1]-1,dom_lo[2],0) == bc_list.get_pout()) ) )
    {

      const int bcv_lo = bct_jlo(dom_lo[0],dom_lo[1]-1,dom_lo[2],1);
      const int bcv_hi = bct_jhi(dom_lo[0],dom_hi[1]+1,dom_lo[2],1);

      const Real  p_lo = m_bc_p_g[bcv_lo];
      const Real  p_hi = m_bc_p_g[bcv_hi];

      delp_dir_loc = 1;
      delp_y = p_lo - p_hi;
    }

  else if( ( (bct_klo(dom_lo[0],dom_lo[1],dom_lo[2]-1,0) == bc_list.get_pinf()) and
             (bct_khi(dom_lo[0],dom_lo[1],dom_hi[2]+1,0) == bc_list.get_pout()) )
           or
           ( (bct_khi(dom_lo[0],dom_lo[1],dom_hi[2]+1,0) == bc_list.get_pinf()) and
             (bct_klo(dom_lo[0],dom_lo[1],dom_lo[2]-1,0) == bc_list.get_pout()) ) )
    {

      const int bcv_lo = bct_klo(dom_lo[0],dom_lo[1],dom_lo[2]-1,1);
      const Real p_lo  = m_bc_p_g[bcv_lo];

      const int bcv_hi = bct_khi(dom_lo[0],dom_lo[1],dom_hi[2]+1,1);
      const Real p_hi  = m_bc_p_g[bcv_hi];

      delp_dir_loc = 2;
      delp_z = p_lo - p_hi;
  }

  // ---------------------------------------------------------------->>>

  //  Make sure that ic_p_g is set if using delp pressure conditions
  for(int icv(0); icv < IC::ic.size(); ++icv)
  {

    if((delp_dir_loc >= 0) and (delp_dir_loc == BC::delp_dir)) {

      if (not IC::ic[icv].fluid.pressure_defined) {
        std::cout << "MUST DEFINE ic_p_g if using the DELP pressure condition" << std::endl;
        exit(0);
      }

    }
    else if((delp_dir_loc >= 0) and (delp_dir_loc != BC::delp_dir)) {
      if( IC::ic[icv].fluid.pressure_defined ) {
        std::cout << "MUST not define ic_p_g if setting p_inflowandp_outflow" << std::endl;
        exit(0);
      }

    } else {

      const amrex::Real gravity_square_module =
        gravity[0]*gravity[0] + gravity[1]*gravity[1] + gravity[2]*gravity[2];

      if( not IC::ic[icv].fluid.pressure_defined or gravity_square_module > tolerance)
        {

          // HACK: This should probably take into consideration
          // variable fluid density.
          const Real ro_g0 = FLUID::ro_g0;

          for (int dim=0; dim<3; dim++){
            if (gravity[dim]*gravity[dim] > tolerance ) {
              gp0[dim] = ro_g0 * gravity[dim];
            }
          }

        return;
      }
    }
  }

  // ---------------------------------------------------------------->>>

  // Here the pressure in each cell is determined from a specified pressure
  // drop across the domain le>. This section requires that the pressure
  // is already defined in all initial condition regions (otherwise this
  // section would be skipped)

  Real xlen = geom[0].ProbHi(0) - geom[0].ProbLo(0);
  Real ylen = geom[0].ProbHi(1) - geom[0].ProbLo(1);
  Real zlen = geom[0].ProbHi(2) - geom[0].ProbLo(2);

  if (std::abs(delp_x) > tolerance){
    gp0[0] = -delp_x / xlen;
  }

  if (std::abs(delp_y) > tolerance){
    gp0[1] = -delp_y / ylen;
  }

  if (std::abs(delp_z) > tolerance){
    gp0[2] = -delp_z / zlen;
  }


  return;
}


