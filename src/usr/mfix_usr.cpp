#include <mfix.H>

#include <MFIX_FLUID_Parms.H>

//
// Purpose: This routine is called before the time loop starts and is
//          user-definable.   This can be used for setting constants
//          and checking errors in data.
//
void
mfix::mfix_usr0 () const
{
}

void
mfix::mfix_usr1 (Real time) const
{
  // const int dim_bc = get_dim_bc();

  // for(unsigned i(1); i <= dim_bc; ++i)
  // {
  //   m_bc_u_g[i] = get_bc_u_g(i);
  //   m_bc_v_g[i] = get_bc_v_g(i);
  //   m_bc_w_g[i] = get_bc_w_g(i);

  //   m_bc_t_g[i] = get_bc_t_g(i);
  //
  //   Loop over species number
  //   m_bc_X_g[n][i] = get_bc_X_g(i,n);

  //   m_bc_ep_g[i] = get_bc_ep_g(i);
  // }
}

void
mfix::mfix_usr2 () const
{
}

//
// Purpose: This routine is called after the time loop ends and is user-definable.
//
void
mfix::mfix_usr3 () const
{
}
