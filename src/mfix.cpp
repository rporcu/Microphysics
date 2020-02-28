#include <mfix_F.H>
#include <mfix_eb_F.H>
#include <mfix.H>
#include <param_mod_F.H>


#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

#include <MFIX_FLUID_Parms.H>

std::string      mfix::particle_init_type   = "AsciiFile";
std::string      mfix::load_balance_type    = "KnapSack";
std::string      mfix::knapsack_weight_type = "RunTimeCosts";
int              mfix::load_balance_fluid   = 1;
int              mfix::knapsack_nmax        = 128;
DragType         mfix::m_drag_type          = DragType::Invalid;
DepositionScheme mfix::m_deposition_scheme;
amrex::Real      mfix::tcoll_ratio          = 50.;
amrex::Real      mfix::m_deposition_diffusion_coeff = -1.0;
amrex::Real      mfix::m_deposition_scale_factor = 1.0;
amrex::Real      mfix::m_max_solids_volume_fraction = 0.64356;

int mfix::nlev  = 1;
int mfix::ntrac = 1;

int  mfix::plot_int        = -1;
Real mfix::plot_per_approx = -1.;
Real mfix::plot_per_exact  = -1.;

EBSupport mfix::m_eb_support_level = EBSupport::full;

RealVect mfix::gravity {0.};
RealVect mfix::gp0     {0.};

// Destructor
mfix::~mfix ()
{
  for (int lev(0); lev < nlev; ++lev)
  {
    // Face-based coefficients b in MAC projection and implicit diffusion solve
    delete bcoeff[lev][0];
    delete bcoeff[lev][1];
    delete bcoeff[lev][2];

    // Boundary conditions types
    delete bc_ilo[lev];
    delete bc_ihi[lev];
    delete bc_jlo[lev];
    delete bc_jhi[lev];
    delete bc_klo[lev];
    delete bc_khi[lev];

    // Level-Set Data
    delete level_sets[lev];

    // RHS for MAC solve
    delete mac_rhs[lev];
  }

  //! EB factory that lives on the fluid grids
  for (int lev(0); lev < ebfactory.size(); lev++)
    delete ebfactory[lev];

  //! EB factory that lives on the particle grids
  for (int lev(0); lev < particle_ebfactory.size(); ++lev)
    delete particle_ebfactory[lev];
};

// Constructor
mfix::mfix ()
  : m_bc_u_g(get_dim_bc()+1, 0)
  , m_bc_v_g(get_dim_bc()+1, 0)
  , m_bc_w_g(get_dim_bc()+1, 0)
  , m_bc_t_g(get_dim_bc()+1, 0)
  , m_bc_ep_g(get_dim_bc()+1, 0)
  , m_bc_p_g(get_dim_bc()+1, 0)
{
    // NOTE: Geometry on all levels has just been defined in the AmrCore
    // constructor. No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    /****************************************************************************
     *                                                                          *
     * Set max number of levels (nlevs)                                         *
     *                                                                          *
     ***************************************************************************/

    nlev = maxLevel() + 1;
    amrex::Print() << "Number of levels: " << nlev << std::endl;


    /****************************************************************************
     *                                                                          *
     * Initialize time steps                                                    *
     *                                                                          *
     ***************************************************************************/

    t_old.resize(nlev,-1.e100);
    t_new.resize(nlev,0.0);


    /****************************************************************************
     *                                                                          *
     * Initialize boundary conditions (used by fill-patch)                      *
     *                                                                          *
     ***************************************************************************/

    bcs_u.resize(3); // one for each velocity component
    bcs_s.resize(2); // density and tracer
    bcs_f.resize(1); // just one

    //___________________________________________________________________________
    // Boundary conditions used for level-sets

    bcs_ls.resize(1);

    // walls (Neumann)
    int bc_lo[] = {FOEXTRAP, FOEXTRAP, FOEXTRAP};
    int bc_hi[] = {FOEXTRAP, FOEXTRAP, FOEXTRAP};

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        // lo-side BCs
        if (bc_lo[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_lo[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_lo[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs_ls[0].setLo(idim, bc_lo[idim]);
        }
        else {
            amrex::Abort("Invalid level-set bc_lo");
        }

        // hi-side BCSs
        if (bc_hi[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_hi[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_hi[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs_ls[0].setHi(idim, bc_hi[idim]);
        }
        else {
            amrex::Abort("Invalid level-set bc_hi");
        }
    }
}

void
mfix::mfix_usr1_cpp (Real time) const
{
  mfix_usr1(&time);

  const int dim_bc = get_dim_bc();

  // for(unsigned i(1); i <= dim_bc; ++i)
  // {
  //   m_bc_u_g[i] = get_bc_u_g(i);
  //   m_bc_v_g[i] = get_bc_v_g(i);
  //   m_bc_w_g[i] = get_bc_w_g(i);

  //   m_bc_t_g[i] = get_bc_t_g(i);

  //   m_bc_ep_g[i] = get_bc_ep_g(i);
  // }
}

void
mfix::usr3 ()
{
    if (FLUID::solve)
    {
       for (int lev = 0; lev < nlev; lev++)
       {
          Real dx = geom[lev].CellSize(0);
          Real dy = geom[lev].CellSize(1);
          Real dz = geom[lev].CellSize(2);

          // We deliberately don't tile this loop
          for (MFIter mfi(*(m_leveldata[lev]->p_g), false); mfi.isValid(); ++mfi)
          {
             mfix_usr3(BL_TO_FORTRAN_ANYD((*m_leveldata[lev]->vel_g)[mfi]),
                       BL_TO_FORTRAN_ANYD((*m_leveldata[lev]->p_g)[mfi]),
                       &dx, &dy, &dz);
          }
       }
    }
}


void
mfix::avgDown (int crse_lev, const MultiFab& S_fine, MultiFab& S_crse)
{
    BL_PROFILE("mfix::avgDown()");

    amrex::EB_average_down(S_fine, S_crse, 0, S_fine.nComp(), refRatio(crse_lev));
}
