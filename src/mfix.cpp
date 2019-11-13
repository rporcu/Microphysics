#include <mfix_F.H>
#include <mfix_eb_F.H>
#include <mfix.H>
#include <param_mod_F.H>
#include <bc_mod_F.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

std::string      mfix::particle_init_type   = "AsciiFile";
std::string      mfix::load_balance_type    = "FixedSize";
std::string      mfix::knapsack_weight_type = "RunTimeCosts";
int              mfix::load_balance_fluid   = 1;
int              mfix::knapsack_nmax        = 128;
DragType         mfix::m_drag_type          = DragType::Invalid;
DepositionScheme mfix::m_deposition_scheme;
amrex::Real      mfix::tcoll_ratio          = 50.;
amrex::Real      mfix::m_deposition_diffusion_coeff = -1.0;
amrex::Real      mfix::m_deposition_scale_factor = 1.0;
amrex::Real      mfix::m_max_solids_volume_fraction = 0.90;

// Define unit vectors for easily convert indices
amrex::IntVect mfix::e_x(1,0,0);
amrex::IntVect mfix::e_y(0,1,0);
amrex::IntVect mfix::e_z(0,0,1);

int mfix::nlev  = 1;
int mfix::ntrac = 1;

EBSupport mfix::m_eb_support_level = EBSupport::full;

Real mfix::gravity[3] {0.0};
Real mfix::gp0[3]     {0.0};

mfix::~mfix () {};

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
mfix::mfix_usr1_cpp(Real time) const
{
  mfix_usr1(&time);

  const int dim_bc = get_dim_bc();

  for(unsigned i(1); i <= dim_bc; ++i)
  {
    m_bc_u_g[i] = get_bc_u_g(i);
    m_bc_v_g[i] = get_bc_v_g(i);
    m_bc_w_g[i] = get_bc_w_g(i);

    m_bc_t_g[i] = get_bc_t_g(i);

    m_bc_ep_g[i] = get_bc_ep_g(i);
  }
}

void
mfix::usr3()
{
    if (solve_fluid)
    {
       for (int lev = 0; lev < nlev; lev++)
       {
          Real dx = geom[lev].CellSize(0);
          Real dy = geom[lev].CellSize(1);
          Real dz = geom[lev].CellSize(2);

          // We deliberately don't tile this loop
          for (MFIter mfi(*p_g[lev]); mfi.isValid(); ++mfi)
          {
             mfix_usr3(BL_TO_FORTRAN_ANYD((*vel_g[lev])[mfi]),
                       BL_TO_FORTRAN_ANYD((  *p_g[lev])[mfi]),
                       &dx, &dy, &dz);
          }
       }
    }
}

void
mfix::mfix_set_bc_type(int lev)
{
    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);
    Real xlen = geom[lev].ProbHi(0) - geom[lev].ProbLo(0);
    Real ylen = geom[lev].ProbHi(1) - geom[lev].ProbLo(1);
    Real zlen = geom[lev].ProbHi(2) - geom[lev].ProbLo(2);
    Box domain(geom[lev].Domain());

    const int dim_bc = get_dim_bc();

    set_bc_type(bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                domain.loVect(),domain.hiVect(),
                &dx, &dy, &dz, &xlen, &ylen, &zlen, &nghost);

    for(unsigned i(1); i <= dim_bc; ++i)
    {
      m_bc_u_g[i] = get_bc_u_g(i);
      m_bc_v_g[i] = get_bc_v_g(i);
      m_bc_w_g[i] = get_bc_w_g(i);
    }
}

void mfix::mfix_set_bc_mod(const int* pID, const int* pType,
                           const amrex::Real* pLo, const amrex::Real* pHi,
                           amrex::Real* pLoc,
                           amrex::Real* pPg,
                           amrex::Real* pVel)
{
  const int dim_bc = get_dim_bc();

  set_bc_mod(pID, pType, pLo, pHi, pLoc, pPg, pVel);

  for(unsigned i(1); i <= dim_bc; ++i)
  {
    m_bc_u_g[i] = get_bc_u_g(i);
    m_bc_v_g[i] = get_bc_v_g(i);
    m_bc_w_g[i] = get_bc_w_g(i);

    m_bc_t_g[i] = get_bc_t_g(i);

    m_bc_ep_g[i] = get_bc_ep_g(i);

    m_bc_p_g[i] = get_bc_p_g(i);
  }
}

void mfix::mfix_set_bc_mod_add_mi(const int* pPlane,
                                  amrex::Real* xLo, amrex::Real* yLo, amrex::Real* zLo,
                                  amrex::Real* xHi, amrex::Real* yHi, amrex::Real* zHi,
                                  amrex::Real* pPg, amrex::Real* pVel)
{
  const int dim_bc = get_dim_bc();

  set_bc_mod_add_mi(pPlane, xLo, yLo, zLo, xHi, yHi, zHi, pPg, pVel);

  for(unsigned i(1); i <= dim_bc; ++i)
  {
    m_bc_u_g[i] = get_bc_u_g(i);
    m_bc_v_g[i] = get_bc_v_g(i);
    m_bc_w_g[i] = get_bc_w_g(i);

    m_bc_t_g[i] = get_bc_t_g(i);

    m_bc_ep_g[i] = get_bc_ep_g(i);

    m_bc_p_g[i] = get_bc_p_g(i);
  }
}

void
mfix::avgDown (int crse_lev, const MultiFab& S_fine, MultiFab& S_crse)
{
    BL_PROFILE("mfix::avgDown()");

    amrex::EB_average_down(S_fine, S_crse, 0, S_fine.nComp(), refRatio(crse_lev));
}
