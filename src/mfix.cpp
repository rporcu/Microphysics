#include <mfix.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>

std::string      mfix::particle_init_type   = "AsciiFile";
std::string      mfix::load_balance_type    = "KnapSack";
std::string      mfix::knapsack_weight_type = "RunTimeCosts";
int              mfix::load_balance_fluid   = 1;
int              mfix::knapsack_nmax        = 128;
int              mfix::m_drag_type          = DragType::Invalid;
int              mfix::m_convection_type    = ConvectionType::Invalid;
DepositionScheme mfix::m_deposition_scheme;
int              mfix::m_reaction_rates_type = ReactionRatesType::RRatesUser;
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
  if (DEM::solve)
    delete pc;

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
  }

  // used if load_balance_type == "KnapSack"
  for (int lev = 0; lev < particle_cost.size(); ++lev)
    delete particle_cost[lev];

  for (int lev = 0; lev < fluid_cost.size(); ++lev)
    delete fluid_cost[lev];

  if (REACTIONS::solve) {
    for (int n(0); n < REACTIONS::nreactions; n++)
      delete m_chemical_reactions[n];
  }
};

// Constructor
mfix::mfix ()
  : m_bc_u_g(50, 0)
  , m_bc_v_g(50, 0)
  , m_bc_w_g(50, 0)
  , m_bc_t_g(50, 0)
  , m_bc_ep_g(50, 0)
  , m_bc_p_g(50, 0)
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
    // This needs to be one bigger than the highest index scalar in mfix_set_scalar_bcs
    bcs_s.resize(6); // density, tracer, ep_g, mu_g, T_g, h_g --> TODO cp_g, k_g
    bcs_X.resize(0); // X_gk, D_gk. TODO this has to be resized on the basis of
                     // FLUID::nspecies. So we do it after parameter parsing
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

    m_vel_g_bc_types["Dirichlet"] = {bc_list.get_minf()};
    m_vel_g_bc_types["Neumann"] = {bc_list.get_pinf(), bc_list.get_pout()};

    m_ro_g_bc_types["Dirichlet"] = {bc_list.get_minf()};
    m_ro_g_bc_types["Neumann"] = {bc_list.get_pinf(), bc_list.get_pout()};

    m_T_g_bc_types["Dirichlet"] = {bc_list.get_minf(), bc_list.get_pinf()};
    m_T_g_bc_types["Neumann"] = {bc_list.get_pout()};

    m_trac_g_bc_types["Dirichlet"] = {bc_list.get_minf()};
    m_trac_g_bc_types["Neumann"] = {bc_list.get_pinf(), bc_list.get_pout()};

    m_X_gk_bc_types["Dirichlet"] = {bc_list.get_minf(), bc_list.get_pinf()};
    m_X_gk_bc_types["Neumann"] = {bc_list.get_pout()};

    Gpu::synchronize();
}

void
mfix::avgDown (int crse_lev, const MultiFab& S_fine, MultiFab& S_crse)
{
    BL_PROFILE("mfix::avgDown()");

    amrex::EB_average_down(S_fine, S_crse, 0, S_fine.nComp(), refRatio(crse_lev));
}

Vector< MultiFab* > mfix::get_ep_g () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->ep_g);
  }
  return r;
}

Vector< MultiFab* > mfix::get_ro_g () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->ro_g);
  }
  return r;
}

Vector< MultiFab* > mfix::get_ro_g_old () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->ro_go);
  }
  return r;
}

Vector< MultiFab* > mfix::get_MW_g () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->MW_g);
  }
  return r;
}

Vector< MultiFab* > mfix::get_trac () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->trac);
  }
  return r;
}

Vector< MultiFab* > mfix::get_trac_old () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->trac_o);
  }
  return r;
}

Vector< MultiFab* > mfix::get_vel_g () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->vel_g);
  }
  return r;
}

Vector< MultiFab* > mfix::get_vel_g_old () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->vel_go);
  }
  return r;
}

Vector< MultiFab* > mfix::get_mu_g () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->mu_g);
  }
  return r;
}

Vector< MultiFab* > mfix::get_T_g () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->T_g);
  }
  return r;
}

Vector< MultiFab* > mfix::get_T_g_old () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->T_go);
  }
  return r;
}

Vector< MultiFab* > mfix::get_cp_g () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->cp_g);
  }
  return r;
}

Vector< MultiFab* > mfix::get_k_g () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->k_g);
  }
  return r;
}

Vector< MultiFab* > mfix::get_h_g () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->h_g);
  }
  return r;
}

Vector< MultiFab* > mfix::get_h_g_old () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->h_go);
  }
  return r;
}

Vector< MultiFab* > mfix::get_T_g_on_eb () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->T_g_on_eb);
  }
  return r;
}

Vector< MultiFab* > mfix::get_k_g_on_eb () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->k_g_on_eb);
  }
  return r;
}

Vector< MultiFab* > mfix::get_X_gk () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->X_gk);
  }
  return r;
}

Vector< MultiFab* > mfix::get_X_gk_old () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->X_gko);
  }
  return r;
}

Vector< MultiFab* > mfix::get_D_gk () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->D_gk);
  }
  return r;
}

Vector< MultiFab* > mfix::get_cp_gk () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->cp_gk);
  }
  return r;
}

Vector< MultiFab* > mfix::get_h_gk () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->h_gk);
  }
  return r;
}

Vector< MultiFab* > mfix::get_txfr () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->txfr);
  }
  return r;
}

Vector< MultiFab* > mfix::get_ro_gk_txfr () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->ro_gk_txfr);
  }
  return r;
}

Vector< MultiFab* > mfix::get_xslopes_u () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->xslopes_u);
  }
  return r;
}

Vector< MultiFab* > mfix::get_yslopes_u () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->yslopes_u);
  }
  return r;
}

Vector< MultiFab* > mfix::get_zslopes_u () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->zslopes_u);
  }
  return r;
}

Vector< MultiFab* > mfix::get_xslopes_s () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->xslopes_s);
  }
  return r;
}

Vector< MultiFab* > mfix::get_yslopes_s () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->yslopes_s);
  }
  return r;
}

Vector< MultiFab* > mfix::get_zslopes_s () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->zslopes_s);
  }
  return r;
}

Vector< MultiFab* > mfix::get_xslopes_X_gk () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->xslopes_X);
  }
  return r;
}

Vector< MultiFab* > mfix::get_yslopes_X_gk () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->yslopes_X);
  }
  return r;
}

Vector< MultiFab* > mfix::get_zslopes_X_gk () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->zslopes_X);
  }
  return r;
}

Vector< MultiFab* > mfix::get_diveu () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->diveu);
  }
  return r;
}

Vector< MultiFab* > mfix::get_mac_phi () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->mac_phi);
  }
  return r;
}

Vector< MultiFab* > mfix::get_mac_rhs () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->mac_rhs);
  }
  return r;
}
