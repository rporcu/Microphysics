#include <mfix.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

#include <mfix_fluid.H>
#include <mfix_species.H>

std::string      mfix::load_balance_type    = "KnapSack";
std::string      mfix::knapsack_weight_type = "RunTimeCosts";
int              mfix::load_balance_fluid   = -1;
int              mfix::knapsack_nmax        = 128;
int              mfix::greedy_dir           = 0;
bool             mfix::greedy_3d            = false;
int              mfix::greedy_min_grid_size = 2;
int              mfix::m_run_type           = RunType::Standard;
int              mfix::m_drag_type          = DragType::Invalid;
amrex::Real      mfix::m_SyamOBrien_coeff_c1 = -1;
amrex::Real      mfix::m_SyamOBrien_coeff_d1 = -1;
int              mfix::m_convection_type    = ConvectionType::Invalid;
int              mfix::m_advection_type     = AdvectionType::Invalid;
DepositionScheme mfix::m_deposition_scheme;
int              mfix::m_reaction_rates_type = ReactionRatesType::RRatesUser;
amrex::Real      mfix::m_deposition_diffusion_coeff = -1.0;
amrex::Real      mfix::m_deposition_scale_factor = 1.0;
amrex::Real      mfix::m_max_solids_volume_fraction = 0.64356;

int mfix::nlev  = 1;
int mfix::ntrac = 1;

EBSupport mfix::m_eb_support_level = EBSupport::full;

RealVect mfix::gravity {0.};
RealVect mfix::gp0     {0.};

// Destructor
mfix::~mfix ()
{
  if ((m_dem.solve() || m_pic.solve()) && pc != nullptr)
    delete pc;

  for (int lev(0); lev < nlev; ++lev)
  {
    // Face-based coefficients b in MAC projection and implicit diffusion solve
    for (int dir(0); dir < 3; ++dir)
      if (bcoeff[lev][dir] != nullptr)
        delete bcoeff[lev][dir];
  }

  // used if load_balance_type == "KnapSack"
  for (int lev = 0; lev < particle_cost.size(); ++lev)
    delete particle_cost[lev];

  for (int lev = 0; lev < particle_proc.size(); ++lev)
    delete particle_proc[lev];

  for (int lev = 0; lev < fluid_cost.size(); ++lev)
    delete fluid_cost[lev];

  for (int lev = 0; lev < fluid_proc.size(); ++lev)
    delete fluid_proc[lev];

  delete m_rw;
  delete m_solids_volume_deposition;
  delete m_interphase_txfr_deposition;
}

// Constructor
mfix::mfix ()
  : m_boundary_conditions(geom, bc_list, m_embedded_boundaries)
  , bc_list(maxLevel() + 1)
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

    // Generic first-order extrapolation
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

    {
      ParmParse pp("mfix");

      if (pp.contains("restart"))
        m_run_type = RunType::Restart;
      else if (ParmParse("pic2dem").contains("convert"))
        m_run_type = RunType::PIC2DEM;
      else
        m_run_type = RunType::Standard;
    }

    m_rw = new MFIXReadWrite(nlev, grids, geom, pc, fluid, m_leveldata,
                             ebfactory, dmap, ooo_debug, level_sets,
                             levelset_refinement, levelset_pad,
                             levelset_eb_refinement, levelset_eb_pad, solids,
                             m_dem, m_pic, reactions, particle_cost,
                             particle_proc, fluid_proc, refRatio(),
                             bc_list, particle_ebfactory,
                             regions);

    m_solids_volume_deposition = new MFIXSolidsVolume;
    m_interphase_txfr_deposition = new MFIXInterphaseTxfr;
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

Vector< MultiFab* > mfix::get_gp () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->gp);
  }
  return r;
}

Vector< MultiFab* > mfix::get_p_g () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->p_g);
  }
  return r;
}

Vector< MultiFab* > mfix::get_p_g_old () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->p_go);
  }
  return r;
}

Vector< Real* > mfix::get_thermodynamic_p_g () noexcept
{
  Vector<Real*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->thermodynamic_p_g);
  }
  return r;
}

Vector< Real* > mfix::get_thermodynamic_p_g_old () noexcept
{
  Vector<Real*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->thermodynamic_p_go);
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

Vector< MultiFab* > mfix::get_txfr () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->txfr);
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

Vector< MultiFab* > mfix::get_divtau () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->divtau_o);
  }
  return r;
}








//////////////////////////////////////









/////////////////////////



Vector< MultiFab const*> mfix::get_ep_g_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->ep_g);
  }
  return r;
}

Vector< MultiFab const*> mfix::get_ro_g_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->ro_g);
  }
  return r;
}

Vector< MultiFab const*> mfix::get_ro_g_old_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->ro_go);
  }
  return r;
}

Vector< MultiFab const*> mfix::get_trac_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->trac);
  }
  return r;
}

Vector< MultiFab const*> mfix::get_trac_old_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->trac_o);
  }
  return r;
}

Vector< MultiFab const*> mfix::get_vel_g_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->vel_g);
  }
  return r;
}

Vector< MultiFab const*> mfix::get_vel_g_old_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->vel_go);
  }
  return r;
}

Vector< MultiFab const*> mfix::get_p_g_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->p_g);
  }
  return r;
}

Vector< MultiFab const*> mfix::get_p_g_old_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->p_go);
  }
  return r;
}

Vector< MultiFab const*> mfix::get_T_g_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->T_g);
  }
  return r;
}

Vector< MultiFab const*> mfix::get_T_g_old_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->T_go);
  }
  return r;
}

Vector< MultiFab const*> mfix::get_h_g_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->h_g);
  }
  return r;
}

Vector< MultiFab const*> mfix::get_h_g_old_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->h_go);
  }
  return r;
}

Vector< Real const*> mfix::get_thermodynamic_p_g_const () const noexcept
{
  Vector<Real const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->thermodynamic_p_g);
  }
  return r;
}

Vector< Real const*> mfix::get_thermodynamic_p_g_old_const () const noexcept
{
  Vector<Real const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->thermodynamic_p_go);
  }
  return r;
}

Vector< MultiFab const*> mfix::get_T_g_on_eb_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->T_g_on_eb);
  }
  return r;
}

Vector< MultiFab const*> mfix::get_X_gk_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->X_gk);
  }
  return r;
}

Vector< MultiFab const*> mfix::get_X_gk_old_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->X_gko);
  }
  return r;
}

Vector< MultiFab const*> mfix::get_txfr_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->txfr);
  }
  return r;
}

Vector< MultiFab const*> mfix::get_diveu_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->diveu);
  }
  return r;
}

Vector< MultiFab const*> mfix::get_mac_phi_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->mac_phi);
  }
  return r;
}

Vector< MultiFab const*> mfix::get_divtau_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->divtau_o);
  }
  return r;
}

void mfix::build_eb_levels_from_chkpt_file () {

   EB2::BuildFromChkptFile(m_rw->geom_chk_file, geom[nlev-1], nlev-1, 100);

   const EB2::IndexSpace& ebis = EB2::IndexSpace::top();
   for (int lev = 0; lev < nlev; lev ++)
   {
      eb_levels[lev] = &(ebis.getLevel(geom[lev]));
      particle_eb_levels[lev] = eb_levels[lev];
   }

   if (nlev == 1)
   {
      if (levelset_refinement == 1) {
           eb_levels[1] = eb_levels[0];
           particle_eb_levels[1] = eb_levels[1];
      }
      else {
         Geometry geom_ls = amrex::refine(geom[0],levelset_refinement);
         EB2::BuildFromChkptFile(m_rw->geom_levelset_chk_file, geom_ls, 0, 100);
         eb_levels[1] = &(EB2::IndexSpace::top().getLevel(geom_ls));
         particle_eb_levels[1] = eb_levels[1];
      }
   }
}
