#ifndef MFIX_RESTARTER_H_
#define MFIX_RESTARTER_H_

#include <mfix.H>
//#include <mfix_solids.H>
#include <mfix_des_parts_gen_K.H>


struct Transfer {

  Transfer (const MFIXSolidsPhase& solids)
    : idx_eps(0)
    , idx_vel(1)
    , idx_temp(4)
    , idx_species(idx_temp + solids.solve_enthalpy())
    , count(idx_species + solids.nspecies()*solids.solve_species())
  {}

  const int idx_eps;
  const int idx_vel;
  const int idx_temp;
  const int idx_species;
  const int count;
};


class MFIXPICDeposition : public MFIXDepositionOp
{
  public:
    MFIXPICDeposition () = default;

    virtual void deposit (int lev,
                          const amrex::Geometry& geom,
                          MFIXParticleContainer* pc,
                          const amrex::MultiFab* volfrac,
                          const amrex::FabArray<amrex::EBCellFlagFab>* flags,
                          amrex::MultiFab* txfr_mf,
                          amrex::MultiFab* eps_mf);

    template <typename F>
    void deposit (F WeightFunc,
                  int lev,
                  const amrex::Geometry& geom,
                  MFIXParticleContainer* pc,
                  const amrex::MultiFab* volfrac,
                  const amrex::FabArray<amrex::EBCellFlagFab>* flags,
                  amrex::MultiFab* txfr_mf,
                  amrex::MultiFab* eps_mf);

    void set_eps_size (const int nlev) { m_eps.resize(nlev, nullptr); }

    Vector< MultiFab* > const& get_eps () const { return m_eps; }

    void clear_eps ()
    {
      for (int lev(0); lev < m_eps.size(); ++lev)
        if (m_eps[lev] != nullptr)
          delete m_eps[lev];
    }

  private:
    Vector< MultiFab* > m_eps;
};


class MFIXRestarter
{
  public:
    MFIXRestarter (const int nlev);

    ~MFIXRestarter ();

    void allocate_coarse_arrays (const mfix* mfix_coarse);

    void allocate_fine_arrays (const mfix* mfix_fine);

    void change_inputs_table () const;

    void set_fine_objects (mfix* mfix_fine,
                           const mfix* mfix_coarse) const;

    void calc_txfr (const mfix* mfix_coarse,
                    const amrex::Vector<amrex::MultiFab*>& avgdPIC,
                    const amrex::Real time);

    void txfr_fluid_data (const mfix* mfix_coarse,
                          mfix* mfix_fine) const;

    void get_particles_radius_and_density (const mfix* mfix_ptr);

    void generate_particles (const mfix* mfix_coarse,
                             mfix* mfix_fine) const;

    void generate_particles (const long particles_count,
                             const amrex::Box& bx,
                             MFIXParticleContainer::ParticleTileType& particles,
                             const Hex_ClosePack* hcp_vector_ptr,
                             const Long* gen_number_ptr,
                             const Long* gen_indexes_ptr,
                             const int id,
                             const int cpu) const;

    void init_particles (const mfix* mfix_coarse,
                         MFIXParticleContainer* pc_fine) const;

    void init_particles_data (mfix* mfix_fine) const;

    MFIXPICDeposition* m_PIC_deposition;

    int m_refinement_ratio;

    amrex::Real m_eps_tolerance;
    amrex::Real m_eps_overflow;

    amrex::Real m_inputs_pdiameter;
    amrex::Real m_inputs_pdensity;

    int nlev;

    amrex::Vector<amrex::MultiFab*> avgdPIC_coarse;
    amrex::Vector<amrex::MultiFab*> avgdPIC_fine;
};


class TxfrAuxiliary
{
  public:
    AMREX_GPU_HOST_DEVICE
    TxfrAuxiliary () = default;

    AMREX_GPU_HOST_DEVICE
    amrex::RealVect get_coordinates_cc (const amrex::GeometryData& geom,
                                        const amrex::IntVect& indexes) const
    {
      const amrex::Real* plo = geom.ProbLo();
      const amrex::Real* dx  = geom.CellSize();

      amrex::RealVect result(0.);

      for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
        result[dir] = (0.5 + amrex::Real(indexes[dir])) * dx[dir] + plo[dir];
      }

      return result;
    }

    AMREX_GPU_HOST_DEVICE
    amrex::IntVect get_indexes_cc (const amrex::GeometryData& geom,
                                   const amrex::RealVect& coords) const
    {
      const amrex::Real* plo = geom.ProbLo();
      const amrex::Real* dx  = geom.CellSize();

      amrex::IntVect result(0);

      for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
        result[dir] = amrex::Math::floor<int>((coords[dir]-plo[dir]) / dx[dir]);
      }

      return result;
    }

    AMREX_GPU_HOST_DEVICE
    void get_weights_cc (amrex::GpuArray<amrex::GpuArray<amrex::GpuArray<amrex::Real,3>,3>,3>& weights,
                         const amrex::Array4<const amrex::Real>& volfrac_coarse_arr,
                         const amrex::Array4<const amrex::Real>& epg_coarse_arr,
                         const amrex::Box& domain,
                         const amrex::IntVect& ijk_coarse,
                         const amrex::Array4<const amrex::EBCellFlag>& flags_coarse_arr,
                         const amrex::IntVect& ijk_fine,
                         const amrex::Array4<const amrex::EBCellFlag>& flags_fine_arr) const
    {
      if ((!flags_fine_arr(ijk_fine).isCovered()) && (flags_coarse_arr(ijk_coarse).isCovered())) {

        amrex::Real sum(0.);

        for (int ii(-1); ii <= 1; ++ii)
        for (int jj(-1); jj <= 1; ++jj)
        for (int kk(-1); kk <= 1; ++kk) {
          amrex::IntVect shift(ii,jj,kk);
          if (!flags_coarse_arr(ijk_coarse+shift).isCovered() &&
              domain.contains(ijk_coarse+shift)) {

            const amrex::Real vfrac = volfrac_coarse_arr(ijk_coarse+shift);
            const amrex::Real epg   = epg_coarse_arr(ijk_coarse+shift);

            sum += vfrac*epg;
          }
        }

        for (int ii(-1); ii <= 1; ++ii)
        for (int jj(-1); jj <= 1; ++jj)
        for (int kk(-1); kk <= 1; ++kk) {
          amrex::IntVect shift(ii,jj,kk);
          if (!flags_coarse_arr(ijk_coarse+shift).isCovered() &&
              domain.contains(ijk_coarse+shift)) {

            const amrex::Real vfrac = volfrac_coarse_arr(ijk_coarse+shift);
            const amrex::Real epg   = epg_coarse_arr(ijk_coarse+shift);

            weights[ii+1][jj+1][kk+1] = (vfrac*epg) / sum;

          } else {

            weights[ii+1][jj+1][kk+1] = 0.;
          }
        }

      } else {

        for (int ii(-1); ii <= 1; ++ii)
        for (int jj(-1); jj <= 1; ++jj)
        for (int kk(-1); kk <= 1; ++kk) {

          weights[ii+1][jj+1][kk+1] = 0.;
        }

        weights[1][1][1] = 1.;
      }
    }

    AMREX_GPU_HOST_DEVICE
    amrex::Real sum_up (const amrex::GpuArray<amrex::GpuArray<amrex::GpuArray<amrex::Real,3>,3>,3>& weights,
                        const amrex::Box& domain,
                        const amrex::Array4<const amrex::EBCellFlag>& flags_arr,
                        const amrex::Array4<const amrex::Real>& mf_arr,
                        const amrex::IntVect& ijk,
                        const int n=0) const
    {
      amrex::Real result(0.);

      for (int ii(-1); ii <= 1; ++ii)
      for (int jj(-1); jj <= 1; ++jj)
      for (int kk(-1); kk <= 1; ++kk) {
        amrex::IntVect shift(ii,jj,kk);
        if (!flags_arr(ijk+shift).isCovered() &&
            domain.contains(ijk+shift)) {
          result += weights[ii+1][jj+1][kk+1] * mf_arr(ijk+shift, n);
        }
      }

      return result;
    }
};

#endif
