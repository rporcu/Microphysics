#ifndef MFIX_DEPOSITION_OP_H_
#define MFIX_DEPOSITION_OP_H_

#include <mfix_pc.H>
#include <mfix_dem.H>
#include <mfix_species.H>
#include <mfix_fluid.H>
#include <mfix_reactions.H>
#include <mfix_algorithm.H>

#include <AMReX.H>
#include <AMReX_Particles.H>


class MFIXDepositionOp
{
  public:
    using ParticleType = MFIXParticleContainer::ParticleType;
    using PairIndex = MFIXParticleContainer::PairIndex;
    using MFIXParIter = MFIXParticleContainer::MFIXParIter;

    // Constructor
    MFIXDepositionOp ()
    {}

    // Destructor
    virtual ~MFIXDepositionOp ()
    {}

    virtual void deposit (int lev,
                          const amrex::Geometry& geom,
                          MFIXParticleContainer* pc,
                          const amrex::MultiFab* volfrac,
                          const amrex::FabArray<amrex::EBCellFlagFab>* flags,
                          amrex::MultiFab* txfr_mf,
                          amrex::MultiFab* eps_mf = nullptr) = 0;
};


class MFIXSolidsVolume : public MFIXDepositionOp
{
  public:
    // Constructor
    MFIXSolidsVolume ()
    {}

    // Destructor
    virtual ~MFIXSolidsVolume ()
    {}

    virtual void deposit (int lev,
                          const amrex::Geometry& geom,
                          MFIXParticleContainer* pc,
                          const amrex::MultiFab* volfrac,
                          const amrex::FabArray<amrex::EBCellFlagFab>* flags,
                          amrex::MultiFab* txfr_mf,
                          amrex::MultiFab* eps_mf = nullptr);

    template <typename F>
    void deposit (F WeightFunc,
                  int lev,
                  const amrex::Geometry& geom,
                  MFIXParticleContainer* pc,
                  const amrex::MultiFab* volfrac,
                  const amrex::FabArray<amrex::EBCellFlagFab>* flags,
                  amrex::MultiFab* txfr_mf,
                  amrex::MultiFab* eps_mf = nullptr);
};


class MFIXInterphaseTxfr : public MFIXDepositionOp
{
  public:
    // Constructor
    MFIXInterphaseTxfr ()
    {}

    // Destructor
    virtual ~MFIXInterphaseTxfr ()
    {}

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

    std::map<PairIndex, amrex::Gpu::DeviceVector<amrex::Real>>& get_aux ()
    { return m_aux; }

  private:
    std::map<PairIndex, amrex::Gpu::DeviceVector<amrex::Real>> m_aux;
};

#endif
