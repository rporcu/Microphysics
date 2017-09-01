#ifndef _MFIXParticleContainer_H_
#define _MFIXParticleContainer_H_


#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_IArrayBox.H>
#include "AMReX_Particles.H"
#include "AMReX_NeighborParticles.H"

using namespace amrex;
using namespace std;

struct realData
{
    enum { // Particle data stored in p.m_data
  radius = 0,
  volume,
  mass,
  density,
  oneOverI,
  velx,
  vely,
  velz,
  omegax,
  omegay,
  omegaz,
  dragx,
  dragy,
  dragz,
  count

    };
};

struct intData
{
    enum { // Particle data stored in p.m_data
  phase = 0,
  state,
  count
    };
};



class MFIXParIter
    : public amrex::ParIter<realData::count,intData::count,0,0>
{
public:
    using amrex::ParIter<realData::count,intData::count,0,0>::ParIter;
};

class MFIXParticleContainer final
    : public amrex::NeighborParticleContainer<realData::count,intData::count,realData::count+2>
{
public:

    // Get number of particles
    int NumberOfParticles(MFIXParIter& pti) { return pti.GetArrayOfStructs().size(); }

    int NumberOfParticles(int& lev, MFIter& mfi)
    {
      const int gridIndex = mfi.index();
      const int tileIndex = mfi.LocalTileIndex();
      auto&     particles = GetParticles(lev)[std::make_pair(gridIndex,tileIndex)];
      auto&     structs   = particles.GetArrayOfStructs();
      return structs.size();
    }

    MFIXParticleContainer (AmrCore* amr_core);

    virtual ~MFIXParticleContainer() {};

    void InitData ();

    void InitParticlesAscii (const std::string& file);

    void Replicate (IntVect& Nrep, amrex::Array<Real>& orig_domain_size);

    void printParticles ();

    void AllocData ();

    void EvolveParticles ( int lev, int nstep, Real dt, Real time );

    void GetParticleAvgProp ( int lev, Real (&avg_dp)[10], Real (&avg_ro)[10]);

    ///
    /// Write out all particles at a level in the format fcompare_par can read
    ///
    void writeAllForComparison(int lev);

    void output ( int lev, int estatus, int finish, int nstep, Real dt, Real time );

    ///
    /// Write out all particles at a level
    ///
    void writeAllAtLevel(int lev = 0);

    ///
    /// Calculate gas volume fraction using a PIC deposition of particle volume
    ///
    void CalcVolumeFraction(amrex::MultiFab& mf,
                            IArrayBox& bc_ilo, IArrayBox& bc_ihi,
                            IArrayBox& bc_jlo, IArrayBox& bc_jhi,
                            IArrayBox& bc_klo, IArrayBox& bc_khi);

    ///
    /// Do the actual PIC deposition for CalcVolumeFraction
    ///
    void PICDeposition(amrex::MultiFab& mf, 
                       IArrayBox& bc_ilo, IArrayBox& bc_ihi,
                       IArrayBox& bc_jlo, IArrayBox& bc_jhi,
                       IArrayBox& bc_klo, IArrayBox& bc_khi,
                       int fortran_particle_comp);

    ///
    /// Calculate particle-induced drag on the fluid
    ///
    void CalcDragOnFluid(amrex::MultiFab& beta_mf_u, 
                         amrex::MultiFab& beta_mf_v,
                         amrex::MultiFab& beta_mf_w,
                         amrex::MultiFab& beta_u_mf,
                         amrex::MultiFab& beta_v_mf,
                         amrex::MultiFab& beta_w_mf,
                         IArrayBox& bc_ilo, IArrayBox& bc_ihi,
                         IArrayBox& bc_jlo, IArrayBox& bc_jhi,
                         IArrayBox& bc_klo, IArrayBox& bc_khi);

    ///
    /// Do the actual PIC deposition for CalcDragOnFluid
    ///
    void PICMultiDeposition(amrex::MultiFab& beta_x_mf, amrex::MultiFab& beta_y_mf, amrex::MultiFab& beta_z_mf, 
                            amrex::MultiFab& beta_u_mf, amrex::MultiFab& beta_v_mf, amrex::MultiFab& beta_w_mf, 
                            IArrayBox& bc_ilo, IArrayBox& bc_ihi,
                            IArrayBox& bc_jlo, IArrayBox& bc_jhi,
                            IArrayBox& bc_klo, IArrayBox& bc_khi,
                            int fortran_beta_comp, int fortran_vel_comp);

private:

    ///
    /// We add a particle to the neighbor list if it's center is within 3 particle radii
    ///
    inline virtual bool check_pair(const ParticleType& p1, const ParticleType& p2) final {
        return AMREX_D_TERM(   (p1.pos(0) - p2.pos(0))*(p1.pos(0) - p2.pos(0)) ,
                             + (p1.pos(1) - p2.pos(1))*(p1.pos(1) - p2.pos(1)) ,
                             + (p1.pos(2) - p2.pos(2))*(p1.pos(2) - p2.pos(2)) ) <= 9.0*p1.rdata(0)*p1.rdata(0);
    }

    static void ReadStaticParameters ();

    static bool use_neighbor_list;
    static bool sort_neighbor_list;

};
#endif
