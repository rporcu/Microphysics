#include <AMReX.H>
#include "AMReX_Particles.H"
#include "AMReX_RealVect.H"
#include <iostream>
#include <MFIXParticleContainer.H>
#include <AMReX_LoadBalanceKD.H>

#include<math.h>

#include "mfix_F.H"

using namespace amrex;
using namespace std;

bool MFIXParticleContainer::use_neighbor_list  {true};
bool MFIXParticleContainer::sort_neighbor_list {false};

MFIXParticleContainer::MFIXParticleContainer (AmrCore* amr_core)
    : NeighborParticleContainer<realData::count,intData::count,realData::count+2>
      (amr_core->GetParGDB(), 1)
{
    ReadStaticParameters();

    this->SetVerbose(0);
}

void MFIXParticleContainer::AllocData ()
{
    reserveData();
    resizeData();
}

void MFIXParticleContainer::InitParticlesAscii(const std::string& file) {

  // only read the file on the IO proc
  if (ParallelDescriptor::IOProcessor())  {
    std::ifstream ifs;
    ifs.open(file.c_str(), std::ios::in);

    if (!ifs.good())
      amrex::FileOpenFailed(file);

    int np = -1;
    ifs >> np >> std::ws;

    // Issue an error if nparticles = 0 is specified
    if ( np == -1 ){
      Abort("\nCannot read number of particles from particle_input.dat: file is corrupt.\
                   \nPerhaps you forgot to specify the number of particles on the first line??? ");
    }

    // we add all the particles to grid 0 and tile 0 and let
    // Redistribute() put them in the right places.
    const int lev  = 0;
    const int grid = 0;
    const int tile = 0;

    auto& particle_tile = GetParticles(lev)[std::make_pair(grid,tile)];

    ParticleType p;
    int        pstate, pphase;
    Real       pradius, pdensity, pvolume, pomoi, pmass, pomega;

    pstate = 1;

    for (int i = 0; i < np; i++) {

      // Read from input file
      ifs >> pphase;
      ifs >> p.pos(0);
      ifs >> p.pos(1);
      ifs >> p.pos(2);
      ifs >> pradius;
      ifs >> pdensity;
      ifs >> p.rdata(realData::velx);
      ifs >> p.rdata(realData::vely);
      ifs >> p.rdata(realData::velz);

      // Set id and cpu for this particle
      p.id()  = ParticleType::NextID();
      p.cpu() = ParallelDescriptor::MyProc();

      // Compute other particle properties
      set_particle_properties( &pstate, &pradius, &pdensity,
                               &pvolume, &pmass, &pomoi, &pomega);

      // Set other particle properties
      p.idata(intData::phase)     = pphase;
      p.idata(intData::state)     = pstate;
      p.rdata(realData::volume)   = pvolume;
      p.rdata(realData::density)  = pdensity;
      p.rdata(realData::mass)     = pmass;
      p.rdata(realData::oneOverI) = pomoi;
      p.rdata(realData::radius)   = pradius;
      p.rdata(realData::omegax)   = pomega;
      p.rdata(realData::omegay)   = pomega;
      p.rdata(realData::omegaz)   = pomega;

      // Initialize these for I/O purposes
      p.rdata(realData::dragx)    = 0.0;
      p.rdata(realData::dragy)    = 0.0;
      p.rdata(realData::dragz)    = 0.0;

      // Add everything to the data structure
      particle_tile.push_back(p);
    }
  }
  Redistribute();
}

void MFIXParticleContainer::Replicate(IntVect& Nrep, Geometry& geom, DistributionMapping& dmap, BoxArray& ba)
{
    int lev = 0;

    Array<Real> orig_domain_size;
    orig_domain_size.resize(BL_SPACEDIM);
    for (int d = 0; d < BL_SPACEDIM; d++) 
        orig_domain_size[d] = (geom.ProbHi(d) - geom.ProbLo(d)) / Nrep[d];

    ParticleType p_rep;

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& particles = pti.GetArrayOfStructs();

        for (const auto& p: particles)
        {
           //
           // Shift the position.
           //
           for (int k = 0; k < Nrep[2]; k++) {
               for (int j = 0; j < Nrep[1]; j++) {
                 for (int i = 0; i < Nrep[0]; i++) {

                   if ( !(i == 0 && j == 0 && k == 0) )
                   {
                    p_rep.m_rdata.pos[0] = p.m_rdata.pos[0] + i * orig_domain_size[0];
                    p_rep.m_rdata.pos[1] = p.m_rdata.pos[1] + j * orig_domain_size[1];
                    p_rep.m_rdata.pos[2] = p.m_rdata.pos[2] + k * orig_domain_size[2];

                    p_rep.rdata(realData::velx)   = p.rdata(realData::velx);
                    p_rep.rdata(realData::vely)   = p.rdata(realData::vely);
                    p_rep.rdata(realData::velz)   = p.rdata(realData::velz);

                    // Set other particle properties
                    p_rep.idata(intData::phase)     = p.idata(intData::phase);
                    p_rep.idata(intData::state)     = p.idata(intData::state);
                    p_rep.rdata(realData::volume)   = p.rdata(realData::volume);
                    p_rep.rdata(realData::density)  = p.rdata(realData::density);
                    p_rep.rdata(realData::mass)     = p.rdata(realData::mass);
                    p_rep.rdata(realData::oneOverI) = p.rdata(realData::oneOverI);
                    p_rep.rdata(realData::radius)   = p.rdata(realData::radius);
                    p_rep.rdata(realData::omegax)   = p.rdata(realData::omegax);
                    p_rep.rdata(realData::omegay)   = p.rdata(realData::omegay);
                    p_rep.rdata(realData::omegaz)   = p.rdata(realData::omegaz);
                    p_rep.rdata(realData::dragx)    = p.rdata(realData::dragx);
                    p_rep.rdata(realData::dragy)    = p.rdata(realData::dragy);
                    p_rep.rdata(realData::dragz)    = p.rdata(realData::dragz);

                    // Set id and cpu for this particle
                    p_rep.id()  = ParticleType::NextID();
                    p_rep.cpu() = ParallelDescriptor::MyProc();

                    // Add everything to the data structure
                    particles.push_back(p_rep);

                   } // not copying itself
                 } // i
              } // j
           } // k
        } // p
    } // pti

    BuildLevelMask(lev,geom,dmap,ba);
    Redistribute();
}

void MFIXParticleContainer:: printParticles()
{
    const int lev = 0;
    const auto& plevel = GetParticles(lev);

    for (const auto& kv : plevel)
    {
       const auto& particles = kv.second.GetArrayOfStructs();

       for (unsigned i = 0; i < particles.numParticles(); ++i)
       {
          std::cout << "Particle ID  = " << i << " " << std::endl;
          std::cout << "X            = " << particles[i].pos(0) << " " << std::endl;
          std::cout << "Y            = " << particles[i].pos(1) << " " << std::endl;
          std::cout << "Z            = " << particles[i].pos(2) << " " << std::endl;
          std::cout << "state        = " << particles[i].idata(intData::state) << " " << std::endl;
          std::cout << "phase        = " << particles[i].idata(intData::phase) << " " << std::endl;
          std::cout << "Real properties = " << std::endl;

          for (int j = 0; j < realData::count; j++)
            std::cout << "property " << j << "  = " << particles[i].rdata(j) << " " << std::endl;

          std::cout << std::endl;
       }
    }
}

void MFIXParticleContainer::ReadStaticParameters ()
{
    static bool initialized = false;

    if (!initialized)
    {
        ParmParse pp("particles");

        do_tiling = true;  // because the default in amrex is false

        pp.query("do_tiling",  do_tiling);

        Array<int> ts(BL_SPACEDIM);

        if (pp.queryarr("tile_size", ts))
            tile_size = IntVect(ts);

        pp.query("use_neighbor_list", use_neighbor_list);
        pp.query("sort_neighbor_list", sort_neighbor_list);

        initialized = true;
    }
}

void
MFIXParticleContainer::InitData()
{
}

void MFIXParticleContainer::EvolveParticles( int lev, int nstep, Real dt, Real time )
{
    BL_PROFILE("mfix_dem::EvolveParticles()");

    bool debug = false;

    Box domain(Geom(lev).Domain());

    Real dx = Geom(lev).CellSize(0);
    Real dy = Geom(lev).CellSize(1);
    Real dz = Geom(lev).CellSize(2);

    Real xlen = Geom(lev).ProbHi(0) - Geom(lev).ProbLo(0);
    Real ylen = Geom(lev).ProbHi(1) - Geom(lev).ProbLo(1);
    Real zlen = Geom(lev).ProbHi(2) - Geom(lev).ProbLo(2);

    int   nsubsteps;
    Real  subdt;

    des_init_time_loop( &time, &dt, &nsubsteps, &subdt );

    int ncoll_total = 0;

    for ( int n = 0; n < nsubsteps; ++n )
    {
      int ncoll = 0;

      if (use_neighbor_list)
      {
          if (n % 25 == 0) {
              clearNeighbors(lev);
              Redistribute();
              fillNeighbors(lev);
              buildNeighborList(lev,sort_neighbor_list);
          } else {
              updateNeighbors(lev);
          }

#ifdef _OPENMP
#pragma omp parallel
#endif
         for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

            // Real particles
            const int np     = NumberOfParticles(pti);
            void* particles  = pti.GetArrayOfStructs().data();

            // Neighbor particles
            PairIndex index(pti.index(), pti.LocalTileIndex());
            int size_ng = neighbors[index].size() / pdata_size;
            int size_nl = neighbor_list[index].size();

            BL_PROFILE_VAR("des_time_loop()", des_time_loop);
            des_time_loop_ops_nl ( &np, particles, &size_ng, neighbors[index].dataPtr(),
                                   &size_nl, neighbor_list[index].dataPtr(),
                                   &subdt, &dx, &dy, &dz,
                                   &xlen, &ylen, &zlen, &n, &ncoll );
            BL_PROFILE_VAR_STOP(des_time_loop);

            if ( des_continuum_coupled () == 0 ) {
              Real stime;
              stime = time + (n+1)*subdt;
              output_manager( &np, &stime, &subdt,  &xlen, &ylen, &zlen,
                                   &n, particles, 0 );

            }

            call_usr2_des( &np, particles );
         }

      } else {

         fillNeighbors(lev);

#ifdef _OPENMP
#pragma omp parallel
#endif
         for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

            // Real particles
            const int np     = NumberOfParticles(pti);
            void* particles  = pti.GetArrayOfStructs().data();

            // Neighbor particles
            PairIndex index(pti.index(), pti.LocalTileIndex());
            int ng = neighbors[index].size() / pdata_size;

            BL_PROFILE_VAR("des_time_loop()", des_time_loop);
            des_time_loop_ops( &np, particles, &ng, neighbors[index].dataPtr(),
                               &subdt, &dx, &dy, &dz,
                               &xlen, &ylen, &zlen, &n, &ncoll );
            BL_PROFILE_VAR_STOP(des_time_loop);

            if ( des_continuum_coupled () == 0 ) {
              Real stime;
              stime = time + (n+1)*subdt;
              output_manager( &np, &stime, &subdt,  &xlen, &ylen, &zlen,
                                   &n, particles, 0 );

            }

            call_usr2_des( &np, particles );
         }

         clearNeighbors(lev);
         Redistribute();

      }

      if (debug) {
         ncoll_total +=  ncoll;
         ParallelDescriptor::ReduceIntSum(ncoll,ParallelDescriptor::IOProcessorNumber());
         Print() << "Number of collisions: " << ncoll << " at step " << n << std::endl;
      }
    }

    clearNeighbors(lev);

    Redistribute();

    if (debug) {
       ParallelDescriptor::ReduceIntSum(ncoll_total,ParallelDescriptor::IOProcessorNumber());
       Print() << "Number of collisions: " << ncoll_total << " in " << nsubsteps << " substeps " << std::endl;
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

      const int np     = NumberOfParticles(pti);
      void* particles  = pti.GetArrayOfStructs().data();

      call_usr3_des( &np, particles );
    }

    if (debug) {
       // Check maxmium particle velocties at each fluid time step
       // with the goal of veriftying a particle cannot travel more than
       // a single cell per fluid time step.
       Real max_vel_x = 0.0;
       Real max_vel_y = 0.0;
       Real max_vel_z = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(max:max_vel_x,max_vel_y,max_vel_z)
#endif
       for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

         auto& particles = pti.GetArrayOfStructs();

         for (const auto& p: particles)
         {
             max_vel_x = std::max(Real(p.rdata(realData::velx)), max_vel_x);
             max_vel_y = std::max(Real(p.rdata(realData::vely)), max_vel_y);
             max_vel_z = std::max(Real(p.rdata(realData::velz)), max_vel_z);
         }
       }

       ParallelDescriptor::ReduceRealMax(max_vel_x,ParallelDescriptor::IOProcessorNumber());
       ParallelDescriptor::ReduceRealMax(max_vel_y,ParallelDescriptor::IOProcessorNumber());
       ParallelDescriptor::ReduceRealMax(max_vel_z,ParallelDescriptor::IOProcessorNumber());

       if (ParallelDescriptor::IOProcessor())
       {
          const Real* dx = Geom(0).CellSize();
          cout << "Maximum possible distance traveled:" << endl;
          cout <<  "x=  " << max_vel_x * dt
               << " y=  " << max_vel_y * dt
               << " z=  " << max_vel_z * dt  << " and note that "
               << " dx= " << dx[0] << endl;
       }
    }

    if ( des_continuum_coupled () != 0 ) {
      nstep = nsubsteps;
      time  = time + nsubsteps * subdt ;
    }

    // Redistribute();
}

void MFIXParticleContainer::CalcVolumeFraction(amrex::MultiFab& mf_to_be_filled,
                                               IArrayBox& bc_ilo, IArrayBox& bc_ihi,
                                               IArrayBox& bc_jlo, IArrayBox& bc_jhi,
                                               IArrayBox& bc_klo, IArrayBox& bc_khi)
{
    int fortran_volume_comp = 5;
    PICDeposition(mf_to_be_filled, bc_ilo, bc_ihi, bc_jlo, bc_jhi, bc_klo,bc_khi,
                  fortran_volume_comp);

    // Now define this mf = (1 - particle_vol)
    mf_to_be_filled.mult(-1.0,mf_to_be_filled.nGrow());
    mf_to_be_filled.plus( 1.0,mf_to_be_filled.nGrow());
}

void MFIXParticleContainer::CalcDragOnFluid(amrex::MultiFab& beta_x_mf,
                                            amrex::MultiFab& beta_y_mf,
                                            amrex::MultiFab& beta_z_mf,
                                            amrex::MultiFab& beta_u_mf,
                                            amrex::MultiFab& beta_v_mf,
                                            amrex::MultiFab& beta_w_mf,
                                            IArrayBox& bc_ilo, IArrayBox& bc_ihi,
                                            IArrayBox& bc_jlo, IArrayBox& bc_jhi,
                                            IArrayBox& bc_klo, IArrayBox& bc_khi)
{
    int fortran_beta_comp = 15;
    int fortran_vel_comp  =  9;
    PICMultiDeposition(beta_x_mf, beta_y_mf, beta_z_mf,
                       beta_u_mf, beta_v_mf, beta_w_mf,
                       bc_ilo, bc_ihi, bc_jlo, bc_jhi, bc_klo,bc_khi,
                       fortran_beta_comp, fortran_vel_comp);
    if (beta_x_mf.contains_nan())
    {
        std::cout << "BETA_X HAS NANS AFTER SOLVE" << std::endl;
        exit(0);
    }
    if (beta_y_mf.contains_nan())
    {
        std::cout << "BETA_Y HAS NANS AFTER SOLVE" << std::endl;
        exit(0);
    }
    if (beta_z_mf.contains_nan())
    {
        std::cout << "BETA_Z HAS NANS AFTER SOLVE" << std::endl;
        exit(0);
    }
    if (beta_u_mf.contains_nan())
    {
        std::cout << "BETA_U HAS NANS AFTER SOLVE" << std::endl;
        exit(0);
    }
    if (beta_v_mf.contains_nan())
    {
        std::cout << "BETA_V HAS NANS AFTER SOLVE" << std::endl;
        exit(0);
    }
    if (beta_w_mf.contains_nan())
    {
        std::cout << "BETA_W HAS NANS AFTER SOLVE" << std::endl;
        exit(0);
    }
}

void MFIXParticleContainer::PICDeposition(amrex::MultiFab& mf_to_be_filled,
                                          IArrayBox& bc_ilo, IArrayBox& bc_ihi,
                                          IArrayBox& bc_jlo, IArrayBox& bc_jhi,
                                          IArrayBox& bc_klo, IArrayBox& bc_khi,
                                          int fortran_particle_comp)
{
    BL_PROFILE("MFIXParticleContainer::PICDeposition()");

    int   lev = 0;
    int ncomp = 1;

    MultiFab* mf_pointer;

    if (OnSameGrids(lev, mf_to_be_filled)) {
      // If we are already working with the internal mf defined on the
      // particle_box_array, then we just work with this.
      mf_pointer = &mf_to_be_filled;
    }
    else {
      // If mf_to_be_filled is not defined on the particle_box_array, then we need
      // to make a temporary here and copy into mf_to_be_filled at the end.
      mf_pointer = new MultiFab(ParticleBoxArray(lev),
        ParticleDistributionMap(lev),
        ncomp, mf_to_be_filled.nGrow());
    }

    // We must have ghost cells for each FAB so that a particle in one grid can spread
    // its effect to an adjacent grid by first putting the value into ghost cells of its
    // own grid.  The mf->sumBoundary call then adds the value from one grid's ghost cell
    // to another grid's valid region.
    if (mf_pointer->nGrow() < 1)
       amrex::Error("Must have at least one ghost cell when in CalcVolumeFraction");

    const Real      strttime    = ParallelDescriptor::second();
    const Geometry& gm          = Geom(lev);
    const Real*     plo         = gm.ProbLo();
    const Real*     dx          = gm.CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*mf_pointer, true); mfi.isValid(); ++mfi) {
        (*mf_pointer)[mfi].setVal(0);
    }

    using ParConstIter = ParConstIter<realData::count,intData::count,0,0>;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox local_vol;
        for (ParConstIter pti(*this, lev); pti.isValid(); ++pti) {
            const auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;
            const long np = pti.numParticles();
            FArrayBox& fab = (*mf_pointer)[pti];
            const Box& box = fab.box();
            Real* data_ptr;
            const int *lo, *hi;
#ifdef _OPENMP
            Box tile_box = pti.tilebox();
            tile_box.grow(1);
            local_vol.resize(tile_box,ncomp);
            local_vol = 0.0;
            data_ptr = local_vol.dataPtr();
            lo = tile_box.loVect();
            hi = tile_box.hiVect();
#else
            data_ptr = fab.dataPtr();
            lo = box.loVect();
            hi = box.hiVect();
#endif

            mfix_deposit_cic(particles.data(), nstride, np, ncomp, data_ptr,
                             lo, hi, plo, dx, &fortran_particle_comp);

#ifdef _OPENMP
            amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_vol),
                                        BL_TO_FORTRAN_3D(fab), ncomp);
#endif

        }
    }

    // Move any field deposited outside the domain back into the domain
    // (at all domain boundaries except periodic)

    Box domain(Geom(lev).Domain());

    for (MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi) {

      const Box& sbx = (*mf_pointer)[mfi].box();

      flip_particle_vol(sbx.loVect(), sbx.hiVect(),
                        (*mf_pointer)[mfi].dataPtr(),
                        bc_ilo.dataPtr(), bc_ihi.dataPtr(),
                        bc_jlo.dataPtr(), bc_jhi.dataPtr(),
                        bc_klo.dataPtr(), bc_khi.dataPtr(),
                        domain.loVect(), domain.hiVect());
    }

    mf_pointer->SumBoundary(gm.periodicity());

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled. I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf_to_be_filled) {
      mf_to_be_filled.copy(*mf_pointer,0,0,ncomp);
      delete mf_pointer;
    }

    if (m_verbose > 1) {
      Real stoptime = ParallelDescriptor::second() - strttime;

      ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

      amrex::Print() << "MFIXParticleContainer::PICDeposition time: " << stoptime << '\n';
    }
}

void MFIXParticleContainer::PICMultiDeposition(amrex::MultiFab& beta_x_mf,
                                               amrex::MultiFab& beta_y_mf,
                                               amrex::MultiFab& beta_z_mf,
                                               amrex::MultiFab& beta_u_mf,
                                               amrex::MultiFab& beta_v_mf,
                                               amrex::MultiFab& beta_w_mf,
                                               IArrayBox& bc_ilo, IArrayBox& bc_ihi,
                                               IArrayBox& bc_jlo, IArrayBox& bc_jhi,
                                               IArrayBox& bc_klo, IArrayBox& bc_khi,
                                               int fortran_beta_comp, int fortran_vel_comp)
{
    BL_PROFILE("MFIXParticleContainer::PICMultiDeposition()");

    int   lev = 0;
    int ncomp = 1+BL_SPACEDIM;

    MultiFab *beta_x_ptr, *beta_y_ptr, *beta_z_ptr;
    MultiFab *beta_u_ptr, *beta_v_ptr, *beta_w_ptr;

    // Make temporaries here and copy into beta_mf and beta_vel_mf at the end.

    BoxArray x_face_ba = ParticleBoxArray(lev);
    x_face_ba.surroundingNodes(0);
    beta_x_ptr = new MultiFab(x_face_ba, ParticleDistributionMap(lev), 1, beta_x_mf.nGrow());
    beta_u_ptr = new MultiFab(x_face_ba, ParticleDistributionMap(lev), 1, beta_x_mf.nGrow());

    BoxArray y_face_ba = ParticleBoxArray(lev);
    y_face_ba.surroundingNodes(1);
    beta_y_ptr = new MultiFab(y_face_ba, ParticleDistributionMap(lev), 1, beta_y_mf.nGrow());
    beta_v_ptr = new MultiFab(y_face_ba, ParticleDistributionMap(lev), 1, beta_y_mf.nGrow());

    BoxArray z_face_ba = ParticleBoxArray(lev);
    z_face_ba.surroundingNodes(2);
    beta_z_ptr = new MultiFab(z_face_ba, ParticleDistributionMap(lev), 1, beta_z_mf.nGrow());
    beta_w_ptr = new MultiFab(z_face_ba, ParticleDistributionMap(lev), 1, beta_z_mf.nGrow());

    const Real      strttime    = ParallelDescriptor::second();
    const Geometry& gm          = Geom(lev);
    const Real*     plo         = gm.ProbLo();
    const Real*     dx          = gm.CellSize();

    beta_x_ptr->setVal(0.0);
    beta_y_ptr->setVal(0.0);
    beta_z_ptr->setVal(0.0);
    beta_u_ptr->setVal(0.0);
    beta_v_ptr->setVal(0.0);
    beta_w_ptr->setVal(0.0);

    using ParConstIter = ParConstIter<realData::count,intData::count,0,0>;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox local_x_vol, local_y_vol, local_z_vol;
        FArrayBox local_u_vol, local_v_vol, local_w_vol;
        for (ParConstIter pti(*this, lev); pti.isValid(); ++pti) {
            const auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;
            const long np = pti.numParticles();
            FArrayBox& beta_x_fab = (*beta_x_ptr)[pti];
            FArrayBox& beta_y_fab = (*beta_y_ptr)[pti];
            FArrayBox& beta_z_fab = (*beta_z_ptr)[pti];
            FArrayBox& beta_u_fab = (*beta_u_ptr)[pti];
            FArrayBox& beta_v_fab = (*beta_v_ptr)[pti];
            FArrayBox& beta_w_fab = (*beta_w_ptr)[pti];

            const Box& x_box = beta_x_fab.box();
            const Box& y_box = beta_y_fab.box();
            const Box& z_box = beta_z_fab.box();
            Real *bx_dataptr, *by_dataptr, *bz_dataptr;
            Real *bu_dataptr, *bv_dataptr, *bw_dataptr;
            const int *lo_x, *hi_x, *lo_y, *hi_y, *lo_z, *hi_z;

#ifdef _OPENMP
            Box tile_xbox = pti.tilebox();
            tile_xbox.surroundingNodes(0);
            tile_xbox.grow(1);

            Box tile_ybox = pti.tilebox();
            tile_ybox.surroundingNodes(1);
            tile_ybox.grow(1);

            Box tile_zbox = pti.tilebox();
            tile_zbox.surroundingNodes(2);
            tile_zbox.grow(1);

            local_x_vol.resize(tile_xbox,ncomp);
            local_u_vol.resize(tile_xbox,ncomp);
            local_y_vol.resize(tile_ybox,ncomp);
            local_v_vol.resize(tile_ybox,ncomp);
            local_z_vol.resize(tile_zbox,ncomp);
            local_w_vol.resize(tile_zbox,ncomp);

            local_x_vol = 0.0; local_y_vol = 0.0; local_z_vol = 0.0;
            local_u_vol = 0.0; local_v_vol = 0.0; local_w_vol = 0.0;

            bx_dataptr = local_x_vol.dataPtr();
            by_dataptr = local_y_vol.dataPtr();
            bz_dataptr = local_z_vol.dataPtr();
            bu_dataptr = local_u_vol.dataPtr();
            bv_dataptr = local_v_vol.dataPtr();
            bw_dataptr = local_w_vol.dataPtr();

            lo_x = tile_xbox.loVect();
            hi_x = tile_xbox.hiVect();
            lo_y = tile_ybox.loVect();
            hi_y = tile_ybox.hiVect();
            lo_z = tile_zbox.loVect();
            hi_z = tile_zbox.hiVect();
#else
            bx_dataptr = beta_x_fab.dataPtr();
            by_dataptr = beta_y_fab.dataPtr();
            bz_dataptr = beta_z_fab.dataPtr();
            bu_dataptr = beta_u_fab.dataPtr();
            bv_dataptr = beta_v_fab.dataPtr();
            bw_dataptr = beta_w_fab.dataPtr();

            lo_x = x_box.loVect();
            hi_x = x_box.hiVect();
            lo_y = y_box.loVect();
            hi_y = y_box.hiVect();
            lo_z = z_box.loVect();
            hi_z = z_box.hiVect();
#endif

            mfix_multi_deposit_cic(particles.data(), nstride, np,
                                   bx_dataptr, by_dataptr, bz_dataptr,
                                   bu_dataptr, bv_dataptr, bw_dataptr,
                                   lo_x, hi_x, lo_y, hi_y, lo_z, hi_z,
                                   plo, dx, &fortran_beta_comp, &fortran_vel_comp);
//          std::cout << "BETAX " << beta_x_fab << std::endl;

#ifdef _OPENMP
            amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_x_vol),
                                        BL_TO_FORTRAN_3D(beta_x_fab), 1);
            amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_y_vol),
                                        BL_TO_FORTRAN_3D(beta_y_fab), 1);
            amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_z_vol),
                                        BL_TO_FORTRAN_3D(beta_z_fab), 1);
            amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_u_vol),
                                        BL_TO_FORTRAN_3D(beta_u_fab), 1);
            amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_v_vol),
                                        BL_TO_FORTRAN_3D(beta_v_fab), 1);
            amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_w_vol),
                                        BL_TO_FORTRAN_3D(beta_w_fab), 1);
#endif

        }
    }

    beta_x_ptr->SumBoundary(gm.periodicity());
    beta_y_ptr->SumBoundary(gm.periodicity());
    beta_z_ptr->SumBoundary(gm.periodicity());
    beta_u_ptr->SumBoundary(gm.periodicity());
    beta_v_ptr->SumBoundary(gm.periodicity());
    beta_w_ptr->SumBoundary(gm.periodicity());

    // Copy back from mf_pointer
    beta_x_mf.copy(*beta_x_ptr,0,0,1);
    beta_y_mf.copy(*beta_y_ptr,0,0,1);
    beta_z_mf.copy(*beta_z_ptr,0,0,1);
    beta_u_mf.copy(*beta_u_ptr,0,0,1);
    beta_v_mf.copy(*beta_v_ptr,0,0,1);
    beta_w_mf.copy(*beta_w_ptr,0,0,1);

    delete beta_x_ptr;
    delete beta_y_ptr; 
    delete beta_z_ptr;
    delete beta_u_ptr; 
    delete beta_v_ptr; 
    delete beta_w_ptr;

    Box domain(Geom(lev).Domain());

    // Make sure there is zero normal drag force on walls, and reflect the tangential
    // drag force around slip walls only
    for (MFIter mfi(beta_u_mf); mfi.isValid(); ++mfi) {

      const Box& xbx = beta_u_mf[mfi].box();
      const Box& ybx = beta_v_mf[mfi].box();
      const Box& zbx = beta_w_mf[mfi].box();

      set_drag_bcs(xbx.loVect(), xbx.hiVect(),
                   ybx.loVect(), ybx.hiVect(),
                   zbx.loVect(), zbx.hiVect(),
                   beta_u_mf[mfi].dataPtr(),
                   beta_v_mf[mfi].dataPtr(),
                   beta_w_mf[mfi].dataPtr(),
                   bc_ilo.dataPtr(), bc_ihi.dataPtr(),
                   bc_jlo.dataPtr(), bc_jhi.dataPtr(),
                   bc_klo.dataPtr(), bc_khi.dataPtr(),
                   domain.loVect(), domain.hiVect());
    }

    if (m_verbose > 1) {
      Real stoptime = ParallelDescriptor::second() - strttime;

      ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

      amrex::Print() << "MFIXParticleContainer::PICMultiDeposition time: " << stoptime << '\n';
    }
}

void MFIXParticleContainer::output(int lev, int estatus, int finish, int nstep, Real dt, Real time)
{

    Real xlen = Geom(lev).ProbHi(0) - Geom(lev).ProbLo(0);
    Real ylen = Geom(lev).ProbHi(1) - Geom(lev).ProbLo(1);
    Real zlen = Geom(lev).ProbHi(2) - Geom(lev).ProbLo(2);

    // Not threaded because its writing output
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {

      //number of particles
      const int     np = NumberOfParticles(pti);
      void* particles  = pti.GetArrayOfStructs().data();

      output_manager( &np, &time, &dt, &xlen, &ylen, &zlen, &nstep,
                           particles, &finish);
    }

}

void MFIXParticleContainer::writeAllAtLevel(int lev)
{
    // Not threaded because its print to terminal
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& particles = pti.GetArrayOfStructs();

        for (const auto& p: particles)
        {
           const IntVect& iv = Index(p, lev);

           RealVect xyz(p.pos(0), p.pos(1), p.pos(2));
           cout << " id " << p.id()
                << " index " << iv
                << " position " << xyz << endl;
       }
    }
}

void MFIXParticleContainer::writeAllForComparison(int lev)
{
  int Np_tot = 0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:Np_tot)
#endif
  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
      Np_tot += pti.numParticles();

  ParallelDescriptor::ReduceIntSum(Np_tot,ParallelDescriptor::IOProcessorNumber());

  cout << Np_tot << std::endl;

  // Not threaded because its print to terminal
  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
  {
      auto& particles = pti.GetArrayOfStructs();

      for (const auto& p: particles)
         cout << p.pos(0)   << " " << p.pos(1)   << " " << p.pos(2) <<  " " <<
                 p.rdata(1) << " " << p.rdata(2) << " " << p.rdata(2) <<  " " <<
                 ParallelDescriptor::MyProc() << " " << p.id() << std::endl;
  }
}

void MFIXParticleContainer::GetParticleAvgProp(int lev,
         Real (&avg_dp)[10], Real (&avg_ro)[10])
{

   // The number of phases was previously hard set at 10, however lowering
   //  this number would make this code faster.
   int num_of_phases_in_use = 10; //Number of different phases being simulated

   // Cycle through the different phases, starting from 1
   for (int phse=1; phse<=num_of_phases_in_use; ++phse){

     Real p_num  = 0.0; //number of particle
     Real p_diam = 0.0; //particle diameters
     Real p_dens = 0.0; //particle density

   #ifdef _OPENMP
   #pragma omp parallel reduction(+:p_num, p_diam, p_dens)
   #endif
     for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

       auto& particles = pti.GetArrayOfStructs();

       for (const auto& p: particles){
         if ( phse==p.idata(intData::phase) ){
           p_num  += 1.0;
           p_diam += p.rdata(realData::radius) * 2.0;
           p_dens += p.rdata(realData::density);
         }
       }
     }

    // A single MPI call passes all three variables
    ParallelDescriptor::ReduceRealSum({p_num,p_diam,p_dens});

   //calculate averages or set = zero if no particles of that phase
   if (p_num==0){
     avg_dp[phse-1] = 0.0;
     avg_ro[phse-1] = 0.0;
   } else {
     avg_dp[phse-1] = p_diam/p_num;
     avg_ro[phse-1] = p_dens/p_num;
   }
 }
}

void
MFIXParticleContainer::BalanceParticleLoad_KDTree()
{
  bool verbose = false;
  BoxArray old_ba = ParticleBoxArray(0);

  if (verbose) 
  {
     amrex::Print() << "Before KDTree BA HAS " << old_ba.size() << " GRIDS " << std::endl;
     if (old_ba.size() < 32) // This is an arbitrary cut-off so we don't spew for large problems
        amrex::Print() << "Before:" << old_ba << std::endl;
  }

  Array<long> num_part;
  num_part = NumberOfParticlesInGrid(0);
  if (verbose) 
  {
     for (int i = 0; i < old_ba.size(); i++)
        amrex::Print() << "NUM PART IN GRID BEFORE " << i << " " << num_part[i] << std::endl;
  }

  BoxArray new_ba;
  Real cell_weight = 0.;
  loadBalanceKD::balance<MFIXParticleContainer>(*this, new_ba, ParallelDescriptor::NProcs(), cell_weight);

  if (verbose) 
  {
     amrex::Print() << "After  KDTree BA HAS " << new_ba.size() << " GRIDS " << std::endl;
     if (new_ba.size() < 32) // This is an arbitrary cut-off so we don't spew for large problems
            amrex::Print() << "After:" << new_ba << std::endl;
  }

  if (verbose) 
  {
     num_part = NumberOfParticlesInGrid(0);
     for (int i = 0; i < old_ba.size(); i++)
        amrex::Print() << "NUM PART IN GRID AFTER  " << i << " " << num_part[i] << std::endl;
  }

  // Create a new DM to go with the new BA
  DistributionMapping new_dm(new_ba);
  
  Regrid(new_dm, new_ba);
}
