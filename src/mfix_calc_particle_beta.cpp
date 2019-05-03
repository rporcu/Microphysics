#include <AMReX_ParmParse.H>
#include <mfix_F.H>
#include <mfix_eb_F.H>
#include <mfix.H>
#include <mfix_des_F.H>
#include <mfix_util_F.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_EBMultiFabUtil.H>

void mfix::mfix_calc_particle_beta()
{
    BL_PROFILE("mfix::mfix_calc_particle_beta()");

    Real velfp[3];

    for (int lev = 0; lev < nlev; lev++)
    {
        bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                             (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

        MultiFab* ep_ptr;
        MultiFab* ro_ptr;
        MultiFab* mu_ptr;
        MultiFab* vel_ptr;

        // This will be temporaries only for the dual grid case
        std::unique_ptr<MultiFab>  ep_g_pba;
        std::unique_ptr<MultiFab>  ro_g_pba;
        std::unique_ptr<MultiFab>  mu_g_pba;
        std::unique_ptr<MultiFab> vel_g_pba;
        std::unique_ptr<MultiFab>  drag_pba;

        if (OnSameGrids)
        {
            ep_ptr    =  ep_g[lev].get();
            ro_ptr    =  ro_g[lev].get();
            mu_ptr    =  mu_g[lev].get();
            vel_ptr   = vel_g[lev].get();
        }
        else
        {
            BoxArray            pba = pc->ParticleBoxArray(lev);
            DistributionMapping pdm = pc->ParticleDistributionMap(lev);

            // Temporary arrays  -- CHECK COPY() ROUTINE
            int ng = ep_g[lev]->nGrow();
            ep_g_pba.reset(new MultiFab(pba,pdm,ep_g[lev]->nComp(),ng));
            ep_g_pba->copy(*ep_g[lev],0,0,1,ng,ng);
            ep_g_pba->FillBoundary(geom[lev].periodicity());

            ng = ro_g[lev]->nGrow();
            ro_g_pba.reset(new MultiFab(pba,pdm,ro_g[lev]->nComp(),ng));
            ro_g_pba->copy(*ro_g[lev],0,0,1,ng,ng);
            ro_g_pba->FillBoundary(geom[lev].periodicity());

            ng = mu_g[lev]->nGrow();
            mu_g_pba.reset(new MultiFab(pba,pdm,mu_g[lev]->nComp(),ng));
            mu_g_pba->copy(*mu_g[lev],0,0,1,ng,ng);
            mu_g_pba->FillBoundary(geom[lev].periodicity());

            EBFArrayBoxFactory ebfactory_loc( * eb_levels[lev], geom[lev], pba, pdm,
                                              {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                               m_eb_full_grow_cells}, EBSupport::basic);

            ng = vel_g[lev]->nGrow();
            vel_g_pba.reset(new MultiFab(pba,pdm,vel_g[lev]->nComp(),ng, MFInfo(), ebfactory_loc));
            vel_g_pba->copy(*vel_g[lev],0,0,vel_g[lev]->nComp(),ng,ng);
            vel_g_pba->FillBoundary(geom[lev].periodicity());

            drag_pba.reset(new MultiFab(pba,pdm,drag[lev]->nComp(),drag[lev]->nGrow()));

            ep_ptr    =  ep_g_pba.get();
            ro_ptr    =  ro_g_pba.get();
            mu_ptr    =  mu_g_pba.get();
            vel_ptr   = vel_g_pba.get();
        }

        // Phi is always on the particles grid
        const MultiFab & phi = * level_sets[lev];

        int band_width = 2;

        Real x0 = geom[lev].ProbLo(0);
        Real y0 = geom[lev].ProbLo(1);
        Real z0 = geom[lev].ProbLo(2);

        Real odx = 1./geom[lev].CellSize(0);
        Real ody = 1./geom[lev].CellSize(1);
        Real odz = 1./geom[lev].CellSize(2);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            // Create fab to host reconstructed velocity field
            FArrayBox vel_r;

            for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
            {
                auto& particles = pti.GetArrayOfStructs();
                const int np = particles.size();

                const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>((*vel_ptr)[pti]);
                const EBFArrayBox&   ep_fab = static_cast<EBFArrayBox const&>(( *ep_ptr)[pti]);
                const EBFArrayBox&   ro_fab = static_cast<EBFArrayBox const&>(( *ro_ptr)[pti]);
                const EBFArrayBox&   mu_fab = static_cast<EBFArrayBox const&>(( *mu_ptr)[pti]);

                // this is to check efficiently if this tile contains any eb stuff
                const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();
                Box bx = pti.tilebox ();

                const auto& vel_array = vel_ptr->array(pti);
                const auto&  ep_array =  ep_ptr->array(pti);
                const auto&  ro_array =  ro_ptr->array(pti);
                const auto&  mu_array =  mu_ptr->array(pti);

                if (flags.getType(amrex::grow(bx,0)) != FabType::covered)
                {
                    if (flags.getType(amrex::grow(bx,1)) == FabType::regular)
                    {
                        // calc_particle_beta( BL_TO_FORTRAN_ANYD((*ep_ptr)[pti]),
                        //                     (*ro_ptr)[pti].dataPtr(),
                        //                     (*mu_ptr)[pti].dataPtr(),
                        //                     BL_TO_FORTRAN_ANYD((*vel_ptr)[pti]),
                        //                     &np, particles.data(),
                        //                     geom[lev].ProbLo(), geom[lev].CellSize());

                        for (int ip = 0; ip < np; ++ip)
                        {
                            MFIXParticleContainer::ParticleType& particle = particles[ip];

                            // This part was in trilinear_interp before
                            // Pick upper cell in the stencil
                            Real lx = (particle.pos(0) - x0)*odx + 0.5;
                            Real ly = (particle.pos(1) - y0)*ody + 0.5;
                            Real lz = (particle.pos(2) - z0)*odz + 0.5;
    
                            int i = std::floor(lx);
                            int j = std::floor(ly);
                            int k = std::floor(lz);
    
                            // If the particle has gone outside the fluid region then we will
                            //    let the wall collision term bring it back and not use any
                            //    extrapolated fluid quantities in a covered region
                            // Weights
                            Real sx_hi = lx - i;  Real sx_lo = 1.0 - sx_hi;
                            Real sy_hi = ly - j;  Real sy_lo = 1.0 - sy_hi;
                            Real sz_hi = lz - k;  Real sz_lo = 1.0 - sz_hi;

                            for (int n = 0; n < 3; n++)
                               velfp[n] = sx_lo*sy_lo*sz_lo*vel_array(i-1, j-1, k-1,n) +
                                          sx_lo*sy_lo*sz_hi*vel_array(i-1, j-1, k  ,n) +
                                          sx_lo*sy_hi*sz_lo*vel_array(i-1, j  , k-1,n) +
                                          sx_lo*sy_hi*sz_hi*vel_array(i-1, j  , k  ,n) +
                                          sx_hi*sy_lo*sz_lo*vel_array(i  , j-1, k-1,n) +
                                          sx_hi*sy_lo*sz_hi*vel_array(i  , j-1, k  ,n) +
                                          sx_hi*sy_hi*sz_lo*vel_array(i  , j  , k-1,n) +
                                          sx_hi*sy_hi*sz_hi*vel_array(i  , j  , k  ,n);

                            // Indices of cell where particle is located
                            int iloc = floor((particle.pos(0) - x0)*odx);
                            int jloc = floor((particle.pos(1) - y0)*ody);
                            int kloc = floor((particle.pos(2) - z0)*odz);

                            Real  ep = ep_array(iloc,jloc,kloc);
                            Real  ro = ro_array(iloc,jloc,kloc);
                            Real  mu = mu_array(iloc,jloc,kloc);

                            Real rad = particle.rdata(realData::radius);
                            Real vol = particle.rdata(realData::volume);
                            Real den = particle.rdata(realData::density);

                            int p_id = particle.id();

                            Real pvel[3]; 
                            pvel[0] = particle.rdata(realData::velx);
                            pvel[1] = particle.rdata(realData::vely);
                            pvel[2] = particle.rdata(realData::velz);

                            Real beta; 
                            des_drag_gp(&p_id, pvel, velfp, &ep, 
                                  &ro, &mu, &beta, &iloc, &jloc, &kloc, &rad, &vol, &den);

                            particle.rdata(realData::dragx) = beta;
                        }
                    }
                    else
                    {
                        Box gbox = amrex::grow(bx,2);
                        vel_r.resize(gbox,3);

                        reconstruct_velocity( BL_TO_FORTRAN_ANYD(vel_r),
                                              BL_TO_FORTRAN_ANYD((*vel_ptr)[pti]),
                                              BL_TO_FORTRAN_ANYD((phi)[pti]),
                                              1,
                                              BL_TO_FORTRAN_ANYD(flags),
                                              geom[lev].ProbLo(), geom[lev].CellSize(),
                                              &band_width);

                        // calc_particle_beta( BL_TO_FORTRAN_ANYD((*ep_ptr)[pti]),
                        //                     (*ro_ptr)[pti].dataPtr(),
                        //                     (*mu_ptr)[pti].dataPtr(),
                        //                     BL_TO_FORTRAN_ANYD(vel_r),
                        //                     &np, particles.data(),
                        //                     geom[lev].ProbLo(), geom[lev].CellSize());

                        const auto& vel_array = vel_r.array();

                        for (int ip = 0; ip < np; ++ip)
                        {
                            MFIXParticleContainer::ParticleType& particle = particles[ip];

                            // This part was in trilinear_interp before
                            // Pick upper cell in the stencil
                            Real lx = (particle.pos(0) - x0)*odx + 0.5;
                            Real ly = (particle.pos(1) - y0)*ody + 0.5;
                            Real lz = (particle.pos(2) - z0)*odz + 0.5;
    
                            int i = std::floor(lx);
                            int j = std::floor(ly);
                            int k = std::floor(lz);
    
                            // If the particle has gone outside the fluid region then we will
                            //    let the wall collision term bring it back and not use any
                            //    extrapolated fluid quantities in a covered region
                            // Weights
                            Real sx_hi = lx - i;  Real sx_lo = 1.0 - sx_hi;
                            Real sy_hi = ly - j;  Real sy_lo = 1.0 - sy_hi;
                            Real sz_hi = lz - k;  Real sz_lo = 1.0 - sz_hi;

                            for (int n = 0; n < 3; n++)
                               velfp[n] = sx_lo*sy_lo*sz_lo*vel_array(i-1, j-1, k-1,n) +
                                          sx_lo*sy_lo*sz_hi*vel_array(i-1, j-1, k  ,n) +
                                          sx_lo*sy_hi*sz_lo*vel_array(i-1, j  , k-1,n) +
                                          sx_lo*sy_hi*sz_hi*vel_array(i-1, j  , k  ,n) +
                                          sx_hi*sy_lo*sz_lo*vel_array(i  , j-1, k-1,n) +
                                          sx_hi*sy_lo*sz_hi*vel_array(i  , j-1, k  ,n) +
                                          sx_hi*sy_hi*sz_lo*vel_array(i  , j  , k-1,n) +
                                          sx_hi*sy_hi*sz_hi*vel_array(i  , j  , k  ,n);

                            // Indices of cell where particle is located
                            int iloc = floor((particle.pos(0) - x0)*odx);
                            int jloc = floor((particle.pos(1) - y0)*ody);
                            int kloc = floor((particle.pos(2) - z0)*odz);

                            if (std::abs(velfp[0]) > 1.e20) 
                            {
                               std::cout << "BAD LOC " << iloc << " " << jloc << " " << kloc << std::endl;

                               if (std::abs(vel_array(i-1,j-1,k-1,0)) > 1.e20) std::cout << "BAD VEL " << i-1 << " " << j-1 << " " << k-1 << std::endl;
                               if (std::abs(vel_array(i-1,j-1,k  ,0)) > 1.e20) std::cout << "BAD VEL " << i-1 << " " << j-1 << " " << k   << std::endl;
                               if (std::abs(vel_array(i-1,j  ,k-1,0)) > 1.e20) std::cout << "BAD VEL " << i-1 << " " << j   << " " << k-1 << std::endl;
                               if (std::abs(vel_array(i-1,j  ,k  ,0)) > 1.e20) std::cout << "BAD VEL " << i-1 << " " << j   << " " << k   << std::endl;

                               if (std::abs(vel_array(i-1,j-1,k-1,0)) > 1.e20) std::cout << "BAD VEL " << i   << " " << j-1 << " " << k-1 << std::endl;
                               if (std::abs(vel_array(i  ,j-1,k  ,0)) > 1.e20) std::cout << "BAD VEL " << i   << " " << j-1 << " " << k   << std::endl;
                               if (std::abs(vel_array(i  ,j  ,k-1,0)) > 1.e20) std::cout << "BAD VEL " << i   << " " << j   << " " << k-1 << std::endl;
                               if (std::abs(vel_array(i  ,j  ,k  ,0)) > 1.e20) std::cout << "BAD VEL " << i   << " " << j   << " " << k   << std::endl;

                               std::cout << " " << std::endl;
                            }

                            Real  ep = ep_array(iloc,jloc,kloc);
                            Real  ro = ro_array(iloc,jloc,kloc);
                            Real  mu = mu_array(iloc,jloc,kloc);

                            Real rad = particle.rdata(realData::radius);
                            Real vol = particle.rdata(realData::volume);
                            Real den = particle.rdata(realData::density);

                            int p_id = particle.id();

                            Real pvel[3]; 
                            pvel[0] = particle.rdata(realData::velx);
                            pvel[1] = particle.rdata(realData::vely);
                            pvel[2] = particle.rdata(realData::velz);

                            Real beta; 
                            des_drag_gp(&p_id, pvel, velfp, &ep, &ro, &mu, 
                                        &beta, &iloc, &jloc, &kloc, &rad, &vol, &den);

                            particle.rdata(realData::dragx) = beta;
                        }
                    }
                }

            }
        }
    } // lev
}
