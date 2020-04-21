#ifndef _MFIX_EB_INTERP_K_H_
#define _MFIX_EB_INTERP_K_H_

#include "MFIXParticleContainer.H"

#include "mfix_eb_newton_solver_K.H"
#include "mfix_eb_node_connectivity_K.H"

/****************************************************************
 *                                                              *
 *                                                              *
 *                                                              *
 *                                                              *
 ***************************************************************/
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void fe_interp(const amrex::RealVect& pos,
               const int ip, const int jp, const int kp,
               const amrex::RealVect& dx,
               const amrex::RealVect& dxi,
               amrex::Array4<const EBCellFlag> const& flags,
               amrex::Array4<const amrex::Real> const& ccent,
               amrex::Array4<const amrex::Real> const& bcent,
               amrex::Array4<const amrex::Real> const& apx,
               amrex::Array4<const amrex::Real> const& apy,
               amrex::Array4<const amrex::Real> const& apz,
               amrex::Array4<amrex::Real const> const& interp_array,
               amrex::Real* interp_vals,
               const int ncomp)
{

    const amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();

    const amrex::Real cdist_x = pos[0] - (ip + 0.5 + ccent(ip,jp,kp,0))*dx[0];
    const amrex::Real cdist_y = pos[1] - (jp + 0.5 + ccent(ip,jp,kp,1))*dx[1];
    const amrex::Real cdist_z = pos[2] - (kp + 0.5 + ccent(ip,jp,kp,2))*dx[2];

    // Distance between particle and cell centoid.
    const amrex::Real cdist =
      std::sqrt(cdist_x*cdist_x + cdist_y*cdist_y + cdist_z*cdist_z);

    // Before doing anything fance, just check that the particle isn't overlapping
    // the cell center. If it is, use the cell value.
    if(cdist < tolerance){

        interp_vals[0] = interp_array(ip,jp,kp,0);
        interp_vals[1] = interp_array(ip,jp,kp,1);
        interp_vals[2] = interp_array(ip,jp,kp,2);
        for(int n=3; n<ncomp; n++)
            interp_vals[n] = interp_array(ip,jp,kp,n);

        return;
    }

    amrex::Real par_dot_EB(2.);
    amrex::Real cent_dot_EB(1.);

    // If the particle is inside a cut-cell, verify that it is on the correct
    // side of the EB before trying to interpolate.
    if(flags(ip,jp,kp).isSingleValued()){

        const Real adx = apx(ip+1,jp  ,kp  ) - apx(ip  ,jp  ,kp  );
        const Real ady = apy(ip  ,jp+1,kp  ) - apy(ip  ,jp  ,kp  );
        const Real adz = apz(ip  ,jp  ,kp+1) - apz(ip  ,jp  ,kp  );

        const Real apnorm = std::sqrt(adx*adx + ady*ady + adz*adz);

        // Temporary sanity check.
        AMREX_ASSERT_WITH_MESSAGE( apnorm > std::numeric_limits<amrex::Real>::epsilon(),
            "Invalid EB normal. Not sure how to interpolate!");

        // EB normal
        const amrex::Real apnorminv = 1. / apnorm;
        const amrex::RealVect normal = {adx*apnorminv, ady*apnorminv, adz*apnorminv};

        // Projection of vector pointing from EB centroid to particle onto EB normal
        par_dot_EB = (pos[0] - (ip + 0.5 + bcent(ip,jp,kp,0))*dx[0])*normal[0] +
                     (pos[1] - (jp + 0.5 + bcent(ip,jp,kp,1))*dx[1])*normal[1] +
                     (pos[2] - (kp + 0.5 + bcent(ip,jp,kp,2))*dx[2])*normal[2];

        // Projection of vector pointing from EB centroid to cell centroid onto EB normal
        cent_dot_EB = (ccent(ip,jp,kp,0) - bcent(ip,jp,kp,0))*dx[0]*normal[0] +
                      (ccent(ip,jp,kp,1) - bcent(ip,jp,kp,1))*dx[1]*normal[1] +
                      (ccent(ip,jp,kp,2) - bcent(ip,jp,kp,2))*dx[2]*normal[2];

        // Temporary sanity check
        AMREX_ASSERT_WITH_MESSAGE( cent_dot_EB > tolerance,
            "cent_dot_EB < tolerance ... this makes no sense!");


    }
    // Particle position relative to cell center [-0.5, 0.5]
    amrex::Real gx = pos[0]*dxi[0] - (ip + 0.5);
    amrex::Real gy = pos[1]*dxi[1] - (jp + 0.5);
    amrex::Real gz = pos[2]*dxi[2] - (kp + 0.5);

    // Use the centoid location of the cell containing the particle
    // to guess the interpolation stencil.
    const int i = (gx < ccent(ip,jp,kp,0)) ? ip : ip + 1;
    const int j = (gy < ccent(ip,jp,kp,1)) ? jp : jp + 1;
    const int k = (gz < ccent(ip,jp,kp,2)) ? kp : kp + 1;

    const int di = i - ip; // 0 or 1
    const int dj = j - jp; // 0 or 1
    const int dk = k - kp; // 0 or 1

    // Count the number of non-conntected cells in the stencil
    int covered = 0;
    for(int kk(-1); kk<1; kk++){
      for(int jj(-1); jj<1; jj++){
        for(int ii(-1); ii<1; ii++){
          if(not flags(ip,jp,kp).isConnected(di+ii,dj+jj,dk+kk))
            covered += 1;
        }
      }
    }

    // A negative value implies that the particle is 'behind' the EB (and
    // therefore inside the wall). Although this shouldn't occur often, it could
    // for fast moving particles.  Set the fluid velocity to zero and use the
    // cell value for remaining items.
    if( flags(ip,jp,kp).isSingleValued() and (par_dot_EB <= tolerance)){

        // Set the velocity to zero
        interp_vals[0] = 0.;
        interp_vals[1] = 0.;
        interp_vals[2] = 0.;
        for(int n=3; n<ncomp; n++){
            interp_vals[n] = interp_array(ip,jp,kp,n);
        }

    // The particle is near the EB. It is either
    // 1) between the cell centroid and the EB, or
    // 2) one or more cells in the 2^3 stencil is covered.
    // Either way, interpolating to the particle is not straight forward so
    // (for now) we do a 1D interpolation to the EB.

    } else if((covered > 0) or (par_dot_EB < cent_dot_EB)){

        amrex::Real scale(1.);

        for(int kk(-1); kk<1; kk++){
            for(int jj(-1); jj<1; jj++){
                for(int ii(-1); ii<1; ii++){
                    if(flags(i+ii,j+jj,k+kk).isSingleValued()) {

                        // EB normal
                        const amrex::Real adx = apx(i+ii+1,j+jj  ,k+kk  ) - apx(i+ii,j+jj,k+kk);
                        const amrex::Real ady = apy(i+ii  ,j+jj+1,k+kk  ) - apy(i+ii,j+jj,k+kk);
                        const amrex::Real adz = apz(i+ii  ,j+jj  ,k+kk+1) - apz(i+ii,j+jj,k+kk);
                        const amrex::Real apnorminv = 1. / std::sqrt(adx*adx + ady*ady + adz*adz);

                        const amrex::RealVect normal = {adx*apnorminv, ady*apnorminv, adz*apnorminv};

                        // Ratio of distance of particle-to-EB and cell-to-EB
                        amrex::Real test_scale = ((pos[0] - (i+ii + 0.5 + bcent(i+ii,j+jj,k+kk,0))*dx[0])*normal[0] +
                                                  (pos[1] - (j+jj + 0.5 + bcent(i+ii,j+jj,k+kk,1))*dx[1])*normal[1] +
                                                  (pos[2] - (k+kk + 0.5 + bcent(i+ii,j+jj,k+kk,2))*dx[2])*normal[2]) /
                                                 ((ccent(ip,jp,kp,0) - (di+ii + bcent(i+ii,j+jj,k+kk,0)))*dx[0]*normal[0] +
                                                  (ccent(ip,jp,kp,1) - (dj+jj + bcent(i+ii,j+jj,k+kk,1)))*dx[1]*normal[1] +
                                                  (ccent(ip,jp,kp,2) - (dk+kk + bcent(i+ii,j+jj,k+kk,2)))*dx[2]*normal[2]);

                        // Save the shortest positive distance.
                        scale = (test_scale > 0.) ? amrex::min(scale,test_scale) : scale;

                    }
                }
            }
        }

        // 1D interpolation from cell center to wall.

        interp_vals[0] = scale*interp_array(ip,jp,kp,0);
        interp_vals[1] = scale*interp_array(ip,jp,kp,1);
        interp_vals[2] = scale*interp_array(ip,jp,kp,2);

        for(int n=3; n<ncomp; n++)
            interp_vals[n] = interp_array(ip,jp,kp,n);

    } else {

        amrex::Real nodes[8][3];
        amrex::Real* values[8];
        for(int it(0); it < 8; it++)
          values[it] = new Real [ncomp];

        // Alias the nodes in the stencil to make the following a little
        // easier to read and debug.
        const amrex::IntVect Node0 = {i-1,j-1,k-1};
        const amrex::IntVect Node1 = {i  ,j-1,k-1};
        const amrex::IntVect Node2 = {i  ,j  ,k-1};
        const amrex::IntVect Node3 = {i-1,j  ,k-1};
        const amrex::IntVect Node4 = {i-1,j-1,k  };
        const amrex::IntVect Node5 = {i  ,j-1,k  };
        const amrex::IntVect Node6 = {i  ,j  ,k  };
        const amrex::IntVect Node7 = {i-1,j  ,k  };

        /*---------------------------------------------------------------------------------*
        *                                                                                  *
        *                                    NODE 0                                        *
        *                                                                                  *
        *----------------------------------------------------------------------------------*/
        nodes[0][0] = (Node0[0] + 0.5 + ccent(Node0,0))*dx[0];
        nodes[0][1] = (Node0[1] + 0.5 + ccent(Node0,1))*dx[1];
        nodes[0][2] = (Node0[2] + 0.5 + ccent(Node0,2))*dx[2];

        for(int n=0; n<ncomp; n++)
            values[0][n] = interp_array(Node0,n);

        /*---------------------------------------------------------------------------------*
        *                                                                                  *
        *                                    NODE 1                                        *
        *                                                                                  *
        *----------------------------------------------------------------------------------*/
        nodes[1][0] = (Node1[0] + 0.5 + ccent(Node1,0))*dx[0];
        nodes[1][1] = (Node1[1] + 0.5 + ccent(Node1,1))*dx[1];
        nodes[1][2] = (Node1[2] + 0.5 + ccent(Node1,2))*dx[2];

        for(int n=0; n<ncomp; n++)
            values[1][n] = interp_array(Node1,n);

        /*---------------------------------------------------------------------------------*
        *                                                                                  *
        *                                    NODE 2                                        *
        *                                                                                  *
        *----------------------------------------------------------------------------------*/
        nodes[2][0] = (Node2[0] + 0.5 + ccent(Node2,0))*dx[0];
        nodes[2][1] = (Node2[1] + 0.5 + ccent(Node2,1))*dx[1];
        nodes[2][2] = (Node2[2] + 0.5 + ccent(Node2,2))*dx[2];

        for(int n=0; n<ncomp; n++)
            values[2][n] = interp_array(Node2,n);

        /*---------------------------------------------------------------------------------*
        *                                                                                  *
        *                                    NODE 3                                        *
        *                                                                                  *
        *----------------------------------------------------------------------------------*/
        nodes[3][0] = (Node3[0] + 0.5 + ccent(Node3,0))*dx[0];
        nodes[3][1] = (Node3[1] + 0.5 + ccent(Node3,1))*dx[1];
        nodes[3][2] = (Node3[2] + 0.5 + ccent(Node3,2))*dx[2];

        for(int n=0; n<ncomp; n++)
            values[3][n] = interp_array(Node3,n);

        /*---------------------------------------------------------------------------------*
        *                                                                                  *
        *                                    NODE 4                                        *
        *                                                                                  *
        *----------------------------------------------------------------------------------*/
        nodes[4][0] = (Node4[0] + 0.5 + ccent(Node4,0))*dx[0];
        nodes[4][1] = (Node4[1] + 0.5 + ccent(Node4,1))*dx[1];
        nodes[4][2] = (Node4[2] + 0.5 + ccent(Node4,2))*dx[2];

        for(int n=0; n<ncomp; n++)
            values[4][n] = interp_array(Node4,n);

        /*---------------------------------------------------------------------------------*
        *                                                                                  *
        *                                    NODE 5                                        *
        *                                                                                  *
        *----------------------------------------------------------------------------------*/
        nodes[5][0] = (Node5[0] + 0.5 + ccent(Node5,0))*dx[0];
        nodes[5][1] = (Node5[1] + 0.5 + ccent(Node5,1))*dx[1];
        nodes[5][2] = (Node5[2] + 0.5 + ccent(Node5,2))*dx[2];

        for(int n=0; n<ncomp; n++)
            values[5][n] = interp_array(Node5,n);

        /*---------------------------------------------------------------------------------*
        *                                                                                  *
        *                                    NODE 6                                        *
        *                                                                                  *
        *----------------------------------------------------------------------------------*/
        nodes[6][0] = (Node6[0] + 0.5 + ccent(Node6,0))*dx[0];
        nodes[6][1] = (Node6[1] + 0.5 + ccent(Node6,1))*dx[1];
        nodes[6][2] = (Node6[2] + 0.5 + ccent(Node6,2))*dx[2];

        for(int n=0; n<ncomp; n++)
            values[6][n] = interp_array(Node6,n);

        /*---------------------------------------------------------------------------------*
        *                                                                                  *
        *                                    NODE 7                                        *
        *                                                                                  *
        *----------------------------------------------------------------------------------*/
        nodes[7][0] = (Node7[0] + 0.5 + ccent(Node7,0))*dx[0];
        nodes[7][1] = (Node7[1] + 0.5 + ccent(Node7,1))*dx[1];
        nodes[7][2] = (Node7[2] + 0.5 + ccent(Node7,2))*dx[2];

        for(int n=0; n<ncomp; n++)
            values[7][n] = interp_array(Node7,n);


        amrex::Real xi;
        amrex::Real eta;
        amrex::Real zeta;

        get_interp_mapping(pos, nodes, dxi, xi, eta, zeta);

// DEBUG START /////////////////////////////////////////////////////////////////
#if 0 && (!defined AMREX_USE_CUDA)
        const amrex::Real lbnd = -0.5;
        const amrex::Real ubnd =  1.5;

        if((xi < lbnd or xi > ubnd) or (eta < lbnd or eta > ubnd) or (zeta < lbnd or zeta > ubnd)) {

            const std::string& myfile = "mapping";
            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << std::endl;

            amrex::PrintToFile(myfile, Print::AllProcs) << "  xi:"  << xi   << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << " eta:"  << eta  << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "zeta:"  << zeta << std::endl;

            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Particle Cell Index:  " << ip << " " << jp << " "<< kp << " " << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Stencil:  " <<  i << " " <<  j << " "<<  k << " " << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "REF:      " << di << " " << dj << " "<< dk << " " << std::endl;

            const Real adx = apx(ip+1,jp  ,kp  ) - apx(ip  ,jp  ,kp  );
            const Real ady = apy(ip  ,jp+1,kp  ) - apy(ip  ,jp  ,kp  );
            const Real adz = apz(ip  ,jp  ,kp+1) - apz(ip  ,jp  ,kp  );

            const Real apnorm = std::sqrt(adx*adx + ady*ady + adz*adz);

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE( apnorm > std::numeric_limits<amrex::Real>::epsilon(),
                "Invalid EB normal. Not sure how to interpolate!");

            // EB normal
            const amrex::Real apnorminv = 1. / apnorm;
            const amrex::RealVect normal = {adx*apnorminv, ady*apnorminv, adz*apnorminv};



            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Particle:                   "
            << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;

            amrex::PrintToFile(myfile, Print::AllProcs) << "Cell Centroid:              "
                << (ip + 0.5 + ccent(ip,jp,kp,0))*dx[0] << " "
                << (jp + 0.5 + ccent(ip,jp,kp,1))*dx[1] << " "
                << (kp + 0.5 + ccent(ip,jp,kp,2))*dx[2] << std::endl;

            amrex::PrintToFile(myfile, Print::AllProcs) << "Boundary Centroid (EB):     "
                << (ip + 0.5 + bcent(ip,jp,kp,0))*dx[0] << " "
                << (jp + 0.5 + bcent(ip,jp,kp,1))*dx[1] << " "
                << (kp + 0.5 + bcent(ip,jp,kp,2))*dx[2] << std::endl;

            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Boundary Normal:            "
                << normal[0] << " " << normal[1] << " " << normal[2] << std::endl;

            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) <<  "Boundary to Particle:      "
                << (pos[0] - (ip + 0.5 + bcent(ip,jp,kp,0))*dx[0]) << "  "
                << (pos[1] - (jp + 0.5 + bcent(ip,jp,kp,1))*dx[1]) << "  "
                << (pos[2] - (kp + 0.5 + bcent(ip,jp,kp,2))*dx[2]) << std::endl;

            amrex::PrintToFile(myfile, Print::AllProcs) <<  "Boundary to Cell Centriod: "
                << (ccent(ip,jp,kp,0) - bcent(ip,jp,kp,0))*dx[0] << "  "
                << (ccent(ip,jp,kp,1) - bcent(ip,jp,kp,1))*dx[1] << "  "
                << (ccent(ip,jp,kp,2) - bcent(ip,jp,kp,2))*dx[2] << std::endl;

            amrex::PrintToFile(myfile, Print::AllProcs) <<  "Cell Centroid to Particle: "
                << (pos[0] - (ip + 0.5 + ccent(ip,jp,kp,0))*dx[0]) << "  "
                << (pos[1] - (jp + 0.5 + ccent(ip,jp,kp,1))*dx[1]) << "  "
                << (pos[2] - (kp + 0.5 + ccent(ip,jp,kp,2))*dx[2]) << std::endl;

           amrex::PrintToFile(myfile, Print::AllProcs) << std::endl
                << "particle loc  dot boundary centroid:  " << par_dot_EB << std::endl
                << "cell centroid dot boundary centroid:  " << cent_dot_EB << std::endl;

            // Connectedness is relative to the cell containing the particle
            const int Node0_isNotConnected = (not flags(ip,jp,kp).isConnected(di-1,dj-1,dk-1));
            const int Node1_isNotConnected = (not flags(ip,jp,kp).isConnected(di  ,dj-1,dk-1));
            const int Node2_isNotConnected = (not flags(ip,jp,kp).isConnected(di  ,dj  ,dk-1));
            const int Node3_isNotConnected = (not flags(ip,jp,kp).isConnected(di-1,dj  ,dk-1));
            const int Node4_isNotConnected = (not flags(ip,jp,kp).isConnected(di-1,dj-1,dk  ));
            const int Node5_isNotConnected = (not flags(ip,jp,kp).isConnected(di  ,dj-1,dk  ));
            const int Node6_isNotConnected = (not flags(ip,jp,kp).isConnected(di  ,dj  ,dk  ));
            const int Node7_isNotConnected = (not flags(ip,jp,kp).isConnected(di-1,dj  ,dk  ));


            // The centroid sits closer to the boundary than the EB.
            if(par_dot_EB > cent_dot_EB){

                {
                    int new_di = (normal[0] < 0.) ? 0 : 1;
                    int new_dj = (normal[1] < 0.) ? 0 : 1;
                    int new_dk = (normal[2] < 0.) ? 0 : 1;

                    amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << "Option 1:"
                        << "NEW REF:     " <<    new_di << " " <<    new_dj << " "<<    new_dk << " " << std::endl
                        << "New Stencil: " << ip+new_di << " " << jp+new_dj << " "<< kp+new_dk << " " << std::endl;
                }

                {

                    int new_di = ((pos[0] - (ip + 0.5 + bcent(ip,jp,kp,0))*dx[0]) < 0.) ? 0 : 1;
                    int new_dj = ((pos[1] - (jp + 0.5 + bcent(ip,jp,kp,1))*dx[1]) < 0.) ? 0 : 1;
                    int new_dk = ((pos[2] - (kp + 0.5 + bcent(ip,jp,kp,2))*dx[2]) < 0.) ? 0 : 1;

                    amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << "Option 2:"
                        << "NEW REF:     " <<    new_di << " " <<    new_dj << " "<<    new_dk << " " << std::endl
                        << "New Stencil: " << ip+new_di << " " << jp+new_dj << " "<< kp+new_dk << " " << std::endl;
                }

            }


            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << "area fractions:  " << std::endl
                << "  x:  " << apx(ip  ,jp  ,kp  ) << "    " << apx(ip+1,jp  ,kp  ) << std::endl
                << "  y:  " << apy(ip  ,jp  ,kp  ) << "    " << apy(ip  ,jp+1,kp  ) << std::endl
                << "  z:  " << apz(ip  ,jp  ,kp  ) << "    " << apz(ip  ,jp  ,kp+1) << std::endl;


            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << std::endl;
            for(int n(0); n<8; n++){
                amrex::PrintToFile(myfile, Print::AllProcs) << " node " << n << ": "
                    << nodes[n][0] << " "
                    << nodes[n][1] << " "
                    << nodes[n][2] << std::endl;
            }



            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << std::endl;
            for(int n(0); n<8; n++){
                amrex::PrintToFile(myfile, Print::AllProcs) << " values " << n << ": "
                    << values[n][0] << " "
                    << values[n][1] << " "
                    << values[n][2] << " "
                    << values[n][3] << std::endl;
            }

            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << "Covered: " << covered << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << "Covered logic:" << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node0 " << Node0_isNotConnected <<  std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node1 " << Node1_isNotConnected <<  std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node2 " << Node2_isNotConnected <<  std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node3 " << Node3_isNotConnected <<  std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node4 " << Node4_isNotConnected <<  std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node5 " << Node5_isNotConnected <<  std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node6 " << Node6_isNotConnected <<  std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node7 " << Node7_isNotConnected <<  std::endl;

            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node 0:" << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Index:  " <<  i-1 << " " <<  j-1 << " "<<  k-1 << " " << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "REF:    " << di-1 << " " << dj-1 << " "<< dk-1 << " " << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is regular?  "  << flags(Node0).isRegular() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is svalued?  "  << flags(Node0).isSingleValued() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is covered?  "  << flags(Node0).isCovered() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Cell Centroids:   "
                            << (Node0[0] + 0.5 + ccent(Node0,0))*dx[0] << " "
                            << (Node0[1] + 0.5 + ccent(Node0,1))*dx[1] << " "
                            << (Node0[2] + 0.5 + ccent(Node0,2))*dx[2] << std::endl;
            if(flags(Node0).isSingleValued()){
                amrex::PrintToFile(myfile, Print::AllProcs) << "Bndry Centroids:  "
                            << (Node0[0] + 0.5 + bcent(Node0,0))*dx[0] << " "
                            << (Node0[1] + 0.5 + bcent(Node0,1))*dx[1] << " "
                            << (Node0[2] + 0.5 + bcent(Node0,2))*dx[2] << std::endl;

            }

            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node 1:" << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Index:  " <<  i   << " " <<  j-1 << " "<<  k-1 << " " << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "REF:    " << di   << " " << dj-1 << " "<< dk-1 << " " << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is regular?  "  << flags(Node1).isRegular() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is svalued?  "  << flags(Node1).isSingleValued() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is covered?  "  << flags(Node1).isCovered() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Cell Centroids:   "
                            << (Node1[0] + 0.5 + ccent(Node1,0))*dx[0] << " "
                            << (Node1[1] + 0.5 + ccent(Node1,1))*dx[1] << " "
                            << (Node1[2] + 0.5 + ccent(Node1,2))*dx[2] << std::endl;
            if(flags(Node1).isSingleValued()){
                amrex::PrintToFile(myfile, Print::AllProcs) << "Bndry Centroids:  "
                            << (Node1[0] + 0.5 + bcent(Node1,0))*dx[0] << " "
                            << (Node1[1] + 0.5 + bcent(Node1,1))*dx[1] << " "
                            << (Node1[2] + 0.5 + bcent(Node1,2))*dx[2] << std::endl;
            }

            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node 2:" << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Index:  " <<  i   << " " <<  j   << " "<<  k-1 << " " << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "REF:    " << di   << " " << dj   << " "<< dk-1 << " " << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is regular?  "  << flags(Node2).isRegular() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is svalued?  "  << flags(Node2).isSingleValued() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is covered?  "  << flags(Node2).isCovered() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Cell Centroids:   "
                            << (Node2[0] + 0.5 + ccent(Node2,0))*dx[0] << " "
                            << (Node2[1] + 0.5 + ccent(Node2,1))*dx[1] << " "
                            << (Node2[2] + 0.5 + ccent(Node2,2))*dx[2] << std::endl;
            if(flags(Node2).isSingleValued()){
                amrex::PrintToFile(myfile, Print::AllProcs) << "Bndry Centroids:  "
                            << (Node2[0] + 0.5 + bcent(Node2,0))*dx[0] << " "
                            << (Node2[1] + 0.5 + bcent(Node2,1))*dx[1] << " "
                            << (Node2[2] + 0.5 + bcent(Node2,2))*dx[2] << std::endl;
            }

            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node 3:" << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Index:  " <<  i-1 << " " <<  j   << " "<<  k-1 << " " << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "REF:    " << di-1 << " " << dj   << " "<< dk-1 << " " << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is regular?  "  << flags(Node3).isRegular() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is svalued?  "  << flags(Node3).isSingleValued() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is covered?  "  << flags(Node3).isCovered() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Cell Centroids:   "
                            << (Node3[0] + 0.5 + ccent(Node3,0))*dx[0] << " "
                            << (Node3[1] + 0.5 + ccent(Node3,1))*dx[1] << " "
                            << (Node3[2] + 0.5 + ccent(Node3,2))*dx[2] << std::endl;
            if(flags(Node3).isSingleValued()){
                amrex::PrintToFile(myfile, Print::AllProcs) << "Bndry Centroids:  "
                            << (Node3[0] + 0.5 + bcent(Node3,0))*dx[0] << " "
                            << (Node3[1] + 0.5 + bcent(Node3,1))*dx[1] << " "
                            << (Node3[2] + 0.5 + bcent(Node3,2))*dx[2] << std::endl;
            }

            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node 4:" << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Index:  " <<  i-1 << " " <<  j-1 << " "<<  k   << " " << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "REF:    " << di-1 << " " << dj-1 << " "<< dk   << " " << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is regular?  "  << flags(Node4).isRegular() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is svalued?  "  << flags(Node4).isSingleValued() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is covered?  "  << flags(Node4).isCovered() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Cell Centroids:   "
                            << (Node4[0] + 0.5 + ccent(Node4,0))*dx[0] << " "
                            << (Node4[1] + 0.5 + ccent(Node4,1))*dx[1] << " "
                            << (Node4[2] + 0.5 + ccent(Node4,2))*dx[2] << std::endl;
            if(flags(Node4).isSingleValued()){
                amrex::PrintToFile(myfile, Print::AllProcs) << "Bndry Centroids:  "
                            << (Node4[0] + 0.5 + bcent(Node4,0))*dx[0] << " "
                            << (Node4[1] + 0.5 + bcent(Node4,1))*dx[1] << " "
                            << (Node4[2] + 0.5 + bcent(Node4,2))*dx[2] << std::endl;
            }

            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node 5:" << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Index:  " <<  i   << " " <<  j-1 << " "<<  k   << " " << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "REF:    " << di   << " " << dj-1 << " "<< dk   << " " << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is regular?  "  << flags(Node5).isRegular() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is svalued?  "  << flags(Node5).isSingleValued() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is covered?  "  << flags(Node5).isCovered() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Cell Centroids:   "
                            << (Node5[0] + 0.5 + ccent(Node5,0))*dx[0] << " "
                            << (Node5[1] + 0.5 + ccent(Node5,1))*dx[1] << " "
                            << (Node5[2] + 0.5 + ccent(Node5,2))*dx[2] << std::endl;
            if(flags(Node5).isSingleValued()){
                amrex::PrintToFile(myfile, Print::AllProcs) << "Bndry Centroids:  "
                            << (Node5[0] + 0.5 + bcent(Node5,0))*dx[0] << " "
                            << (Node5[1] + 0.5 + bcent(Node5,1))*dx[1] << " "
                            << (Node5[2] + 0.5 + bcent(Node5,2))*dx[2] << std::endl;
            }

            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node 6:" << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Index:  " <<  i   << " " <<  j   << " "<<  k   << " " << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "REF:    " << di   << " " << dj   << " "<< dk   << " " << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is regular?  "  << flags(Node6).isRegular() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is svalued?  "  << flags(Node6).isSingleValued() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is covered?  "  << flags(Node6).isCovered() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Cell Centroids:   "
                            << (Node6[0] + 0.5 + ccent(Node6,0))*dx[0] << " "
                            << (Node6[1] + 0.5 + ccent(Node6,1))*dx[1] << " "
                            << (Node6[2] + 0.5 + ccent(Node6,2))*dx[2] << std::endl;
            if(flags(Node6).isSingleValued()){
                amrex::PrintToFile(myfile, Print::AllProcs) << "Bndry Centroids:  "
                            << (Node6[0] + 0.5 + bcent(Node6,0))*dx[0] << " "
                            << (Node6[1] + 0.5 + bcent(Node6,1))*dx[1] << " "
                            << (Node6[2] + 0.5 + bcent(Node6,2))*dx[2] << std::endl;
            }

            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node 7:" << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Index:  " <<  i-1 << " " <<  j   << " "<<  k   << " " << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "REF:    " << di-1 << " " << dj   << " "<< dk   << " " << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is regular?  "  << flags(Node7).isRegular() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is svalued?  "  << flags(Node7).isSingleValued() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node is covered?  "  << flags(Node7).isCovered() << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Cell Centroids:   "
                            << (Node7[0] + 0.5 + ccent(Node7,0))*dx[0] << " "
                            << (Node7[1] + 0.5 + ccent(Node7,1))*dx[1] << " "
                            << (Node7[2] + 0.5 + ccent(Node7,2))*dx[2] << std::endl;
            if(flags(Node7).isSingleValued()){
                amrex::PrintToFile(myfile, Print::AllProcs) << "Bndry Centroids:  "
                            << (Node7[0] + 0.5 + bcent(Node7,0))*dx[0] << " "
                            << (Node7[1] + 0.5 + bcent(Node7,1))*dx[1] << " "
                            << (Node7[2] + 0.5 + bcent(Node7,2))*dx[2] << std::endl;
            }

        }

        AMREX_ASSERT_WITH_MESSAGE( lbnd <= xi  &&  xi   <= ubnd,
            "Invalid interpolation: xi out of range!");
        AMREX_ASSERT_WITH_MESSAGE( lbnd <= eta &&  eta  <= ubnd,
            "Invalid interpolation: eta out of range!");
        AMREX_ASSERT_WITH_MESSAGE( lbnd <= zeta && zeta <= ubnd,
            "Invalid interpolation: zeta out of range!");
#endif
// DEBUG END   ////////////////////////////////////////////////////////////////

        for(int n=0; n < ncomp; n++) {
            interp_vals[n] = (1 - xi - eta - zeta + xi*eta + xi*zeta + eta*zeta - xi*eta*zeta) * values[0][n]
                           + (    xi              - xi*eta - xi*zeta            + xi*eta*zeta) * values[1][n]
                           + (                      xi*eta                      - xi*eta*zeta) * values[2][n]
                           + (         eta        - xi*eta           - eta*zeta + xi*eta*zeta) * values[3][n]
                           + (               zeta          - xi*zeta - eta*zeta + xi*eta*zeta) * values[4][n]
                           + (                               xi*zeta            - xi*eta*zeta) * values[5][n]
                           + (                                                    xi*eta*zeta) * values[6][n]
                           + (                                         eta*zeta - xi*eta*zeta) * values[7][n];
        }


        for(int it(0); it < 8; it++)
          delete values[it];
    }
}



#endif
