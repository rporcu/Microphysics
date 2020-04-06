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
               const amrex::Array4<const amrex::Real> ccent,
               const amrex::Array4<const amrex::Real> bcent,
               const amrex::Array4<const amrex::Real> apx,
               const amrex::Array4<const amrex::Real> apy,
               const amrex::Array4<const amrex::Real> apz,
               amrex::Array4<amrex::Real const> const& vel_g,
               amrex::Array4<amrex::Real const> const& ep_g,
               amrex::Real* vel_gp, amrex::Real& ep_gp)
{

    // Particle position relative to cell center [-0.5, 0.5]
    Real gx = pos[0]*dxi[0] - (ip + 0.5);
    Real gy = pos[1]*dxi[1] - (jp + 0.5);
    Real gz = pos[2]*dxi[2] - (kp + 0.5);

    // Use the centoid location of the cell containing the particle
    // to determine the interpolation stencil. If the particle is
    // on the low side, then the high side stencil is the index of the
    // cell contianing the particle, otherwise the particle is in
    // the low side cell.
    int i = (gx < ccent(ip,jp,kp,0)) ? ip : ip + 1;
    int j = (gy < ccent(ip,jp,kp,1)) ? jp : jp + 1;
    int k = (gz < ccent(ip,jp,kp,2)) ? kp : kp + 1;

    amrex::Real nodes[8][3];
    amrex::Real values[8][4];

    int di = i - ip; // -1 or 0
    int dj = j - jp; // -1 or 0
    int dk = k - kp; // -1 or 0

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


    amrex::Real par_dot_EB(1.);
    amrex::Real cent_dot_EB(1.);

    // If the particle is inside a cut-cell, verify that it is on the correct side of the EB
    // before trying to interpolate.
    if(flags(ip,jp,kp).isSingleValued()){

        const Real adx = apx(ip+1,jp  ,kp  ) - apx(ip  ,jp  ,kp  );
        const Real ady = apy(ip  ,jp+1,kp  ) - apy(ip  ,jp  ,kp  );
        const Real adz = apz(ip  ,jp  ,kp+1) - apz(ip  ,jp  ,kp  );

        const Real apnorm = std::sqrt(adx*adx + ady*ady + adz*adz);

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( apnorm > std::numeric_limits<amrex::Real>::epsilon(),
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

    }

    // A negative value implies that the particle is 'behind' the EB (and therefore inside
    // the wall). Although this shouldn't occur often, it could for fast moving particles.
    // Set the fluid velocity to zero and use the cut-cells volume fraction.
    if( par_dot_EB <= std::numeric_limits<amrex::Real>::epsilon()){

        // Void fraction is constant
        ep_gp = ep_g(ip,jp,kp);

        // Set the velocity to zero
        vel_gp[0] = 0.;
        vel_gp[1] = 0.;
        vel_gp[2] = 0.;

// DEBUG START ///////////////////////////////////////////////////////////////////////////////////////////////
#if(0)
            const std::string& myfile = "behind";

            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << "particle: "
                << pos[0] << " " << pos[1] << " " << pos[2] << std::endl << std::endl;

            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl  << "Cell Centroids:   "
                << (ip + 0.5 + ccent(ip,jp,kp,0))*dx[0] << " "
                << (jp + 0.5 + ccent(ip,jp,kp,1))*dx[1] << " "
                << (kp + 0.5 + ccent(ip,jp,kp,2))*dx[2] << std::endl;

            if(flags(ip,jp,kp).isSingleValued()){
                amrex::PrintToFile(myfile, Print::AllProcs) << "Bndry Centroids:  "
                    << (ip + 0.5 + bcent(ip,jp,kp,0))*dx[0] << " "
                    << (jp + 0.5 + bcent(ip,jp,kp,1))*dx[1] << " "
                    << (kp + 0.5 + bcent(ip,jp,kp,2))*dx[2] << std::endl;
            }
#endif
// DEBUG END   ///////////////////////////////////////////////////////////////////////////////////////////////

    // Four or more cells in the 2^3 stencil are covered implying that the particle
    //  is very close to the EB. We could do the gymnatics to figure out how to
    // interpolate from the few cells available, but we don't. We use the cell
    // containing the particle and ignore the rest.
    } else if(covered > 4) {

        // Linear interpolation from cell center to wall.
        const amrex::Real scale = (par_dot_EB / cent_dot_EB);
        vel_gp[0] = scale*vel_g(ip,jp,kp,0);
        vel_gp[1] = scale*vel_g(ip,jp,kp,1);
        vel_gp[2] = scale*vel_g(ip,jp,kp,2);


// DEBUG START ///////////////////////////////////////////////////////////////////////////////////////////////
#if(0)
        const std::string& myfile = "covered";

        amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << "particle: "
            << pos[0] << " " << pos[1] << " " << pos[2] << std::endl << std::endl;

        amrex::PrintToFile(myfile, Print::AllProcs) << std::endl  << "Cell Centroids:   "
            << (ip + 0.5 + ccent(ip,jp,kp,0))*dx[0] << " "
            << (jp + 0.5 + ccent(ip,jp,kp,1))*dx[1] << " "
            << (kp + 0.5 + ccent(ip,jp,kp,2))*dx[2] << std::endl;

        if(flags(ip,jp,kp).isSingleValued()){
            amrex::PrintToFile(myfile, Print::AllProcs) << "Bndry Centroids:  "
                << (ip + 0.5 + bcent(ip,jp,kp,0))*dx[0] << " "
                << (jp + 0.5 + bcent(ip,jp,kp,1))*dx[1] << " "
                << (kp + 0.5 + bcent(ip,jp,kp,2))*dx[2] << std::endl;
        }
        amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << "EB Normal:   "
            << normal[0] << " "
            << normal[1] << " "
            << normal[2] << std::endl;

        amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << "Result:" << std::endl
                            << "u_g  " << vel_gp[0] << "  "
                            << "v_g  " << vel_gp[1] << "  "
                            << "w_g  " << vel_gp[2] << "  " << std::endl
                            << "ep_g " << ep_gp << std::endl;

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( covered  <= 4, "Covered cells > 4");
#endif
// DEBUG END   ///////////////////////////////////////////////////////////////////////////////////////////////

    } else {

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

        // Connectedness is relative to the cell containing the particle
        const int Node0_isNotConnected = (not flags(ip,jp,kp).isConnected(di-1,dj-1,dk-1));
        const int Node1_isNotConnected = (not flags(ip,jp,kp).isConnected(di  ,dj-1,dk-1));
        const int Node2_isNotConnected = (not flags(ip,jp,kp).isConnected(di  ,dj  ,dk-1));
        const int Node3_isNotConnected = (not flags(ip,jp,kp).isConnected(di-1,dj  ,dk-1));
        const int Node4_isNotConnected = (not flags(ip,jp,kp).isConnected(di-1,dj-1,dk  ));
        const int Node5_isNotConnected = (not flags(ip,jp,kp).isConnected(di  ,dj-1,dk  ));
        const int Node6_isNotConnected = (not flags(ip,jp,kp).isConnected(di  ,dj  ,dk  ));
        const int Node7_isNotConnected = (not flags(ip,jp,kp).isConnected(di-1,dj  ,dk  ));


        /*---------------------------------------------------------------------------------*
        *                                                                                  *
        *                                    NODE 0                                        *
        *                                                                                  *
        *----------------------------------------------------------------------------------*/
        if(Node0_isNotConnected) {
            //amrex::Print() << std::endl << "Node 0 is not connected." << std::endl;
            amrex::IntVect NodeEB = Node0_EB(covered,
                                            Node0, Node0_isNotConnected,
                                            Node1, Node1_isNotConnected,
                                            Node2, Node2_isNotConnected,
                                            Node3, Node3_isNotConnected,
                                            Node4, Node4_isNotConnected,
                                            Node5, Node5_isNotConnected,
                                            Node6, Node6_isNotConnected,
                                            Node7, Node7_isNotConnected);


            nodes[0][0] = (NodeEB[0] + 0.5 + bcent(NodeEB,0))*dx[0];
            nodes[0][1] = (NodeEB[1] + 0.5 + bcent(NodeEB,1))*dx[1];
            nodes[0][2] = (NodeEB[2] + 0.5 + bcent(NodeEB,2))*dx[2];

            values[0][0] = 0.;
            values[0][1] = 0.;
            values[0][2] = 0.;
            values[0][3] = ep_g(NodeEB);

        } else {

            nodes[0][0] = (Node0[0] + 0.5 + ccent(Node0,0))*dx[0];
            nodes[0][1] = (Node0[1] + 0.5 + ccent(Node0,1))*dx[1];
            nodes[0][2] = (Node0[2] + 0.5 + ccent(Node0,2))*dx[2];

            values[0][0] = vel_g(Node0,0);
            values[0][1] = vel_g(Node0,1);
            values[0][2] = vel_g(Node0,2);
            values[0][3] =  ep_g(Node0);

        }

        /*---------------------------------------------------------------------------------*
        *                                                                                  *
        *                                    NODE 1                                        *
        *                                                                                  *
        *----------------------------------------------------------------------------------*/
        if(Node1_isNotConnected) {
            //amrex::Print() << std::endl << "Node is not connected." << std::endl;
            amrex::IntVect NodeEB = Node1_EB(covered,
                                            Node0, Node0_isNotConnected,
                                            Node1, Node1_isNotConnected,
                                            Node2, Node2_isNotConnected,
                                            Node3, Node3_isNotConnected,
                                            Node4, Node4_isNotConnected,
                                            Node5, Node5_isNotConnected,
                                            Node6, Node6_isNotConnected,
                                            Node7, Node7_isNotConnected);

            nodes[1][0] = (NodeEB[0] + 0.5 + bcent(NodeEB,0))*dx[0];
            nodes[1][1] = (NodeEB[1] + 0.5 + bcent(NodeEB,1))*dx[1];
            nodes[1][2] = (NodeEB[2] + 0.5 + bcent(NodeEB,2))*dx[2];

            values[1][0] = 0.;
            values[1][1] = 0.;
            values[1][2] = 0.;
            values[1][3] = ep_g(NodeEB);

        } else {

            nodes[1][0] = (Node1[0] + 0.5 + ccent(Node1,0))*dx[0];
            nodes[1][1] = (Node1[1] + 0.5 + ccent(Node1,1))*dx[1];
            nodes[1][2] = (Node1[2] + 0.5 + ccent(Node1,2))*dx[2];

            values[1][0] = vel_g(Node1,0);
            values[1][1] = vel_g(Node1,1);
            values[1][2] = vel_g(Node1,2);
            values[1][3] = ep_g(Node1);

        }

        /*---------------------------------------------------------------------------------*
        *                                                                                  *
        *                                    NODE 2                                        *
        *                                                                                  *
        *----------------------------------------------------------------------------------*/
        if(Node2_isNotConnected) {
            //amrex::Print() << std::endl << "Node 2 is not connected." << std::endl;
            amrex::IntVect NodeEB = Node2_EB(covered,
                                            Node0, Node0_isNotConnected,
                                            Node1, Node1_isNotConnected,
                                            Node2, Node2_isNotConnected,
                                            Node3, Node3_isNotConnected,
                                            Node4, Node4_isNotConnected,
                                            Node5, Node5_isNotConnected,
                                            Node6, Node6_isNotConnected,
                                            Node7, Node7_isNotConnected);

            nodes[2][0] = (NodeEB[0] + 0.5 + bcent(NodeEB,0))*dx[0];
            nodes[2][1] = (NodeEB[1] + 0.5 + bcent(NodeEB,1))*dx[1];
            nodes[2][2] = (NodeEB[2] + 0.5 + bcent(NodeEB,2))*dx[2];

            values[2][0] = 0.;
            values[2][1] = 0.;
            values[2][2] = 0.;
            values[2][3] = ep_g(NodeEB);

        } else {

            nodes[2][0] = (Node2[0] + 0.5 + ccent(Node2,0))*dx[0];
            nodes[2][1] = (Node2[1] + 0.5 + ccent(Node2,1))*dx[1];
            nodes[2][2] = (Node2[2] + 0.5 + ccent(Node2,2))*dx[2];

            values[2][0] = vel_g(Node2,0);
            values[2][1] = vel_g(Node2,1);
            values[2][2] = vel_g(Node2,2);
            values[2][3] =  ep_g(Node2);

        }

        /*---------------------------------------------------------------------------------*
        *                                                                                  *
        *                                    NODE 3                                        *
        *                                                                                  *
        *----------------------------------------------------------------------------------*/
        if(Node3_isNotConnected) {
            //amrex::Print() << std::endl << "Node 3 is not connected." << std::endl;
            amrex::IntVect NodeEB = Node3_EB(covered,
                                            Node0, Node0_isNotConnected,
                                            Node1, Node1_isNotConnected,
                                            Node2, Node2_isNotConnected,
                                            Node3, Node3_isNotConnected,
                                            Node4, Node4_isNotConnected,
                                            Node5, Node5_isNotConnected,
                                            Node6, Node6_isNotConnected,
                                            Node7, Node7_isNotConnected);

            nodes[3][0] = (NodeEB[0] + 0.5 + bcent(NodeEB,0))*dx[0];
            nodes[3][1] = (NodeEB[1] + 0.5 + bcent(NodeEB,1))*dx[1];
            nodes[3][2] = (NodeEB[2] + 0.5 + bcent(NodeEB,2))*dx[2];

            values[3][0] = 0.;
            values[3][1] = 0.;
            values[3][2] = 0.;
            values[3][3] = ep_g(NodeEB);

            } else {

            nodes[3][0] = (Node3[0] + 0.5 + ccent(Node3,0))*dx[0];
            nodes[3][1] = (Node3[1] + 0.5 + ccent(Node3,1))*dx[1];
            nodes[3][2] = (Node3[2] + 0.5 + ccent(Node3,2))*dx[2];

            values[3][0] = vel_g(Node3,0);
            values[3][1] = vel_g(Node3,1);
            values[3][2] = vel_g(Node3,2);
            values[3][3] =  ep_g(Node3);
        }

        /*---------------------------------------------------------------------------------*
        *                                                                                  *
        *                                    NODE 4                                        *
        *                                                                                  *
        *----------------------------------------------------------------------------------*/
        if(Node4_isNotConnected) {
            //amrex::Print() << std::endl << "Node 4 is not connected." << std::endl;
            amrex::IntVect NodeEB = Node4_EB(covered,
                                            Node0, Node0_isNotConnected,
                                            Node1, Node1_isNotConnected,
                                            Node2, Node2_isNotConnected,
                                            Node3, Node3_isNotConnected,
                                            Node4, Node4_isNotConnected,
                                            Node5, Node5_isNotConnected,
                                            Node6, Node6_isNotConnected,
                                            Node7, Node7_isNotConnected);

        nodes[4][0] = (NodeEB[0] + 0.5 + bcent(NodeEB,0))*dx[0];
        nodes[4][1] = (NodeEB[1] + 0.5 + bcent(NodeEB,1))*dx[1];
        nodes[4][2] = (NodeEB[2] + 0.5 + bcent(NodeEB,2))*dx[2];

        values[4][0] = 0.;
        values[4][1] = 0.;
        values[4][2] = 0.;
        values[4][3] = ep_g(NodeEB);

        } else {

        nodes[4][0] = (Node4[0] + 0.5 + ccent(Node4,0))*dx[0];
        nodes[4][1] = (Node4[1] + 0.5 + ccent(Node4,1))*dx[1];
        nodes[4][2] = (Node4[2] + 0.5 + ccent(Node4,2))*dx[2];

        values[4][0] = vel_g(Node4,0);
        values[4][1] = vel_g(Node4,1);
        values[4][2] = vel_g(Node4,2);
        values[4][3] =  ep_g(Node4);

        }

        /*---------------------------------------------------------------------------------*
        *                                                                                  *
        *                                    NODE 5                                        *
        *                                                                                  *
        *----------------------------------------------------------------------------------*/
        if(Node5_isNotConnected) {
            //amrex::Print() << std::endl << "Node 5 is not connected." << std::endl;
            amrex::IntVect NodeEB = Node5_EB(covered,
                                            Node0, Node0_isNotConnected,
                                            Node1, Node1_isNotConnected,
                                            Node2, Node2_isNotConnected,
                                            Node3, Node3_isNotConnected,
                                            Node4, Node4_isNotConnected,
                                            Node5, Node5_isNotConnected,
                                            Node6, Node6_isNotConnected,
                                            Node7, Node7_isNotConnected);

            nodes[5][0] = (NodeEB[0] + 0.5 + bcent(NodeEB,0))*dx[0];
            nodes[5][1] = (NodeEB[1] + 0.5 + bcent(NodeEB,1))*dx[1];
            nodes[5][2] = (NodeEB[2] + 0.5 + bcent(NodeEB,2))*dx[2];

            values[5][0] = 0.;
            values[5][1] = 0.;
            values[5][2] = 0.;
            values[5][3] = ep_g(NodeEB);

        } else {

            nodes[5][0] = (Node5[0] + 0.5 + ccent(Node5,0))*dx[0];
            nodes[5][1] = (Node5[1] + 0.5 + ccent(Node5,1))*dx[1];
            nodes[5][2] = (Node5[2] + 0.5 + ccent(Node5,2))*dx[2];

            values[5][0] = vel_g(Node5,0);
            values[5][1] = vel_g(Node5,1);
            values[5][2] = vel_g(Node5,2);
            values[5][3] =  ep_g(Node5);
        }

        /*---------------------------------------------------------------------------------*
        *                                                                                  *
        *                                    NODE 6                                        *
        *                                                                                  *
        *----------------------------------------------------------------------------------*/
        if(Node6_isNotConnected) {
            //amrex::Print() << std::endl << "Node 6 is not connected." << std::endl;
            amrex::IntVect NodeEB = Node6_EB(covered,
                                            Node0, Node0_isNotConnected,
                                            Node1, Node1_isNotConnected,
                                            Node2, Node2_isNotConnected,
                                            Node3, Node3_isNotConnected,
                                            Node4, Node4_isNotConnected,
                                            Node5, Node5_isNotConnected,
                                            Node6, Node6_isNotConnected,
                                            Node7, Node7_isNotConnected);

            nodes[6][0] = (NodeEB[0] + 0.5 + bcent(NodeEB,0))*dx[0];
            nodes[6][1] = (NodeEB[1] + 0.5 + bcent(NodeEB,1))*dx[1];
            nodes[6][2] = (NodeEB[2] + 0.5 + bcent(NodeEB,2))*dx[2];

            values[6][0] = 0.;
            values[6][1] = 0.;
            values[6][2] = 0.;
            values[6][3] = ep_g(NodeEB);

        } else {

            nodes[6][0] = (Node6[0] + 0.5 + ccent(Node6,0))*dx[0];
            nodes[6][1] = (Node6[1] + 0.5 + ccent(Node6,1))*dx[1];
            nodes[6][2] = (Node6[2] + 0.5 + ccent(Node6,2))*dx[2];

            values[6][0] = vel_g(Node6,0);
            values[6][1] = vel_g(Node6,1);
            values[6][2] = vel_g(Node6,2);
            values[6][3] =  ep_g(Node6);

        }

        /*---------------------------------------------------------------------------------*
        *                                                                                  *
        *                                    NODE 7                                        *
        *                                                                                  *
        *----------------------------------------------------------------------------------*/
        if(Node7_isNotConnected) {
            //amrex::Print() << std::endl << "Node 7 is not connected." << std::endl;
            amrex::IntVect NodeEB = Node7_EB(covered,
                                            Node0, Node0_isNotConnected,
                                            Node1, Node1_isNotConnected,
                                            Node2, Node2_isNotConnected,
                                            Node3, Node3_isNotConnected,
                                            Node4, Node4_isNotConnected,
                                            Node5, Node5_isNotConnected,
                                            Node6, Node6_isNotConnected,
                                            Node7, Node7_isNotConnected);

            nodes[7][0] = (NodeEB[0] + 0.5 + bcent(NodeEB,0))*dx[0];
            nodes[7][1] = (NodeEB[1] + 0.5 + bcent(NodeEB,1))*dx[1];
            nodes[7][2] = (NodeEB[2] + 0.5 + bcent(NodeEB,2))*dx[2];

            values[7][0] = 0.;
            values[7][1] = 0.;
            values[7][2] = 0.;
            values[7][3] = ep_g(NodeEB);

        } else {

            nodes[7][0] = (Node7[0] + 0.5 + ccent(Node7,0))*dx[0];
            nodes[7][1] = (Node7[1] + 0.5 + ccent(Node7,1))*dx[1];
            nodes[7][2] = (Node7[2] + 0.5 + ccent(Node7,2))*dx[2];

            values[7][0] = vel_g(Node7,0);
            values[7][1] = vel_g(Node7,1);
            values[7][2] = vel_g(Node7,2);
            values[7][3] =  ep_g(Node7);

        }





#if(0)
        amrex::Print() << std::endl << std::endl;
        for(int n(0); n<8; n++){
        amrex::Print() << " node " << n << ": "
        << nodes[n][0] << " "
        << nodes[n][1] << " "
        << nodes[n][2] << std::endl;
        }

        amrex::Print() << std::endl << std::endl;
        for(int n(0); n<8; n++){
        amrex::Print() << " values " << n << ": "
        << values[n][0] << " "
        << values[n][1] << " "
        << values[n][2] << " "
        << values[n][3] << std::endl;
        }
#endif

        amrex::Real xi;
        amrex::Real eta;
        amrex::Real zeta;

        get_interp_mapping(pos, nodes, dxi, xi, eta, zeta);




        if((xi < -0.1 or xi > 1.1) or (eta < -0.1 or eta > 1.1) or (zeta < -0.1 or zeta > 1.1)) {

// DEBUG START ///////////////////////////////////////////////////////////////////////////////////////////////
#if(0)
            const std::string& myfile = "mapping";
            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << std::endl;

            amrex::PrintToFile(myfile, Print::AllProcs) << "  xi:"  << xi   << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << " eta:"  << eta  << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "zeta:"  << zeta << std::endl;

            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << std::endl;
            for(int n(0); n<8; n++){
            amrex::PrintToFile(myfile, Print::AllProcs) << " node " << n << ": "
            << nodes[n][0] << " "
            << nodes[n][1] << " "
            << nodes[n][2] << std::endl;
            }

            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << "particle: "
            << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;


            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << std::endl;
            for(int n(0); n<8; n++){
            amrex::PrintToFile(myfile, Print::AllProcs) << " values " << n << ": "
            << values[n][0] << " "
            << values[n][1] << " "
            << values[n][2] << " "
            << values[n][3] << std::endl;
            }

            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << std::endl << "Covered: " << covered << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << std::endl << "Covered logic:" << std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node0 " << Node0_isNotConnected <<  std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node1 " << Node1_isNotConnected <<  std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node2 " << Node2_isNotConnected <<  std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node3 " << Node3_isNotConnected <<  std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node4 " << Node4_isNotConnected <<  std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node5 " << Node5_isNotConnected <<  std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node6 " << Node6_isNotConnected <<  std::endl;
            amrex::PrintToFile(myfile, Print::AllProcs) << "Node7 " << Node7_isNotConnected <<  std::endl;

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE( -0.1 <= xi  &&  xi   <= 1.1, "Invalid xi for interpolation");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE( -0.1 <= eta &&  eta  <= 1.1, "Invalid eta for interpolation");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE( -0.1 <= zeta && zeta <= 1.1, "Invalid zeta for interpolation");


#endif
// DEBUG END   ///////////////////////////////////////////////////////////////////////////////////////////////


        } else {

            amrex::Real val_gp[4];
            for(int lc(0); lc<4; lc++){
                val_gp[lc] = (1 - xi - eta - zeta + xi*eta + xi*zeta + eta*zeta - xi*eta*zeta) * values[0][lc]
                           + (    xi              - xi*eta - xi*zeta            + xi*eta*zeta) * values[1][lc]
                           + (                      xi*eta                      - xi*eta*zeta) * values[2][lc]
                           + (         eta        - xi*eta           - eta*zeta + xi*eta*zeta) * values[3][lc]
                           + (               zeta          - xi*zeta - eta*zeta + xi*eta*zeta) * values[4][lc]
                           + (                               xi*zeta            - xi*eta*zeta) * values[5][lc]
                           + (                                                    xi*eta*zeta) * values[6][lc]
                           + (                                         eta*zeta - xi*eta*zeta) * values[7][lc];
            }

            vel_gp[0] = val_gp[0];
            vel_gp[1] = val_gp[1];
            vel_gp[2] = val_gp[2];

            ep_gp = val_gp[3];
        }



// DEBUG START ///////////////////////////////////////////////////////////////////////////////////////////////
#if(0)
        {
            const std::string& myfile = "result";
            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << std::endl << "particle: "
                        << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;


            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << std::endl;
            for(int n(0); n<8; n++){
                amrex::PrintToFile(myfile, Print::AllProcs) << " values " << n << ": "
                                                            << values[n][0] << " "
                                                            << values[n][1] << " "
                                                            << values[n][2] << " "
                                                            << values[n][3] << std::endl;
            }
            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl;
            for(int n(0); n<8; n++){
                amrex::PrintToFile(myfile, Print::AllProcs) << " node " << n << ": "
                                                            << nodes[n][0] << " "
                                                            << nodes[n][1] << " "
                                                            << nodes[n][2] << std::endl;
            }

            amrex::PrintToFile(myfile, Print::AllProcs) << std::endl
                                << "u_g  " << vel_gp[0] << "  "
                                << "v_g  " << vel_gp[1] << "  "
                                << "w_g  " << vel_gp[2] << "  " << std::endl
                                << "ep_g " << ep_gp << std::endl;

        }
#endif
// DEBUG END   ///////////////////////////////////////////////////////////////////////////////////////////////

    }
}



#endif