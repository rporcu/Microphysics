#ifndef _MFIX_EB_INTERP_K_H_
#define _MFIX_EB_INTERP_K_H_

#include "MFIXParticleContainer.H"

/****************************************************************
 *                                                              *
 *                                                              *
 *                                                              *
 *                                                              *
 ***************************************************************/
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real f (const int dir,
               const amrex::RealVect& pos,
               const amrex::Real nodes[8][3],
               const amrex::Real& xi,
               const amrex::Real& eta,
               const amrex::Real& zeta)
{
    amrex::Real a0 = (nodes[0][dir] - pos[dir]);
    amrex::Real a1 = (nodes[1][dir] - nodes[0][dir]);
    amrex::Real a2 = (nodes[3][dir] - nodes[0][dir]);
    amrex::Real a3 = (nodes[4][dir] - nodes[0][dir]);
    amrex::Real a4 = (nodes[0][dir] - nodes[1][dir] + nodes[2][dir] - nodes[3][dir]);
    amrex::Real a5 = (nodes[0][dir] - nodes[1][dir] - nodes[4][dir] + nodes[5][dir]);
    amrex::Real a6 = (nodes[0][dir] - nodes[3][dir] - nodes[4][dir] + nodes[7][dir]);
    amrex::Real a7 = (nodes[1][dir] - nodes[2][dir] + nodes[3][dir] + nodes[4][dir] - nodes[5][dir] + nodes[6][dir] - nodes[7][dir] - nodes[0][dir]);

    return a0 + a1*xi + a2*eta + a3*zeta + a4*xi*eta + a5*xi*zeta + a6*eta*zeta + a7*xi*eta*zeta;
}



AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real dfdxi (const int dir,
                   const amrex::Real nodes[8][3],
                   const amrex::Real& xi,
                   const amrex::Real& eta,
                   const amrex::Real& zeta)
{
    amrex::Real a1 = (nodes[1][dir] - nodes[0][dir]);
    amrex::Real a4 = (nodes[0][dir] - nodes[1][dir] + nodes[2][dir] - nodes[3][dir]);
    amrex::Real a5 = (nodes[0][dir] - nodes[1][dir] - nodes[4][dir] + nodes[5][dir]);
    amrex::Real a7 = (nodes[1][dir] - nodes[2][dir] + nodes[3][dir] + nodes[4][dir] - nodes[5][dir] + nodes[6][dir] - nodes[7][dir] - nodes[0][dir]);

    return a1 + a4*eta + a5*zeta + a7*eta*zeta;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real dfdeta (const int dir,
                    const amrex::Real nodes[8][3],
                    const amrex::Real& xi,
                    const amrex::Real& eta,
                    const amrex::Real& zeta)
{

    amrex::Real a2 = (nodes[3][dir] - nodes[0][dir]);
    amrex::Real a4 = (nodes[0][dir] - nodes[1][dir] + nodes[2][dir] - nodes[3][dir]);
    amrex::Real a6 = (nodes[0][dir] - nodes[3][dir] - nodes[4][dir] + nodes[7][dir]);
    amrex::Real a7 = (nodes[1][dir] - nodes[2][dir] + nodes[3][dir] + nodes[4][dir] - nodes[5][dir] + nodes[6][dir] - nodes[7][dir] - nodes[0][dir]);

    return a2 + a4*xi + a6*zeta + a7*xi*zeta;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real dfdzeta (const int dir,
                     const amrex::Real nodes[8][3],
                     const amrex::Real& xi,
                     const amrex::Real& eta,
                     const amrex::Real& zeta)
{

    amrex::Real a3 = (nodes[4][dir] - nodes[0][dir]);
    amrex::Real a5 = (nodes[0][dir] - nodes[1][dir] - nodes[4][dir] + nodes[5][dir]);
    amrex::Real a6 = (nodes[0][dir] - nodes[3][dir] - nodes[4][dir] + nodes[7][dir]);
    amrex::Real a7 = (nodes[1][dir] - nodes[2][dir] + nodes[3][dir] + nodes[4][dir] - nodes[5][dir] + nodes[6][dir] - nodes[7][dir] - nodes[0][dir]);

    return a3 + a5*xi + a6*eta + a7*xi*eta;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_interp_mapping (const amrex::RealVect& pos,
                         const amrex::Real nodes[8][3],
                         const amrex::RealVect& dxi,
                         amrex::Real& xi,
                         amrex::Real& eta,
                         amrex::Real& zeta)
{
    xi   = (pos[0] - nodes[0][0])*dxi[0];
    eta  = (pos[1] - nodes[0][1])*dxi[1];
    zeta = (pos[2] - nodes[0][2])*dxi[2];

#if(0)
    amrex::Print() << 0
                   << "  xi: " << xi << "  "
                   << " eta: " << eta << "  "
                   << "zeta: " << zeta << std::endl;
#endif

    int lc(0);
    amrex::Real err(1.0);

    while(err > 1.0e-3 && lc < 10) {

        amrex::Real f0 = f(0, pos, nodes, xi, eta, zeta);
        amrex::Real f1 = f(1, pos, nodes, xi, eta, zeta);
        amrex::Real f2 = f(2, pos, nodes, xi, eta, zeta);

        amrex::Real df0dxi   = dfdxi(0, nodes, xi, eta, zeta);
        amrex::Real df0deta  = dfdeta(0, nodes, xi, eta, zeta);
        amrex::Real df0dzeta = dfdzeta(0, nodes, xi, eta, zeta);

        amrex::Real df1dxi   = dfdxi(1, nodes, xi, eta, zeta);
        amrex::Real df1deta  = dfdeta(1, nodes, xi, eta, zeta);
        amrex::Real df1dzeta = dfdzeta(1, nodes, xi, eta, zeta);

        amrex::Real df2dxi   = dfdxi(2, nodes, xi, eta, zeta);
        amrex::Real df2deta  = dfdeta(2, nodes, xi, eta, zeta);
        amrex::Real df2dzeta = dfdzeta(2, nodes, xi, eta, zeta);

        amrex::Real detJ = df0dxi   * (df1deta*df2dzeta - df1dzeta*df2deta)
                         - df0deta  * (df1dxi*df2dzeta  - df1dzeta*df2dxi)
                         + df0dzeta * (df1dxi*df2deta   - df1deta*df2dxi);

        amrex::Real detJ_xi = f0       * (df1deta*df2dzeta - df1dzeta*df2deta)
                            - df0deta  * (f1*df2dzeta      - df1dzeta*f2)
                            + df0dzeta * (f1*df2deta       - df1deta*f2);

        amrex::Real detJ_eta = df0dxi   * (f1*df2dzeta     - df1dzeta*f2)
                             - f0       * (df1dxi*df2dzeta - df1dzeta*df2dxi)
                             + df0dzeta * (df1dxi*f2       - f1*df2dxi);

        amrex::Real detJ_zeta = df0dxi   * (df1deta*f2     - f1*df2deta)
                              - df0deta  * (df1dxi*f2      - f1*df2dxi)
                              + f0       * (df1dxi*df2deta - df1deta*df2dxi);

#if(0)
        amrex::Print() << std::endl << std::endl;
        amrex::Print() << "detJ      " << detJ << std::endl;
        amrex::Print() << "detJ_xi   " << detJ_xi << std::endl;
        amrex::Print() << "detJ_eta  " << detJ_eta << std::endl;
        amrex::Print() << "detJ_zeta " << detJ_zeta << std::endl;
        amrex::Print() << std::endl << std::endl;
#endif
        amrex::Real new_xi   = xi   - detJ_xi   / detJ;
        amrex::Real new_eta  = eta  - detJ_eta  / detJ;
        amrex::Real new_zeta = zeta - detJ_zeta / detJ;

        lc += 1;
        err = amrex::max(std::abs(xi-new_xi), std::abs(eta-new_eta), std::abs(zeta-new_zeta));

        xi   = new_xi;
        eta  = new_eta;
        zeta = new_zeta;

#if(0)
    amrex::Print() << lc
                   << "  xi: " << xi << "  "
                   << " eta: " << eta << "  "
                   << "zeta: " << zeta << "  "
                   << " err: " << err << std::endl;
#endif
    }

}


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

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( covered  <= 4, "Too many convered cells");

    const amrex::IntVect Node0 = {i-1,j-1,k-1};
    const amrex::IntVect Node1 = {i  ,j-1,k-1};
    const amrex::IntVect Node2 = {i  ,j  ,k-1};
    const amrex::IntVect Node3 = {i-1,j  ,k-1};
    const amrex::IntVect Node4 = {i-1,j-1,k  };
    const amrex::IntVect Node5 = {i  ,j-1,k  };
    const amrex::IntVect Node6 = {i  ,j  ,k  };
    const amrex::IntVect Node7 = {i-1,j  ,k  };

    const int Node0_isNotConnected = (not flags(ip,jp,kp).isConnected(di-1,dj-1,dk-1));
    const int Node1_isNotConnected = (not flags(ip,jp,kp).isConnected(di  ,dj-1,dk-1));
    const int Node2_isNotConnected = (not flags(ip,jp,kp).isConnected(di  ,dj  ,dk-1));
    const int Node3_isNotConnected = (not flags(ip,jp,kp).isConnected(di-1,dj  ,dk-1));
    const int Node4_isNotConnected = (not flags(ip,jp,kp).isConnected(di-1,dj-1,dk  ));
    const int Node5_isNotConnected = (not flags(ip,jp,kp).isConnected(di  ,dj-1,dk  ));
    const int Node6_isNotConnected = (not flags(ip,jp,kp).isConnected(di  ,dj  ,dk  ));
    const int Node7_isNotConnected = (not flags(ip,jp,kp).isConnected(di-1,dj  ,dk  ));


    //amrex::Print() <<  std::endl <<  std::endl;
    //amrex::Print() << "Particle Index: " << ip << " " << jp << " "<< kp << std::endl;
    //amrex::Print() << "Stencil Index:  " << i  << " " << j  << " "<< k  << std::endl;
    //amrex::Print() << "Covered:        " << covered << std::endl;

    /*----------------------------------------------------------------------------------*
    *                                                                                  *
    *                                    NODE 0                                        *
    *                                                                                  *
    *----------------------------------------------------------------------------------*/
    if(Node0_isNotConnected) {
    //amrex::Print() << std::endl << "Node 0 is not connected." << std::endl;
    amrex::IntVect NodeEB;

    if(covered == 1) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        NodeEB = Node6;

    } else if(covered == 2) {
        if(Node1_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        NodeEB = Node7;
        } else if(Node3_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        NodeEB = Node5;
        } else if(Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        NodeEB = Node2;
        } else {
          amrex::Abort("Invalid interpolation for Node 0 with 2 covered cells.");
        }

    } else if(covered == 3){
        if(Node2_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        NodeEB = Node4;
        } else if(Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        NodeEB = Node1;
        } else if(Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        NodeEB = Node3;
        } else if(Node1_isNotConnected and Node3_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        NodeEB = Node2;
        } else if(Node1_isNotConnected and Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        NodeEB = Node5;
        } else if(Node3_isNotConnected and Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        NodeEB = Node7;
        } else {
          amrex::Abort("Invalid interpolation for Node 0 with 3 covered cells.");
        }

    } else if(covered == 4){
        if(Node1_isNotConnected and Node3_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        NodeEB = Node4;
        } else if(Node1_isNotConnected and Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        NodeEB = Node3;
        } else if(Node3_isNotConnected and Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        NodeEB = Node1;
        } else {
            amrex::Abort("Invalid interpolation for Node 0 with 4 covered cells.");
        }
    }

      nodes[0][0] = (NodeEB[0] + 0.5 + bcent(NodeEB,0))*dx[0];
      nodes[0][1] = (NodeEB[1] + 0.5 + bcent(NodeEB,1))*dx[1];
      nodes[0][2] = (NodeEB[2] + 0.5 + bcent(NodeEB,2))*dx[2];

      values[0][0] = 0.;
      values[0][1] = 0.;
      values[0][2] = 0.;
      values[0][3] = ep_g(NodeEB);

    } else {

      nodes[0][0] = (i - 0.5 + ccent(Node0,0))*dx[0];
      nodes[0][1] = (j - 0.5 + ccent(Node0,1))*dx[1];
      nodes[0][2] = (k - 0.5 + ccent(Node0,2))*dx[2];

      values[0][0] = vel_g(Node0,0);
      values[0][1] = vel_g(Node0,1);
      values[0][2] = vel_g(Node0,2);
      values[0][3] =  ep_g(Node0);

    }

    /*----------------------------------------------------------------------------------*
    *                                                                                  *
    *                                    NODE 1                                        *
    *                                                                                  *
    *----------------------------------------------------------------------------------*/
    if(Node1_isNotConnected) {
    //amrex::Print() << std::endl << "Node is not connected." << std::endl;
    amrex::IntVect NodeEB;

    if(covered == 1){
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        NodeEB = Node7;

    } else if(covered == 2) {
        if(Node0_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        NodeEB = Node7;
        } else if(Node2_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        NodeEB = Node7;
        } else if(Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        NodeEB = Node7;
        } else {
            amrex::Abort("Invalid interpolation for Node 1 with 2 covered cells.");
        }

    } else if(covered == 3){
        if(Node3_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        NodeEB = Node5;
        } else if(Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        NodeEB = Node2;
        } else if(Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        NodeEB = Node0;
        } else if(Node0_isNotConnected and Node2_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        NodeEB = Node3;
        } else if(Node0_isNotConnected and Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        NodeEB = Node4;
        } else if(Node2_isNotConnected and Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        NodeEB = Node6;
        } else {
            amrex::Abort("Invalid interpolation for Node 1 with 3 covered cells.");
        }

    } else if(covered == 4){
        if(Node0_isNotConnected and Node2_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        NodeEB = Node5;
        } else if(Node0_isNotConnected and Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        NodeEB = Node2;
        } else if(Node2_isNotConnected and Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        NodeEB = Node0;
        } else {
            amrex::Abort("Invalid interpolation for Node 1 with 4 covered cells.");
        }
    }

      nodes[1][0] = (NodeEB[0] + 0.5 + bcent(NodeEB,0))*dx[0];
      nodes[1][1] = (NodeEB[1] + 0.5 + bcent(NodeEB,1))*dx[1];
      nodes[1][2] = (NodeEB[2] + 0.5 + bcent(NodeEB,2))*dx[2];

      values[1][0] = 0.;
      values[1][1] = 0.;
      values[1][2] = 0.;
      values[1][3] = ep_g(NodeEB);

    } else {

      nodes[1][0] = (i + 0.5 + ccent(Node1,0))*dx[0];
      nodes[1][1] = (j - 0.5 + ccent(Node1,1))*dx[1];
      nodes[1][2] = (k - 0.5 + ccent(Node1,2))*dx[2];

      values[1][0] = vel_g(Node1,0);
      values[1][1] = vel_g(Node1,1);
      values[1][2] = vel_g(Node1,2);
      values[1][3] = ep_g(Node1);

    }

    /*----------------------------------------------------------------------------------*
    *                                                                                  *
    *                                    NODE 2                                        *
    *                                                                                  *
    *----------------------------------------------------------------------------------*/
    if(Node2_isNotConnected) {
    //amrex::Print() << std::endl << "Node 2 is not connected." << std::endl;
    amrex::IntVect NodeEB;

    if(covered == 1){
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        NodeEB = Node4;

    } else if(covered == 2) {
        if(Node1_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        NodeEB = Node7;
        } else if(Node3_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        NodeEB = Node5;
        } else if(Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        NodeEB = Node0;
        } else {
            amrex::Abort("Invalid interpolation for Node 2 with 2 covered cells.");
        }

    } else if(covered == 3){
        if(Node0_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        NodeEB = Node6;
        } else if(Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        NodeEB = Node3;
        } else if(Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        NodeEB = Node1;
        } else if(Node1_isNotConnected and Node3_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        NodeEB = Node0;
        } else if(Node1_isNotConnected and Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        NodeEB = Node5;
        } else if(Node3_isNotConnected and Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        NodeEB = Node7;
        } else {
            amrex::Abort("Invalid interpolation for Node 2 with 3 covered cells.");
        }

    } else if(covered == 4){
        if(Node1_isNotConnected and Node3_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        NodeEB = Node6;
        } else if(Node1_isNotConnected and Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        NodeEB = Node3;
        } else if(Node3_isNotConnected and Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        NodeEB = Node1;
        } else {
            amrex::Abort("Invalid interpolation for Node 2 with 4 covered cells.");
        }
    }

      nodes[2][0] = (NodeEB[0] + 0.5 + bcent(NodeEB,0))*dx[0];
      nodes[2][1] = (NodeEB[1] + 0.5 + bcent(NodeEB,1))*dx[1];
      nodes[2][2] = (NodeEB[2] + 0.5 + bcent(NodeEB,2))*dx[2];

      values[2][0] = 0.;
      values[2][1] = 0.;
      values[2][2] = 0.;
      values[2][3] = ep_g(NodeEB);

    } else {

      nodes[2][0] = (i + 0.5 + ccent(Node2,0))*dx[0];
      nodes[2][1] = (j + 0.5 + ccent(Node2,1))*dx[1];
      nodes[2][2] = (k - 0.5 + ccent(Node2,2))*dx[2];

      values[2][0] = vel_g(Node2,0);
      values[2][1] = vel_g(Node2,1);
      values[2][2] = vel_g(Node2,2);
      values[2][3] =  ep_g(Node2);

    }

    /*----------------------------------------------------------------------------------*
    *                                                                                  *
    *                                    NODE 3                                        *
    *                                                                                  *
    *----------------------------------------------------------------------------------*/
    if(Node3_isNotConnected) {
    //amrex::Print() << std::endl << "Node 3 is not connected." << std::endl;
    amrex::IntVect NodeEB;

    if(covered == 1){
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        NodeEB = Node5;

    } else if(covered == 2) {
        if(Node0_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        NodeEB = Node6;
        } else if(Node2_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        NodeEB = Node4;
        } else if(Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        NodeEB = Node1;
        } else {
            amrex::Abort("Invalid interpolation for Node 3 with 2 covered cells.");
        }

    } else if(covered == 3){
        if(Node1_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        NodeEB = Node7;
        } else if(Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        NodeEB = Node2;
        } else if(Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        NodeEB = Node0;
        } else if(Node0_isNotConnected and Node2_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        NodeEB = Node1;
        } else if(Node0_isNotConnected and Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        NodeEB = Node4;
        } else if(Node2_isNotConnected and Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        NodeEB = Node6;
        } else {
            amrex::Abort("Invalid interpolation for Node 3 with 3 covered cells.");
        }

    } else if(covered == 4){
        if(Node0_isNotConnected and Node2_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        NodeEB = Node7;
        } else if(Node0_isNotConnected and Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        NodeEB = Node2;
        } else if(Node2_isNotConnected and Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        NodeEB = Node0;
        } else {
            amrex::Abort("Invalid interpolation for Node 3 with 4 covered cells.");
        }
    }

      nodes[3][0] = (NodeEB[0] + 0.5 + bcent(NodeEB,0))*dx[0];
      nodes[3][1] = (NodeEB[1] + 0.5 + bcent(NodeEB,1))*dx[1];
      nodes[3][2] = (NodeEB[2] + 0.5 + bcent(NodeEB,2))*dx[2];

      values[3][0] = 0.;
      values[3][1] = 0.;
      values[3][2] = 0.;
      values[3][3] = ep_g(NodeEB);

    } else {

      nodes[3][0] = (i - 0.5 + ccent(Node3,0))*dx[0];
      nodes[3][1] = (j + 0.5 + ccent(Node3,1))*dx[1];
      nodes[3][2] = (k - 0.5 + ccent(Node3,2))*dx[2];

      values[3][0] = vel_g(Node3,0);
      values[3][1] = vel_g(Node3,1);
      values[3][2] = vel_g(Node3,2);
      values[3][3] =  ep_g(Node3);
    }

    /*----------------------------------------------------------------------------------*
    *                                                                                  *
    *                                    NODE 4                                        *
    *                                                                                  *
    *----------------------------------------------------------------------------------*/
    if(Node4_isNotConnected) {
    //amrex::Print() << std::endl << "Node 4 is not connected." << std::endl;
    amrex::IntVect NodeEB;

    if(covered == 1){
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        NodeEB = Node2;

    } else if(covered == 2) {
        if(Node0_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        NodeEB = Node6;
        } else if(Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        NodeEB = Node3;
        } else if(Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        NodeEB = Node1;
        } else {
            amrex::Abort("Invalid interpolation for Node 4 with 2 covered cells.");
        }

    } else if(covered == 3){
        if(Node3_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        NodeEB = Node5;
        } else if(Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        NodeEB = Node0;
        } else if(Node1_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        NodeEB = Node7;
        } else if(Node0_isNotConnected and Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        NodeEB = Node1;
        } else if(Node0_isNotConnected and Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        NodeEB = Node3;
        } else if(Node5_isNotConnected and Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        NodeEB = Node6;
        } else {
            amrex::Abort("Invalid interpolation for Node 4 with 3 covered cells.");
        }

    } else if(covered == 4){
        if(Node0_isNotConnected and Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        NodeEB = Node7;
        } else if(Node0_isNotConnected and Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        NodeEB = Node5;
        } else if(Node5_isNotConnected and Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        NodeEB = Node0;
        } else {
            amrex::Abort("Invalid interpolation for Node 4 with 4 covered cells.");
        }
    }

      nodes[4][0] = (NodeEB[0] + 0.5 + bcent(NodeEB,0))*dx[0];
      nodes[4][1] = (NodeEB[1] + 0.5 + bcent(NodeEB,1))*dx[1];
      nodes[4][2] = (NodeEB[2] + 0.5 + bcent(NodeEB,2))*dx[2];

      values[4][0] = 0.;
      values[4][1] = 0.;
      values[4][2] = 0.;
      values[4][3] = ep_g(NodeEB);

    } else {

      nodes[4][0] = (i - 0.5 + ccent(Node4,0))*dx[0];
      nodes[4][1] = (j - 0.5 + ccent(Node4,1))*dx[1];
      nodes[4][2] = (k + 0.5 + ccent(Node4,2))*dx[2];

      values[4][0] = vel_g(Node4,0);
      values[4][1] = vel_g(Node4,1);
      values[4][2] = vel_g(Node4,2);
      values[4][3] =  ep_g(Node4);

    }

    /*----------------------------------------------------------------------------------*
    *                                                                                  *
    *                                    NODE 5                                        *
    *                                                                                  *
    *----------------------------------------------------------------------------------*/
    if(Node5_isNotConnected) {
    //amrex::Print() << std::endl << "Node 5 is not connected." << std::endl;
    amrex::IntVect NodeEB;

    if(covered == 1) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        NodeEB = Node3;

    } else if(covered == 2) {
        if(Node1_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        NodeEB = Node7;
        } else if(Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        NodeEB = Node2;
        } else if(Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        NodeEB = Node0;
        } else {
            amrex::Abort("Invalid interpolation for Node 5 with 2 covered cells.");
        }

    } else if(covered == 3){
        if(Node0_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        NodeEB = Node6;
        } else if(Node2_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        NodeEB = Node4;
        } else if(Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        NodeEB = Node1;
        } else if(Node1_isNotConnected and Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        NodeEB = Node0;
        } else if(Node1_isNotConnected and Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        NodeEB = Node2;
        } else if(Node4_isNotConnected and Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        NodeEB = Node7;
        } else {
            amrex::Abort("Invalid interpolation for Node 5 with 3 covered cells.");
        }

    } else if(covered == 4){
        if(Node1_isNotConnected and Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        NodeEB = Node6;
        } else if(Node1_isNotConnected and Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        NodeEB = Node4;
        } else if(Node4_isNotConnected and Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        NodeEB = Node1;
        } else {
            amrex::Abort("Invalid interpolation for Node 5 with 4 covered cells.");
        }
    }

      nodes[5][0] = (NodeEB[0] + 0.5 + bcent(NodeEB,0))*dx[0];
      nodes[5][1] = (NodeEB[1] + 0.5 + bcent(NodeEB,1))*dx[1];
      nodes[5][2] = (NodeEB[2] + 0.5 + bcent(NodeEB,2))*dx[2];

      values[5][0] = 0.;
      values[5][1] = 0.;
      values[5][2] = 0.;
      values[5][3] = ep_g(NodeEB);

    } else {

      nodes[5][0] = (i + 0.5 + ccent(Node5,0))*dx[0];
      nodes[5][1] = (j - 0.5 + ccent(Node5,1))*dx[1];
      nodes[5][2] = (k + 0.5 + ccent(Node5,2))*dx[2];

      values[5][0] = vel_g(Node5,0);
      values[5][1] = vel_g(Node5,1);
      values[5][2] = vel_g(Node5,2);
      values[5][3] =  ep_g(Node5);
    }

    /*----------------------------------------------------------------------------------*
    *                                                                                  *
    *                                    NODE 6                                        *
    *                                                                                  *
    *----------------------------------------------------------------------------------*/
    if(Node6_isNotConnected) {
    //amrex::Print() << std::endl << "Node 6 is not connected." << std::endl;
    amrex::IntVect NodeEB;

    if(covered == 1){
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        NodeEB = Node0;

    } else if(covered == 2) {
        if(Node2_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        NodeEB = Node4;
        } else if(Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        NodeEB = Node3;
        } else if(Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        NodeEB = Node1;
        } else {
            amrex::Abort("Invalid interpolation for Node 6 with 2 covered cells.");
        }

    } else if(covered == 3){
        if(Node1_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        NodeEB = Node7;
        } else if(Node3_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        NodeEB = Node5;
        } else if(Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        NodeEB = Node2;
        } else if(Node2_isNotConnected and Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        NodeEB = Node1;
        } else if(Node2_isNotConnected and Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        NodeEB = Node3;
        } else if(Node5_isNotConnected and Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        NodeEB = Node4;
        } else {
            amrex::Abort("Invalid interpolation for Node 6 with 3 covered cells.");
        }

    } else if(covered == 4){
        if(Node2_isNotConnected and Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        NodeEB = Node7;
        } else if(Node2_isNotConnected and Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        NodeEB = Node5;
        } else if(Node5_isNotConnected and Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        NodeEB = Node2;
        } else {
            amrex::Abort("Invalid interpolation for Node 6 with 4 covered cells.");
        }
    }

      nodes[6][0] = (NodeEB[0] + 0.5 + bcent(NodeEB,0))*dx[0];
      nodes[6][1] = (NodeEB[1] + 0.5 + bcent(NodeEB,1))*dx[1];
      nodes[6][2] = (NodeEB[2] + 0.5 + bcent(NodeEB,2))*dx[2];

      values[6][0] = 0.;
      values[6][1] = 0.;
      values[6][2] = 0.;
      values[6][3] = ep_g(NodeEB);

    } else {

      nodes[6][0] = (i + 0.5 + ccent(Node6,0))*dx[0];
      nodes[6][1] = (j + 0.5 + ccent(Node6,1))*dx[1];
      nodes[6][2] = (k + 0.5 + ccent(Node6,2))*dx[2];

      values[6][0] = vel_g(Node6,0);
      values[6][1] = vel_g(Node6,1);
      values[6][2] = vel_g(Node6,2);
      values[6][3] =  ep_g(Node6);

    }

    /*----------------------------------------------------------------------------------*
    *                                                                                  *
    *                                    NODE 7                                        *
    *                                                                                  *
    *----------------------------------------------------------------------------------*/
    if(Node7_isNotConnected) {

    //amrex::Print() << std::endl << "Node 7 is not connected." << std::endl;
    amrex::IntVect NodeEB;

    if(covered == 1) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        NodeEB = Node1;

    } else if(covered == 2) {
        if(Node3_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        NodeEB = Node5;
        } else if(Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        NodeEB = Node2;
        } else if(Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        NodeEB = Node0;
        } else {
            amrex::Abort("Invalid interpolation for Node 7 with 2 covered cells.");
        }

    } else if(covered == 3){
        if(Node0_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        NodeEB = Node6;
        } else if(Node2_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        NodeEB = Node4;
        } else if(Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        NodeEB = Node3;
        } else if(Node3_isNotConnected and Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        NodeEB = Node0;
        } else if(Node3_isNotConnected and Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        NodeEB = Node2;
        } else if(Node4_isNotConnected and Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        NodeEB = Node5;
        } else {
            amrex::Abort("Invalid interpolation for Node 7 with 3 covered cells.");
        }

    } else if(covered == 4){
        if(Node3_isNotConnected and Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        NodeEB = Node6;
        } else if(Node3_isNotConnected and Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        NodeEB = Node4;
        } else if(Node4_isNotConnected and Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        NodeEB = Node3;
        } else {
            amrex::Abort("Invalid interpolation for Node 7 with 4 covered cells.");
        }
    }

      nodes[7][0] = (NodeEB[0] + 0.5 + bcent(NodeEB,0))*dx[0];
      nodes[7][1] = (NodeEB[1] + 0.5 + bcent(NodeEB,1))*dx[1];
      nodes[7][2] = (NodeEB[2] + 0.5 + bcent(NodeEB,2))*dx[2];

      values[7][0] = 0.;
      values[7][1] = 0.;
      values[7][2] = 0.;
      values[7][3] = ep_g(NodeEB);

    } else {

      nodes[7][0] = (i - 0.5 + ccent(Node7,0))*dx[0];
      nodes[7][1] = (j + 0.5 + ccent(Node7,1))*dx[1];
      nodes[7][2] = (k + 0.5 + ccent(Node7,2))*dx[2];

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






/////////////////////////////////////////////////////////////////////////////////////////////////
            if(xi < -0.1 or xi > 1.1) {
                const std::string& myfile = "xi";

                amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << std::endl;

                amrex::PrintToFile(myfile, Print::AllProcs) << "  xi:"  << xi   << std::endl;
                amrex::PrintToFile(myfile, Print::AllProcs) << " eta:"  << eta  << std::endl;
                amrex::PrintToFile(myfile, Print::AllProcs) << "zeta:"  << zeta << std::endl;

                amrex::PrintToFile(myfile, Print::AllProcs) << std::endl << std::endl;
                for(int n(0); n<8; n++){
                amrex::PrintToFile(myfile, Print::AllProcs) << "  node " << n << ": "
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
            }

            if(eta < -0.1 or eta > 1.1) {
                const std::string& myfile = "eta";
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
            }

            if(zeta < -1.0 or zeta > 1.1) {
                const std::string& myfile = "zeta";
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
            }


/////////////////////////////////////////////////////////////////////////////////////////////////


    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( -0.1 <= xi  &&  xi   <= 1.1, "Invalid xi for interpolation");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( -0.1 <= eta &&  eta  <= 1.1, "Invalid eta for interpolation");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( -0.1 <= zeta && zeta <= 1.1, "Invalid zeta for interpolation");

    //amrex::Print() << std::endl << std::endl << "Interpolated values:" << std::endl;
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
    //amrex::Print() << "value " << lc << " = " << val_gp[lc] << std::endl;
    }
    //amrex::Print() << std::endl << std::endl;

    vel_gp[0] = val_gp[0];
    vel_gp[1] = val_gp[1];
    vel_gp[2] = val_gp[2];

    ep_gp = val_gp[3];
}



#endif