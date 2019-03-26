#include "MFIXParticleContainer.H"

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void cross_product (const Real* a, const Real* b, Real* c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real dot_product (const Real* a, const Real* b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void cfrelvel (const MFIXParticleContainer::ParticleType& p1, 
               const MFIXParticleContainer::ParticleType& p2,
               Real& vrn, Real* vslip, const Real* normal, 
               const Real dist_mag)
{
    // translational relative velocity
    Real vreltrans[3];
    vreltrans[0] =  p1.rdata(realData::velx) - p2.rdata(realData::velx);
    vreltrans[1] =  p1.rdata(realData::vely) - p2.rdata(realData::vely);
    vreltrans[2] =  p1.rdata(realData::velz) - p2.rdata(realData::velz);

    // distance from the contact point to the particle centers
    Real dist_c1 = (dist_mag*dist_mag 
                    + p1.rdata(realData::radius)*p1.rdata(realData::radius) 
                    - p2.rdata(realData::radius)*p2.rdata(realData::radius)) / (2.0*dist_mag);

    Real dist_c2 = dist_mag - dist_c1;

    Real omega_sum[3];
    omega_sum[0] = p1.rdata(realData::omegax) * dist_c1 + p2.rdata(realData::omegax) * dist_c2;
    omega_sum[1] = p1.rdata(realData::omegay) * dist_c1 + p2.rdata(realData::omegay) * dist_c2;
    omega_sum[2] = p1.rdata(realData::omegaz) * dist_c1 + p2.rdata(realData::omegaz) * dist_c2;

    // the rotational relative velocity
    Real v_rot[3];
    cross_product(omega_sum, normal, v_rot);

    // total relative velocity
    vreltrans[0] += v_rot[0];
    vreltrans[1] += v_rot[1];
    vreltrans[2] += v_rot[2];

    // normal component of relative velocity
    vrn = dot_product(vreltrans, normal);

    vslip[0] = vreltrans[0] - vrn*normal[0];
    vslip[1] = vreltrans[1] - vrn*normal[1];
    vslip[2] = vreltrans[2] - vrn*normal[2];    
}
