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

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void interp_level_set (MFIXParticleContainer::ParticleType const& p,
                       amrex::Real& ls_value, const int n_refine,
                       amrex::Array4<amrex::Real const> const& phi,
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
{

    amrex::Real scaled_dxi[3];
    scaled_dxi[0] = n_refine * dxi[0];
    scaled_dxi[1] = n_refine * dxi[1];
    scaled_dxi[2] = n_refine * dxi[2];
    
    amrex::Real x = (p.pos(0) - plo[0]) * scaled_dxi[0];
    amrex::Real y = (p.pos(1) - plo[1]) * scaled_dxi[1];
    amrex::Real z = (p.pos(2) - plo[2]) * scaled_dxi[2];

    int i = std::floor(x);
    int j = std::floor(y);
    int k = std::floor(z);
    
    amrex::Real wx_hi = x - i;
    amrex::Real wy_hi = y - j;
    amrex::Real wz_hi = z - k;
    
    amrex::Real wx_lo = 1.0 - wx_hi;
    amrex::Real wy_lo = 1.0 - wy_hi;
    amrex::Real wz_lo = 1.0 - wz_hi;

    ls_value = phi(i,   j  , k  ) * wx_lo * wy_lo * wz_lo;
    ls_value = phi(i+1, j  , k  ) * wx_hi * wy_lo * wz_lo;
    ls_value = phi(i,   j+1, k  ) * wx_lo * wy_hi * wz_lo;
    ls_value = phi(i+1, j+1, k  ) * wx_hi * wy_hi * wz_lo;
    ls_value = phi(i,   j  , k+1) * wx_lo * wy_lo * wz_hi;
    ls_value = phi(i+1, j  , k+1) * wx_hi * wy_lo * wz_hi;
    ls_value = phi(i  , j+1, k+1) * wx_lo * wy_hi * wz_hi;
    ls_value = phi(i+1, j+1, k+1) * wx_hi * wy_hi * wz_hi;
}
