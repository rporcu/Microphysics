#include "MFIXParticleContainer.H"

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void cross_product (const amrex::Real* a, const amrex::Real* b, amrex::Real* c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real dot_product (const amrex::Real* a, const amrex::Real* b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void cfrelvel (const MFIXParticleContainer::ParticleType& p1, 
               const MFIXParticleContainer::ParticleType& p2,
               amrex::Real& vrn, amrex::Real* vslip,
               const amrex::Real* normal, 
               const amrex::Real dist_mag)
{
    // translational relative velocity
    amrex::Real vreltrans[3];
    vreltrans[0] =  p1.rdata(realData::velx) - p2.rdata(realData::velx);
    vreltrans[1] =  p1.rdata(realData::vely) - p2.rdata(realData::vely);
    vreltrans[2] =  p1.rdata(realData::velz) - p2.rdata(realData::velz);

    // distance from the contact point to the particle centers
    amrex::Real dist_c1 = (dist_mag*dist_mag 
                    + p1.rdata(realData::radius)*p1.rdata(realData::radius) 
                    - p2.rdata(realData::radius)*p2.rdata(realData::radius)) / (2.0*dist_mag);

    amrex::Real dist_c2 = dist_mag - dist_c1;

    amrex::Real omega_sum[3];
    omega_sum[0] = p1.rdata(realData::omegax) * dist_c1 + p2.rdata(realData::omegax) * dist_c2;
    omega_sum[1] = p1.rdata(realData::omegay) * dist_c1 + p2.rdata(realData::omegay) * dist_c2;
    omega_sum[2] = p1.rdata(realData::omegaz) * dist_c1 + p2.rdata(realData::omegaz) * dist_c2;

    // the rotational relative velocity
    amrex::Real v_rot[3];
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
amrex::Real interp_level_set (MFIXParticleContainer::ParticleType const& p,
                              const int n_refine,
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

    int i = static_cast<int>(amrex::Math::floor(x));
    int j = static_cast<int>(amrex::Math::floor(y));
    int k = static_cast<int>(amrex::Math::floor(z));
    
    amrex::Real wx_hi = x - i;
    amrex::Real wy_hi = y - j;
    amrex::Real wz_hi = z - k;
    
    amrex::Real wx_lo = 1.0 - wx_hi;
    amrex::Real wy_lo = 1.0 - wy_hi;
    amrex::Real wz_lo = 1.0 - wz_hi;

    amrex::Real ls_value = 0;
    ls_value += phi(i,   j  , k  ) * wx_lo * wy_lo * wz_lo;
    ls_value += phi(i+1, j  , k  ) * wx_hi * wy_lo * wz_lo;
    ls_value += phi(i,   j+1, k  ) * wx_lo * wy_hi * wz_lo;
    ls_value += phi(i+1, j+1, k  ) * wx_hi * wy_hi * wz_lo;
    ls_value += phi(i,   j  , k+1) * wx_lo * wy_lo * wz_hi;
    ls_value += phi(i+1, j  , k+1) * wx_hi * wy_lo * wz_hi;
    ls_value += phi(i  , j+1, k+1) * wx_lo * wy_hi * wz_hi;
    ls_value += phi(i+1, j+1, k+1) * wx_hi * wy_hi * wz_hi;
    
    return ls_value;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void level_set_normal (MFIXParticleContainer::ParticleType const& p,
                       const int n_refine,
                       amrex::Real* normal,
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

    int i = static_cast<int>(amrex::Math::floor(x));
    int j = static_cast<int>(amrex::Math::floor(y));
    int k = static_cast<int>(amrex::Math::floor(z));
    
    amrex::Real wx_hi = x - i;
    amrex::Real wy_hi = y - j;
    amrex::Real wz_hi = z - k;
    
    amrex::Real wx_lo = 1.0 - wx_hi;
    amrex::Real wy_lo = 1.0 - wy_hi;
    amrex::Real wz_lo = 1.0 - wz_hi;

    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 0.0;

    normal[0] -= phi(i,   j  , k  ) * scaled_dxi[0] * wy_lo * wz_lo;
    normal[0] += phi(i+1, j  , k  ) * scaled_dxi[0] * wy_lo * wz_lo;
    normal[0] -= phi(i,   j+1, k  ) * scaled_dxi[0] * wy_hi * wz_lo;
    normal[0] += phi(i+1, j+1, k  ) * scaled_dxi[0] * wy_hi * wz_lo;
    normal[0] -= phi(i,   j  , k+1) * scaled_dxi[0] * wy_lo * wz_hi;
    normal[0] += phi(i+1, j  , k+1) * scaled_dxi[0] * wy_lo * wz_hi;
    normal[0] -= phi(i  , j+1, k+1) * scaled_dxi[0] * wy_hi * wz_hi;
    normal[0] += phi(i+1, j+1, k+1) * scaled_dxi[0] * wy_hi * wz_hi;

    normal[1] -= phi(i,   j  , k  ) * scaled_dxi[1] * wx_lo * wz_lo;
    normal[1] += phi(i  , j+1, k  ) * scaled_dxi[1] * wx_lo * wz_lo;
    normal[1] -= phi(i+1, j  , k  ) * scaled_dxi[1] * wx_hi * wz_lo;
    normal[1] += phi(i+1, j+1, k  ) * scaled_dxi[1] * wx_hi * wz_lo;
    normal[1] -= phi(i,   j  , k+1) * scaled_dxi[1] * wx_lo * wz_hi;
    normal[1] += phi(i  , j+1, k+1) * scaled_dxi[1] * wx_lo * wz_hi;
    normal[1] -= phi(i+1, j  , k+1) * scaled_dxi[1] * wx_hi * wz_hi;
    normal[1] += phi(i+1, j+1, k+1) * scaled_dxi[1] * wx_hi * wz_hi;

    normal[2] -= phi(i  , j  , k  ) * scaled_dxi[2] * wx_lo * wy_lo;
    normal[2] += phi(i  , j  , k+1) * scaled_dxi[2] * wx_lo * wy_lo;
    normal[2] -= phi(i+1, j  , k  ) * scaled_dxi[2] * wx_hi * wy_lo;
    normal[2] += phi(i+1, j  , k+1) * scaled_dxi[2] * wx_hi * wy_lo;
    normal[2] -= phi(i,   j+1, k  ) * scaled_dxi[2] * wx_lo * wy_hi;
    normal[2] += phi(i  , j+1, k+1) * scaled_dxi[2] * wx_lo * wy_hi;
    normal[2] -= phi(i+1, j+1, k  ) * scaled_dxi[2] * wx_hi * wy_hi;
    normal[2] += phi(i+1, j+1, k+1) * scaled_dxi[2] * wx_hi * wy_hi;
    
    amrex::Real inv_norm = 1.0 / sqrt(dot_product(normal, normal));
    normal[0] *= inv_norm;
    normal[1] *= inv_norm;
    normal[2] *= inv_norm;
}
