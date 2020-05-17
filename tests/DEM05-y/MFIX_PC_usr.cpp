#include <MFIXParticleContainer.H>

void MFIXParticleContainer::usr2_des (int np,
                                      NeighborParticleContainer<16, 2>::ParticleType*& particles)
{
// Purpose: This routine is called within the discrete phase time loop 
// after the source terms are applied and the time step updated. 

    // Move particles (Fortran ordering) 63-93 below particles 32-62 to fake a wall.
    for (unsigned i = 62; i < 94; ++i)
    {
        particles[i].pos(1)   = 0.0475;
        particles[i].pos(0)   = particles[i-31].pos(0);
        particles[i].pos(2)   = particles[i-31].pos(2);

        particles[i].rdata(realData::velx)  = 0.0;
        particles[i].rdata(realData::vely)  = 0.0;
        particles[i].rdata(realData::velz)  = 0.0;

        particles[i].rdata(realData::omegax)  = 0.0;
        particles[i].rdata(realData::omegay)  = 0.0;
        particles[i].rdata(realData::omegaz)  = 0.0;
    }
}

void MFIXParticleContainer::usr3_des (int np, void*&)
{
}
