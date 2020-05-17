#include <MFIXParticleContainer.H>

//  This routine is called before the discrete phase time loop
void MFIXParticleContainer::usr0_des ()
{
}

//  This routine is called within the discrete phase time loop 
//  after the source terms have been calculated but before they are applied.
void MFIXParticleContainer::usr1_des ()
{
}

void MFIXParticleContainer::usr2_des (int np,
                                      NeighborParticleContainer<16, 2>::ParticleType*&)
{
}

void MFIXParticleContainer::usr3_des (int np, void*&)
{
}
