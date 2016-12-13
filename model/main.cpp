#include <fstream>
#include <iomanip>

#include <ParmParse.H>
#include <Geometry.H>
#include <VisMF.H>
#include <iMultiFab.H>

#include <mfix_F.H>

int main (int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    Real strt_time = ParallelDescriptor::second();

    mfix_MAIN();

    Real end_time = ParallelDescriptor::second() - strt_time;

    if (ParallelDescriptor::IOProcessor())
       std::cout << "Time spent in main " << end_time << std::endl;

    BoxLib::Finalize();
    return 0;
}
