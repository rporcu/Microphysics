#include <MultiFabUtil.H>
#include <PlotFileUtil.H>

#include <AmrCore.H>

#include "buildInfo.H"

#include "mfix_level.H"


namespace
{
    const std::string level_prefix {"Level_"};
}



void
mfix_level::WriteMfixHeader(const std::string& name) const
{
    if (ParallelDescriptor::IOProcessor())
    {
    	std::string HeaderFileName(name + "/MfixHeader");
    	std::ofstream HeaderFile(HeaderFileName.c_str(), std::ofstream::out   |
				 std::ofstream::trunc |
				 std::ofstream::binary);

    	if ( ! HeaderFile.good() )
	{
    	    BoxLib::FileOpenFailed(HeaderFileName);
    	}

    	HeaderFile.precision(17);

    	VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
    	HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
	
    	HeaderFile << "Checkpoint version: 1\n";

    	const int nlevels = finestLevel()+1;
    	HeaderFile << nlevels << "\n";


    	// Geometry
    	for (int i = 0; i < BL_SPACEDIM; ++i) 
	{
            HeaderFile << Geometry::ProbLo(i) << ' ';
    	}
        HeaderFile << '\n';

        for (int i = 0; i < BL_SPACEDIM; ++i) 
	{
            HeaderFile << Geometry::ProbHi(i) << ' ';
    	}
        HeaderFile << '\n';

    	// BoxArray
    	for (int lev = 0; lev < nlevels; ++lev)
	{
    	    boxArray(lev).writeOn(HeaderFile);
    	    HeaderFile << '\n';
    	}

    }
}




void
mfix_level::WriteCheckPointFile( int nstep, Real dt, Real time )  const
{
    BL_PROFILE("mfix_level::WriteCheckPointFile()");

    // Return if it's not time to dump checkfile yet or 
    // checkfiles have not been enabled ( check_int < 1 )
    if ( (check_int < 1) || ( nstep % check_int != 0 ) )  return;

    const std::string& checkpointname = BoxLib::Concatenate( check_file, nstep );

    if (ParallelDescriptor::IOProcessor()) {
    	std::cout << "  Writing checkpoint " << checkpointname << std::endl;
    }

    const int nlevels = finestLevel()+1;
    BoxLib::PreBuildDirectorHierarchy(checkpointname, level_prefix, nlevels, true);

    WriteMfixHeader(checkpointname);
    
    WriteJobInfo(checkpointname);

    for (int lev = 0; lev < nlevels; ++lev)
    {
    	VisMF::Write(*u_g[0], BoxLib::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ug"));
    	// VisMF::Write(*Efield[lev][1],
    	// 	     BoxLib::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ey"));
    	// VisMF::Write(*Efield[lev][2],
    	// 	     BoxLib::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ez"));
    	// VisMF::Write(*Bfield[lev][0],
    	// 	     BoxLib::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Bx"));
    	// VisMF::Write(*Bfield[lev][1],
    	// 	     BoxLib::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "By"));
    	// VisMF::Write(*Bfield[lev][2],
    	// 	     BoxLib::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Bz"));
    }

}



// void
// mfix_level::InitFromCheckpoint ()
// {
//     BL_PROFILE("mfix_level::InitFromCheckpoint()");

//     if (ParallelDescriptor::IOProcessor()) {
// 	std::cout << "  Restarting from checkpoint " << restart_chkfile << std::endl;
//     }

//     const int checkpoint_nfiles = 64;  // could make this parameter
//     VisMF::SetNOutFiles(checkpoint_nfiles);
    
//     // Header
//     {
// 	std::string File(restart_chkfile + "/MfixHeader");

// 	VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

// 	Array<char> fileCharPtr;
// 	ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
// 	std::string fileCharPtrString(fileCharPtr.dataPtr());
// 	std::istringstream is(fileCharPtrString, std::istringstream::in);

// 	std::string line, word;

// 	std::getline(is, line);

// 	int nlevs;
// 	is >> nlevs;
// 	GotoNextLine(is);
// 	finest_level = nlevs-1;

// 	GotoNextLine(is);

// 	Real prob_lo[BL_SPACEDIM];
// 	std::getline(is, line);
// 	{
// 	    std::istringstream lis(line);
// 	    int i = 0;
// 	    while (lis >> word) {
// 		prob_lo[i++] = std::stod(word);
// 	    }
// 	}
	
// 	Real prob_hi[BL_SPACEDIM];
// 	std::getline(is, line);
// 	{
// 	    std::istringstream lis(line);
// 	    int i = 0;
// 	    while (lis >> word) {
// 		prob_hi[i++] = std::stod(word);
// 	    }
// 	}

// 	Geometry::ProbDomain(RealBox(prob_lo,prob_hi));

// 	for (int lev = 0; lev < nlevs; ++lev) {
// 	    BoxArray ba;
// 	    ba.readFrom(is);
// 	    GotoNextLine(is);
// 	    DistributionMapping dm { ba, ParallelDescriptor::NProcs() };
// 	    MakeNewLevel(lev, ba, dm);
// 	}

// 	// mypc->ReadHeader(is);
//     }

//     // // Initialize the field data
//     // for (int lev = 0, nlevs=finestLevel()+1; lev < nlevs; ++lev)
//     // {
//     // 	for (int i = 0; i < 3; ++i) {
//     // 	    Efield[lev][i]->setVal(0.0);
//     // 	    Bfield[lev][i]->setVal(0.0);
//     // 	    current[lev][i]->setVal(0.0);
//     // 	}

//     // 	// xxxxx This will be done differently in amrex!
//     // 	{
//     // 	    MultiFab mf;
//     // 	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ex"));
//     // 	    Efield[lev][0]->copy(mf, 0, 0, 1, 0, 0);
//     // 	}
//     // 	{
//     // 	    MultiFab mf;
//     // 	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ey"));
//     // 	    Efield[lev][1]->copy(mf, 0, 0, 1, 0, 0);
//     // 	}
//     // 	{
//     // 	    MultiFab mf;
//     // 	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ez"));
//     // 	    Efield[lev][2]->copy(mf, 0, 0, 1, 0, 0);
//     // 	}
//     // 	{
//     // 	    MultiFab mf;
//     // 	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Bx"));
//     // 	    Bfield[lev][0]->copy(mf, 0, 0, 1, 0, 0);
//     // 	}
//     // 	{
//     // 	    MultiFab mf;
//     // 	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "By"));
//     // 	    Bfield[lev][1]->copy(mf, 0, 0, 1, 0, 0);
//     // 	}
//     // 	{
//     // 	    MultiFab mf;
//     // 	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Bz"));
//     // 	    Bfield[lev][2]->copy(mf, 0, 0, 1, 0, 0);
//     // 	}
//     // }

//     // // Initilize particles
//     // mypc->AllocData();
//     // mypc->Restart(restart_chkfile, "particle");
// }




void
mfix_level::WriteJobInfo (const std::string& dir) const
{
    if (ParallelDescriptor::IOProcessor())
    {
	// job_info file with details about the run
	std::ofstream jobInfoFile;
	std::string FullPathJobInfoFile = dir;
	std::string PrettyLine = "===============================================================================\n";

	FullPathJobInfoFile += "/Mfix_job_info";
	jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

	// job information
	jobInfoFile << PrettyLine;
	jobInfoFile << " MFIx Job Information\n";
	jobInfoFile << PrettyLine;

	jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
	jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif

	jobInfoFile << "\n\n";

        // build information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Build Information\n";
	jobInfoFile << PrettyLine;

	jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
	jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
	jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
	jobInfoFile << "BoxLib dir:    " << buildInfoGetBoxlibDir() << "\n";

	jobInfoFile << "\n";

	jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
	jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";
	jobInfoFile << "FCOMP:         " << buildInfoGetFcomp() << "\n";
	jobInfoFile << "FCOMP version: " << buildInfoGetFcompVersion() << "\n";

	jobInfoFile << "\n";

	const char* githash1 = buildInfoGetGitHash(1);
	const char* githash2 = buildInfoGetGitHash(2);
	if (strlen(githash1) > 0) 
	{
	    jobInfoFile << "MFIX git hash: " << githash1 << "\n";
	}
	if (strlen(githash2) > 0)
	{
	    jobInfoFile << "BoxLib git hash: " << githash2 << "\n";
	}

	jobInfoFile << "\n\n";

	// grid information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Grid Information\n";
        jobInfoFile << PrettyLine;

        for (int i = 0; i <= finest_level; i++)
	{
            jobInfoFile << " level: " << i << "\n";
            jobInfoFile << "   number of boxes = " << grids[i].size() << "\n";
            jobInfoFile << "   maximum zones   = ";
            for (int n = 0; n < BL_SPACEDIM; n++)
	    {
                jobInfoFile << geom[i].Domain().length(n) << " ";
	    }
            jobInfoFile << "\n\n";
	}


	// Add here the info on BCs
        // jobInfoFile << " Boundary conditions\n";

        // jobInfoFile << "   -x: " << "interior" << "\n";
        // jobInfoFile << "   +x: " << "interior" << "\n";
        // if (BL_SPACEDIM >= 2) {
	//     jobInfoFile << "   -y: " << "interior" << "\n";
	//     jobInfoFile << "   +y: " << "interior" << "\n";
        // }
        // if (BL_SPACEDIM == 3) {
	//     jobInfoFile << "   -z: " << "interior" << "\n";
	//     jobInfoFile << "   +z: " << "interior" << "\n";
        // }

        jobInfoFile << "\n\n";


	// runtime parameters
	jobInfoFile << PrettyLine;
	jobInfoFile << " Inputs File Parameters\n";
	jobInfoFile << PrettyLine;

	ParmParse::dumpTable(jobInfoFile, true);

	jobInfoFile.close();
    }
}



