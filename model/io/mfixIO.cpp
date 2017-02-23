#include <MultiFabUtil.H>
#include <PlotFileUtil.H>

#include <AmrCore.H>

#include "buildInfo.H"

#include "mfix_level.H"


namespace
{
    const std::string level_prefix {"Level_"};
}


// This function iniatializes the attributes vectorVars, scalarVars, vecVarsName and
// scaVarsName.				
// If new variables need to be added to the output/checkpoint, simply
// add them here and the IO routines will automatically take care of them.
void
mfix_level::InitIOData ()
{
    // Define the list of vector variables that need to be written
    // to plotfile/checkfile.
    // These are the variables that need interpolation from face to 
    // node when written to plotfile
    vecVarsName = {"u_g", "v_g", "w_g"}; 
    vectorVars  = { &u_g, &v_g, &w_g };

    // Define the list of scalar variables that need to be written
    // to plotfile/checkfile.
    scaVarsName = {"ep_g", "p_g", "ro_g", "rop_g",  "mu_g"}; 
    scalarVars  = { &ep_g, &p_g, &ro_g,  &rop_g,  &mu_g};
}


void
mfix_level::WriteMfixHeader(const std::string& name, int nstep, Real dt, Real time) const
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

	// Time stepping controls
    	HeaderFile << nstep << "\n";
    	HeaderFile << dt << "\n";
    	HeaderFile << time << "\n";

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
    	std::cout << "\n\t Writing checkpoint " << checkpointname << std::endl;
    }

    const int nlevels = finestLevel()+1;
    BoxLib::PreBuildDirectorHierarchy(checkpointname, level_prefix, nlevels, true);

    WriteMfixHeader(checkpointname, nstep, dt, time);
    
    WriteJobInfo(checkpointname);

    for (int lev = 0; lev < nlevels; ++lev) {

	// Velocities
    	VisMF::Write(*u_g[lev], BoxLib::MultiFabFileFullPrefix(lev, checkpointname, 
    							       level_prefix, "u_g"));
    	VisMF::Write(*u_go[lev], BoxLib::MultiFabFileFullPrefix(lev, checkpointname, 
    								level_prefix, "u_go"));
    	VisMF::Write(*u_gt[lev], BoxLib::MultiFabFileFullPrefix(lev, checkpointname, 
    								level_prefix, "u_gt"));
    	VisMF::Write(*v_g[lev], BoxLib::MultiFabFileFullPrefix(lev, checkpointname, 
    							       level_prefix, "v_g"));
    	VisMF::Write(*v_go[lev], BoxLib::MultiFabFileFullPrefix(lev, checkpointname, 
    								level_prefix, "v_go"));
    	VisMF::Write(*v_gt[lev], BoxLib::MultiFabFileFullPrefix(lev, checkpointname, 
    								level_prefix, "v_gt"));
    	VisMF::Write(*w_g[lev], BoxLib::MultiFabFileFullPrefix(lev, checkpointname, 
    							       level_prefix, "w_g"));
    	VisMF::Write(*w_go[lev], BoxLib::MultiFabFileFullPrefix(lev, checkpointname, 
    								level_prefix, "w_go"));
    	VisMF::Write(*w_gt[lev], BoxLib::MultiFabFileFullPrefix(lev, checkpointname, 
    								level_prefix, "w_gt"));

	// Material properties
    	VisMF::Write(*ep_g[lev], BoxLib::MultiFabFileFullPrefix(lev, checkpointname, 
    							       level_prefix, "ep_g"));
    	VisMF::Write(*ep_go[lev], BoxLib::MultiFabFileFullPrefix(lev, checkpointname, 
    								level_prefix, "ep_go"));
    	VisMF::Write(*p_g[lev], BoxLib::MultiFabFileFullPrefix(lev, checkpointname, 
    								level_prefix, "p_g"));
    	VisMF::Write(*p_go[lev], BoxLib::MultiFabFileFullPrefix(lev, checkpointname, 
    								level_prefix, "p_go"));
    	VisMF::Write(*ro_g[lev], BoxLib::MultiFabFileFullPrefix(lev, checkpointname, 
    								level_prefix, "ro_g"));
    	VisMF::Write(*ro_go[lev], BoxLib::MultiFabFileFullPrefix(lev, checkpointname, 
    								level_prefix, "ro_go"));
    	VisMF::Write(*rop_g[lev], BoxLib::MultiFabFileFullPrefix(lev, checkpointname, 
    								level_prefix, "rop_g"));
    	VisMF::Write(*rop_go[lev], BoxLib::MultiFabFileFullPrefix(lev, checkpointname, 
    								level_prefix, "rop_go"));
    	VisMF::Write(*mu_g[lev], BoxLib::MultiFabFileFullPrefix(lev, checkpointname, 
    								level_prefix, "mu_g"));
    }

    mypc->Checkpoint(checkpointname, "particle", true);
}



int
mfix_level::IsRestartEnabled () const
{
    return ( restart_chkfile.compare("no_restart") != 0 );
}


void
mfix_level::InitFromCheckpoint (int *nstep, Real *dt, Real *time) const
{
    BL_PROFILE("mfix_level::InitFromCheckpoint()");

    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "  Restarting from checkpoint " << restart_chkfile << std::endl;
    }


    // Header
    {
    	std::string File(restart_chkfile + "/MfixHeader");

    	VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    	Array<char> fileCharPtr;
    	ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    	std::string fileCharPtrString(fileCharPtr.dataPtr());
    	std::istringstream is(fileCharPtrString, std::istringstream::in);

    	std::string line, word;

    	std::getline(is, line);

    	int  nlevs;
	int  int_tmp;
	Real real_tmp;

    	is >> nlevs;
    	GotoNextLine(is);
    	// finest_level = nlevs-1;

	// Time stepping controls
	is >> int_tmp;
	*nstep = int_tmp;
    	GotoNextLine(is);

	is >> real_tmp;
	*dt = real_tmp;
    	GotoNextLine(is);

	is >> real_tmp;
	*time = real_tmp;
     	GotoNextLine(is);

//     	Real prob_lo[BL_SPACEDIM];
//     	std::getline(is, line);
//     	{
//     	    std::istringstream lis(line);
//     	    int i = 0;
//     	    while (lis >> word) {
//     		prob_lo[i++] = std::stod(word);
//     	    }
//     	}
	
//     	Real prob_hi[BL_SPACEDIM];
//     	std::getline(is, line);
//     	{
//     	    std::istringstream lis(line);
//     	    int i = 0;
//     	    while (lis >> word) {
//     		prob_hi[i++] = std::stod(word);
//     	    }
//     	}

//     	// Geometry::ProbDomain(RealBox(prob_lo,prob_hi));

//     	// for (int lev = 0; lev < nlevs; ++lev) {
//     	//     BoxArray ba;
//     	//     ba.readFrom(is);
//     	//     GotoNextLine(is);
//     	//     DistributionMapping dm { ba, ParallelDescriptor::NProcs() };
//     	//     MakeNewLevel(lev, ba, dm);
//     	// }

// //	mypc->ReadHeader(is);
    }

    int nghost;
    if (ParallelDescriptor::NProcs() == 1) {
       nghost = 1;
    } else {
       nghost = 2;
    }

    // Initialize the field data
    for (int lev = 0, nlevs=finestLevel()+1; lev < nlevs; ++lev)
    {
	// >>>>>>>>>>>>>>>>>>>>>>>>>>
	// Velocities
	// >>>>>>>>>>>>>>>>>>>>>>>>>>

    	// xxxxx This will be done differently in amrex!
    	{
    	    MultiFab mf;
	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "u_g"));
	    u_g[lev]->copy(mf, 0, 0, 1, 0, 0);
    	}
    	{
    	    MultiFab mf;
    	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "u_go"));
    	    u_go[lev]->copy(mf, 0, 0, 1, 0, 0);
    	}
    	{
    	    MultiFab mf;
    	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "u_gt"));
    	    u_gt[lev]->copy(mf, 0, 0, 1, 0, 0);
    	}
    	{
    	    MultiFab mf;
    	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "v_g"));
    	    v_g[lev]->copy(mf, 0, 0, 1, 0, 0);
    	}
    	{
    	    MultiFab mf;
    	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "v_go"));
    	    v_go[lev]->copy(mf, 0, 0, 1, 0, 0);
    	}
    	{
    	    MultiFab mf;
    	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "v_gt"));
    	    v_gt[lev]->copy(mf, 0, 0, 1, 0, 0);
    	}
    	{
    	    MultiFab mf;
    	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "w_g"));
    	    w_g[lev]->copy(mf, 0, 0, 1, 0, 0);
    	}
    	{
    	    MultiFab mf;
    	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "w_go"));
    	    w_go[lev]->copy(mf, 0, 0, 1, 0, 0);
    	}
    	{
    	    MultiFab mf;
    	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "w_gt"));
    	    w_gt[lev]->copy(mf, 0, 0, 1, 0, 0);
    	}

    	{
    	    MultiFab mf;
    	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "ep_g"));
    	    ep_g[lev]->copy(mf, 0, 0, 1, 0, 0);
    	}
    	{
    	    MultiFab mf;
    	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "ep_go"));
    	    ep_go[lev]->copy(mf, 0, 0, 1, 0, 0);
    	}
    	{
    	    MultiFab mf;
    	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "p_g"));
    	    p_g[lev]->copy(mf, 0, 0, 1, 0, 0);
    	}
    	{
    	    MultiFab mf;
    	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "p_go"));
    	    p_go[lev]->copy(mf, 0, 0, 1, 0, 0);
    	}
    	{
    	    MultiFab mf;
    	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "ro_g"));
    	    ro_g[lev]->copy(mf, 0, 0, 1, 0, 0);
    	}
    	{
    	    MultiFab mf;
    	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "ro_go"));
    	    ro_go[lev]->copy(mf, 0, 0, 1, 0, 0);
    	}
    	{
    	    MultiFab mf;
    	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "rop_g"));
    	    rop_g[lev]->copy(mf, 0, 0, 1, 0, 0);
    	}
    	{
    	    MultiFab mf;
    	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "rop_go"));
    	    rop_go[lev]->copy(mf, 0, 0, 1, 0, 0);
    	}
    	{
    	    MultiFab mf;
    	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "mu_g"));
    	    mu_g[lev]->copy(mf, 0, 0, 1, 0, 0);
    	}

    }

    // Initilize particles
    mypc->Restart(restart_chkfile, "particle");
}


void
mfix_level::GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}

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


void
mfix_level::WritePlotFile ( int nstep, Real dt, Real time ) const
{
    BL_PROFILE("mfix_vele::WritePlotFile()");

    // Return if it's not time to dump plotfile yet or 
    // plotfile have not been enabled ( check_int < 1 )
    // Condition nstep !=0 ensures that plotfile is dumped if WritePlotFile is 
    // called without arguments (useful for steady state case).
    if ( (plot_int < 0) || ( nstep % plot_int != 0 && nstep != 0) )  return;

    const std::string& plotfilename = BoxLib::Concatenate(plot_file,nstep);

    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "  Writing plotfile " << plotfilename << std::endl;
    }
    
    {
	Array< std::unique_ptr<MultiFab> > mf(finest_level+1);
    
	for (int lev = 0; lev <= finest_level; ++lev) {

	    const int ncomp = vectorVars.size() + scalarVars.size();
	    const int ngrow = 0;

	    mf[lev].reset(new MultiFab(grids[lev], ncomp, ngrow, dmap[lev]));

	    // Vector variables
	    int dcomp = 0;
	    std::vector<MultiFab*> srcmf(3);

	    for( dcomp = 0; dcomp < vectorVars.size(); dcomp=dcomp+3 ) {

		srcmf[0] = (*vectorVars[dcomp])[lev].get();
		srcmf[1] = (*vectorVars[dcomp+1])[lev].get();
		srcmf[2] = (*vectorVars[dcomp+2])[lev].get();
		BoxLib::average_face_to_cellcenter(*mf[lev], dcomp, srcmf);

	    };

	    // Scalar variables
	    for( int i = 0; i < scalarVars.size(); i++ ) {

		MultiFab::Copy(*mf[lev], *((*scalarVars[i])[lev].get()), 0, dcomp, 1, 0);
		dcomp++;

	    }		

	}
    
	Array<const MultiFab*> mf2(finest_level+1);

	for (int lev = 0; lev <= finest_level; ++lev) {
	    mf2[lev] = mf[lev].get();
	}

	// Concatenate scalar and vector var names
	Array<std::string>  names;
	names.insert( names.end(),vecVarsName.begin(), vecVarsName.end());
	names.insert( names.end(),scaVarsName.begin(), scaVarsName.end());

	BoxLib::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf2, names,
					Geom(), time, istep, refRatio());
    }

    mypc->Checkpoint(plotfilename, "particle", false);

    WriteJobInfo(plotfilename);

    WriteMfixHeader(plotfilename, nstep, dt, time);
}



