#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_AmrCore.H>

#include "AMReX_buildInfo.H"

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
    // Define the list of vector variables on faces that need to be written
    // to plotfile/checkfile.
    vecVarsName = {"u_g", "v_g", "w_g"}; 
    vectorVars  = { &u_g, &v_g, &w_g };

    // Define the list of scalar variables at cell centers that need to be written
    // to plotfile/checkfile.
    scaVarsName = {"ep_g", "p_g", "ro_g", "rop_g",  "mu_g"}; 
    scalarVars  = { &ep_g, &p_g, &ro_g,  &rop_g,  &mu_g};
}

void
mfix_level::WriteHeader(const std::string& name, int nstep, Real dt, Real time) const
{
    if (ParallelDescriptor::IOProcessor())
    {
    	std::string HeaderFileName(name + "/Header");
    	std::ofstream HeaderFile(HeaderFileName.c_str(), std::ofstream::out   |
				 std::ofstream::trunc |
				 std::ofstream::binary);

    	if ( ! HeaderFile.good() )
	{
    	    amrex::FileOpenFailed(HeaderFileName);
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
mfix_level::WriteCheckPointFile(std::string& check_file, int nstep, Real dt, Real time ) const 
{
    BL_PROFILE("mfix_level::WriteCheckPointFile()");

    const std::string& checkpointname = amrex::Concatenate( check_file, nstep );

    if (ParallelDescriptor::IOProcessor()) {
    	std::cout << "\n\t Writing checkpoint " << checkpointname << std::endl;
    }

    const int nlevels = finestLevel()+1;
    amrex::PreBuildDirectorHierarchy(checkpointname, level_prefix, nlevels, true);

    WriteHeader(checkpointname, nstep, dt, time);
    
    WriteJobInfo(checkpointname);

    for (int lev = 0; lev < nlevels; ++lev) {

	// Write vector variables
	for (int i = 0; i < vectorVars.size(); i++ ) {
	    VisMF::Write( *((*vectorVars[i])[lev]),
			  amrex::MultiFabFileFullPrefix(lev, checkpointname, 
							level_prefix, vecVarsName[i]));
	}

	// Write scalar variables
	for (int i = 0; i < scalarVars.size(); i++ ) {
	    VisMF::Write( *((*scalarVars[i])[lev]),
			  amrex::MultiFabFileFullPrefix(lev, checkpointname, 
							level_prefix, scaVarsName[i]));
	}

    }

    if ( solve_dem )
    {
        Array<std::string> real_comp_names;
        Array<std::string>  int_comp_names;
        real_comp_names.push_back("radius");
        real_comp_names.push_back("volume");
        real_comp_names.push_back("mass");
        real_comp_names.push_back("density");
        real_comp_names.push_back("omoi");
        real_comp_names.push_back("velx");
        real_comp_names.push_back("vely");
        real_comp_names.push_back("velz");
        real_comp_names.push_back("omegax");
        real_comp_names.push_back("omegay");
        real_comp_names.push_back("omegaz");
        real_comp_names.push_back("dragx");
        real_comp_names.push_back("dragy");
        real_comp_names.push_back("dragz");
         int_comp_names.push_back("phase");
         int_comp_names.push_back("state");

       bool is_checkpoint = true;  
       pc -> Checkpoint(checkpointname, "particles", is_checkpoint, real_comp_names, int_comp_names);
    }
}

void
mfix_level::Restart (std::string& restart_file, int *nstep, Real *dt, Real *time, 
                     IntVect& Nrep) 
{
    BL_PROFILE("mfix_level::Restart()");

    if (ParallelDescriptor::IOProcessor()) 
	std::cout << "  Restarting from checkpoint " << restart_file << std::endl;

    if (ParallelDescriptor::IOProcessor()) 
	std::cout << "  Replication " << Nrep << std::endl;

    Real prob_lo[BL_SPACEDIM];
    Real prob_hi[BL_SPACEDIM];

    // Header
    {
    	std::string File(restart_file + "/Header");

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

       	std::getline(is, line);
       	{
       	    std::istringstream lis(line);
       	    int i = 0;
       	    while (lis >> word) {
       		prob_lo[i++] = std::stod(word);
       	    }
       	}
 
       	std::getline(is, line);
       	{
       	    std::istringstream lis(line);
       	    int i = 0;
       	    while (lis >> word) {
       		prob_hi[i++] = std::stod(word);
       	    }
       	}

     	
        if (Nrep != IntVect::TheUnitVector())
        {
           for (int d = 0; d < BL_SPACEDIM; d++)
           {
              prob_lo[d] = Nrep[d]*prob_lo[d];
              prob_hi[d] = Nrep[d]*prob_hi[d];
           }
        } 
        Geometry::ProbDomain(RealBox(prob_lo,prob_hi));

       	for (int lev = 0; lev < nlevs; ++lev) {
       	    BoxArray orig_ba,ba;
       	    orig_ba.readFrom(is);
       	    GotoNextLine(is);
            Box orig_domain(orig_ba.minimalBox());
            if (ParallelDescriptor::IOProcessor())
               std::cout << " OLD BA HAS " << orig_ba.size() << " GRIDS " << std::endl;

            BoxList bl;
            for (int nb = 0; nb < orig_ba.size(); nb++) {
             for (int k = 0; k < Nrep[2]; k++) {
                 for (int j = 0; j < Nrep[1]; j++) {
                   for (int i = 0; i < Nrep[0]; i++) {
                      Box b(orig_ba[nb]);
                      IntVect lo(b.smallEnd());
                      IntVect hi(b.bigEnd());
                      IntVect shift_vec(i*orig_domain.length(0),
                                        j*orig_domain.length(1),
                                        k*orig_domain.length(2));
                      b.shift(shift_vec);
                      bl.push_back(b);
                   }
                 }
               }
            }
            ba.define(bl);

            if (ParallelDescriptor::IOProcessor())
               std::cout << " NEW BA HAS " << ba.size() << " GRIDS " << std::endl;
            SetBoxArray(lev, ba);
       	    DistributionMapping dm { ba, ParallelDescriptor::NProcs() };
            SetDistributionMap(lev, dm);
            AllocateArrays(lev);
       	}
    }

    // Initialize the field data
    for (int lev = 0, nlevs=finestLevel()+1; lev < nlevs; ++lev)
    {
	// Read vector variables
	for (int i = 0; i < vectorVars.size(); i++ ) {
    	    MultiFab mf;
	    VisMF::Read(mf, amrex::MultiFabFileFullPrefix(lev, restart_file, level_prefix,
  					  vecVarsName[i]));

            if (Nrep == IntVect::TheUnitVector())
            {

   	       // Simply copy mf into vectorVars
   	       (*vectorVars[i])[lev] -> copy(mf, 0, 0, 1, 0, 0);

            } else {

               if (mf.boxArray().size() > 1) 
                   amrex::Abort("Replication only works if one initial grid");

               mf.FillBoundary(geom[lev].periodicity());

               Box edge_bx = mf.boxArray()[0];
               edge_bx.surroundingNodes(i);

               FArrayBox single_fab(edge_bx,1);
   	       mf.copyTo(single_fab);

   	       // Copy and replicate mf into vectorVars
	       for (MFIter mfi( *(*vectorVars[i])[lev] ); mfi.isValid(); ++mfi)
	       {
                   int ib = mfi.index();
	           (*(*vectorVars[i])[lev])[ib].copy(single_fab,single_fab.box(),0,mfi.validbox(),0,1);
	       }
            }
	}

	// Read scalar variables
	for (int i = 0; i < scalarVars.size(); i++ ) {
    	    MultiFab mf;
	    VisMF::Read(mf, amrex::MultiFabFileFullPrefix(lev, restart_file, level_prefix,
  					  scaVarsName[i]));
            if (Nrep == IntVect::TheUnitVector())
            {

   	       // Simply copy mf into scalarVars
	       (*scalarVars[i])[lev] -> copy(mf, 0, 0, 1, 0, 0);

            } else {

               if (mf.boxArray().size() > 1) 
                   amrex::Abort("Replication only works if one initial grid");

               mf.FillBoundary(geom[lev].periodicity());

               FArrayBox single_fab(mf.boxArray()[0],1);
   	       mf.copyTo(single_fab);

   	       // Copy and replicate mf into scalarVars
	       for (MFIter mfi( *(*scalarVars[i])[lev] ); mfi.isValid(); ++mfi)
	       {
                   int ib = mfi.index();
	           (*(*scalarVars[i])[lev])[ib].copy(single_fab,single_fab.box(),0,mfi.validbox(),0,1);
	       }
            }
	}
    }

    // Initialize particles
    pc->Restart(restart_file, "particles");

    // This call to Replicate adds the new particles, then calls BuildLevelMask and Redistribute
    int lev = 0;
    pc->Replicate(Nrep,geom[lev],dmap[lev],grids[lev]);
}

void
mfix_level::GotoNextLine (std::istream& is)
{ 
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}

void mfix_level::WriteJobInfo (const std::string& dir) const
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
	jobInfoFile << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

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
	    jobInfoFile << "AMReX git hash: " << githash2 << "\n";
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


void mfix_level::WritePlotFile (std::string& plot_file, int nstep, Real dt, Real time ) const
{
    BL_PROFILE("mfix_level::WritePlotFile()");

    const std::string& plotfilename = amrex::Concatenate(plot_file,nstep);

    if (ParallelDescriptor::IOProcessor())  
	std::cout << "  Writing plotfile " << plotfilename << std::endl;
    {
	Array< std::unique_ptr<MultiFab> > mf(finest_level+1);
    
	for (int lev = 0; lev <= finest_level; ++lev) {

	    const int ncomp = vectorVars.size() + scalarVars.size();
	    const int ngrow = 0;

	    mf[lev].reset(new MultiFab(grids[lev], dmap[lev], ncomp, ngrow));

	    // Vector variables
	    int dcomp = 0;
	    Array<const MultiFab*> srcmf(3);

	    for( dcomp = 0; dcomp < vectorVars.size(); dcomp=dcomp+3 ) {
		srcmf[0] = (*vectorVars[dcomp])[lev].get();
		srcmf[1] = (*vectorVars[dcomp+1])[lev].get();
		srcmf[2] = (*vectorVars[dcomp+2])[lev].get();
		amrex::average_face_to_cellcenter(*mf[lev], dcomp, srcmf);
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
	names.insert( names.end(), vecVarsName.begin(), vecVarsName.end());
	names.insert( names.end(), scaVarsName.begin(), scaVarsName.end());

	amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf2, names,
				       Geom(), time, istep, refRatio());
    }

    WriteJobInfo(plotfilename);

    if ( solve_dem )
    {
        Array<std::string> real_comp_names;
        Array<std::string>  int_comp_names;
        real_comp_names.push_back("radius");
        real_comp_names.push_back("volume");
        real_comp_names.push_back("mass");
        real_comp_names.push_back("density");
        real_comp_names.push_back("omoi");
        real_comp_names.push_back("velx");
        real_comp_names.push_back("vely");
        real_comp_names.push_back("velz");
        real_comp_names.push_back("omegax");
        real_comp_names.push_back("omegay");
        real_comp_names.push_back("omegaz");
        real_comp_names.push_back("dragx");
        real_comp_names.push_back("dragy");
        real_comp_names.push_back("dragz");
         int_comp_names.push_back("phase");
         int_comp_names.push_back("state");

       bool is_checkpoint = true;
       pc -> Checkpoint(plotfilename, "particles", is_checkpoint, real_comp_names, int_comp_names);
    }
}

void 
mfix_level::WriteParticleAscii ( std::string& par_ascii_file, int nstep ) const 
{
    BL_PROFILE("mfix_level::WriteParticleASCII()");

    const std::string& par_filename = amrex::Concatenate(par_ascii_file,nstep);
    
    pc -> WriteAsciiFile(par_filename);
}

